# -*- coding: utf-8 -*-
"""
natural_drainage.py
--------------------
Modulo OPCIONAL de red de drenaje natural. Es independiente del flujo
principal de calles/rios: solo se invoca si el usuario carga la capa de
drenaje natural.

Los cauces naturales NO se modelan como conductos ni las cuencas naturales
como subcuencas SWMM. Cada outlet (extremo aguas abajo) que cae cerca de una
calle aporta al modelo unicamente como un INFLOW (Direct Inflow, FLOW) en el
nodo de calle vinculado. La serie temporal de ese INFLOW es el hidrograma
unitario sintetico triangular del SCS de la cuenca de aporte (delimitada con
el DEM) convolucionado con la forma del hietograma efectivo:

    Q[n] = sum_m Pe[m] * U[n-m]          (Pe en mm por paso, U en m3/s/mm)
    serie = Q / Pe_total   (m3/s por mm de lluvia efectiva total)
    Scale Factor del INFLOW = Pe_total (mm)   ->  serie * Scale = Q

Funciones publicas:
    load_cn_table()                -> tabla clase de cobertura -> CN
    prepare_drainage_segments()    -> limpia/noda/orienta los ejes, detecta
                                      outlets, conecta a calles y PURGA lo que
                                      no desemboca en ninguna calle
    parse_hyetograph()             -> hietograma (volumen por paso) desde la
                                      serie de un pluviometro
    delineate_natural_basins()     -> cuenca de aporte (D8 sobre DEM), L, H, Tc,
                                      CN, lluvia efectiva y serie de caudal
    assemble_node_inflows()        -> agrupa por nodo y arma [INFLOWS] y las
                                      series para [TIMESERIES]
    tc_california_culverts(), scs_triangular_uh(), scs_effective_rain()

Solo usa la API de PyQGIS y numpy (sin GRASS/SAGA/GDAL para hidrologia).
"""

import csv
import heapq
import math
from datetime import datetime, timedelta

import numpy as np

from qgis.core import (
    QgsCoordinateTransform,
    QgsFeature,
    QgsGeometry,
    QgsPointXY,
    QgsRectangle,
    QgsSpatialIndex,
    QgsWkbTypes,
)

# ----------------------------------------------------------------------
# Constantes
# ----------------------------------------------------------------------
BURN_DEPTH_M = 1.0          # cuanto se baja el DEM bajo los ejes de drenaje (burn)
SNAP_RADIUS_CELLS = 2       # radio (celdas) para ajustar el outlet a la celda de mayor acumulacion
FILL_EPS = 1e-5             # gradiente minimo al rellenar depresiones (m)
DEFAULT_CN = 75.0           # CN si la cobertura de la cuenca no esta en la tabla
MAX_WINDOW_CELLS = 40_000_000
SCS_QP_COEF = 0.208         # Qp[m3/s/mm] = 0.208 * A[km2] / Tp[h]
SCS_BASE_FACTOR = 2.67      # Tb = 2.67 * Tp
SCS_LAG_FACTOR = 0.6        # t_lag = 0.6 * Tc
M_PER_MILE = 1609.344
M_PER_FOOT = 0.3048

# vecinos D8: (dfila, dcolumna)
_NBR = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]

_NP_BY_NAME = {
    "Byte": np.uint8, "Int8": np.int8, "UInt16": np.uint16, "Int16": np.int16,
    "UInt32": np.uint32, "Int32": np.int32, "Float32": np.float32, "Float64": np.float64,
}


# ----------------------------------------------------------------------
# Tabla CN
# ----------------------------------------------------------------------
def load_cn_table(path):
    """CSV de cobertura (landuse_imperviousness.csv) con columnas class_value y cn.
    Devuelve {clase(int): CN(float)}; vacio si la columna 'cn' no existe."""
    table = {}
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or "cn" not in [n.strip().lower() for n in reader.fieldnames]:
            return {}
        for row in reader:
            row = {str(k).strip().lower(): v for k, v in row.items() if k}
            try:
                table[int(round(float(row["class_value"])))] = float(row["cn"])
            except (KeyError, ValueError, TypeError):
                continue
    return table


# ----------------------------------------------------------------------
# Helpers geometricos
# ----------------------------------------------------------------------
def union_geometries(layer_or_source):
    """Union de todas las geometrias de una capa (o None si esta vacia)."""
    geoms = []
    for feat in layer_or_source.getFeatures():
        g = feat.geometry()
        if g is None or g.isEmpty():
            continue
        if not g.isGeosValid():
            g = g.makeValid()
        geoms.append(g)
    if not geoms:
        return None
    return QgsGeometry.unaryUnion(geoms)


def _line_parts(geom):
    """Descompone una geometria en LineStrings simples (descarta puntos/poligonos)."""
    if geom is None or geom.isEmpty():
        return []
    out = []
    for g in geom.asGeometryCollection():
        if g.type() == QgsWkbTypes.LineGeometry and not g.isEmpty():
            out.append(g)
    return out


def _endpoints(geom):
    pl = geom.asPolyline()
    return pl[0], pl[-1]


def _points_of(geom):
    """Puntos contenidos en el resultado de una interseccion (Point/MultiPoint/Collection)."""
    pts = []
    if geom is None or geom.isEmpty():
        return pts
    for g in geom.asGeometryCollection():
        if g.type() == QgsWkbTypes.PointGeometry:
            if g.isMultipart():
                pts.extend(g.asMultiPoint())
            else:
                pts.append(g.asPoint())
    return pts


def _split_at_distances(geom, dists, tol):
    """Corta una linea en las distancias dadas (a lo largo de la linea)."""
    total = geom.length()
    ds = sorted(d for d in dists if tol < d < total - tol)
    clean = []
    for d in ds:
        if not clean or d - clean[-1] > tol:
            clean.append(d)
    if not clean:
        return [geom]
    bounds = [0.0] + clean + [total]
    curve = geom.constGet()
    pieces = []
    for a, b in zip(bounds[:-1], bounds[1:]):
        if b - a <= 0:
            continue
        pieces.append(QgsGeometry(curve.curveSubstring(a, b)))
    return pieces


def _node_lines(pieces, tol):
    """Parte cada linea en los cruces con otras y donde el extremo de otra
    linea la toca (confluencias en T). Devuelve lista de (geom, etiqueta)."""
    idx = QgsSpatialIndex()
    for i, (g, _t) in enumerate(pieces):
        f = QgsFeature()
        f.setId(i)
        f.setGeometry(g)
        idx.addFeature(f)

    out = []
    for i, (a, tname) in enumerate(pieces):
        box = a.boundingBox()
        box.grow(tol)
        dists = []
        for j in idx.intersects(box):
            if j == i:
                continue
            b = pieces[j][0]
            for p in _endpoints(b):
                pg = QgsGeometry.fromPointXY(p)
                if a.distance(pg) <= tol:
                    dists.append(a.lineLocatePoint(pg))
            inter = a.intersection(b)
            for p in _points_of(inter):
                dists.append(a.lineLocatePoint(QgsGeometry.fromPointXY(p)))
        for piece in _split_at_distances(a, dists, tol):
            if piece.length() >= tol:
                out.append((piece, tname))
    return out


def _reverse(geom):
    pl = geom.asPolyline()
    pl.reverse()
    return QgsGeometry.fromPolylineXY(pl)


def _street_end_points(layer):
    pts = []
    for feat in layer.getFeatures():
        g = feat.geometry()
        if g is None or g.isEmpty():
            continue
        pl = g.asPolyline() if not g.isMultipart() else g.asMultiPolyline()[0]
        if len(pl) >= 2:
            pts.append(QgsPointXY(pl[0]))
            pts.append(QgsPointXY(pl[-1]))
    return pts


def _point_index(points):
    idx = QgsSpatialIndex()
    for i, p in enumerate(points):
        f = QgsFeature()
        f.setId(i)
        f.setGeometry(QgsGeometry.fromPointXY(p))
        idx.addFeature(f)
    return idx


def _box(p, d):
    return QgsRectangle(p.x() - d, p.y() - d, p.x() + d, p.y() + d)


class SinkRef:
    """Referencia de SALIDAS del modelo (rios y storages recortados al AOI).
    distance(punto, maxd) -> distancia minima al sumidero mas cercano si es
    <= maxd; si no hay ninguno dentro de maxd devuelve None."""

    def __init__(self, geoms):
        self._geoms = {}
        self._idx = QgsSpatialIndex()
        k = 0
        for g in geoms:
            if g is None or g.isEmpty():
                continue
            k += 1
            f = QgsFeature()
            f.setId(k)
            f.setGeometry(g)
            self._idx.insertFeature(f)
            self._geoms[k] = g

    def __len__(self):
        return len(self._geoms)

    def distance(self, pt, maxd):
        if not self._geoms:
            return None
        p = QgsPointXY(pt)
        pg = QgsGeometry.fromPointXY(p)
        best = None
        for i in self._idx.intersects(_box(p, maxd)):
            d = self._geoms[i].distance(pg)
            if d <= maxd and (best is None or d < best):
                best = d
        return best


# ----------------------------------------------------------------------
# 1) Segmentos de drenaje, outlets, conexion a calles y PURGA
# ----------------------------------------------------------------------
def prepare_drainage_segments(drain_source, target_crs, transform_context, aoi_geom,
                              blocks_union, streets_layer, sink_ref, elev_fn,
                              tol, connect_buffer, feedback):
    """
    Devuelve un dict:
        "segments": lista de tramos CONSERVADOS
                    {"geom" (linea aguas arriba -> aguas abajo), "seg_id", "length"}
        "purged":   lista de tramos ELIMINADOS (no llegan a ninguna calle)
        "outlets":  lista de outlets conectados
                    {"point": QgsPointXY (extremo aguas abajo),
                     "target": QgsPointXY (extremo de calle/rio al que se vincula),
                     "seg_ids": [seg_id, ...], "link": QgsGeometry | None}
        "to_sink":  tramos que drenan a un rio/storage (descartados, sin aporte)
    La red se recorta al AOI y se le descuentan las manzanas (area urbana ya
    controlada por el drenaje de calles). El outlet es el extremo aguas abajo
    que ningun otro tramo continua. Los tramos intermedios NO se cortan.
    Se purga todo tramo cuyo recorrido aguas abajo no termina en un outlet
    conectado a una calle dentro de `connect_buffer`.
    Si el outlet esta a <= `connect_buffer` de un rio/storage (`sink_ref`) y mas
    cerca que el extremo de calle mas cercano, desemboca en una SALIDA del
    modelo: se descarta (sin INFLOW) y sus tramos aguas arriba se devuelven en
    "to_sink" (se dibujan 'a_salida'). Devuelve ademas "n_sink_outlets".
    """
    to_target = None
    if drain_source.sourceCrs() != target_crs:
        to_target = QgsCoordinateTransform(drain_source.sourceCrs(), target_crs, transform_context)

    # --- 1. Leer, recortar al AOI, descontar manzanas, explotar multiparte ---
    pieces = []
    n_outside = 0
    for feat in drain_source.getFeatures():
        g = feat.geometry()
        if g is None or g.isEmpty():
            continue
        if to_target is not None:
            g.transform(to_target)
        if not g.isGeosValid():
            g = g.makeValid()
        if not g.intersects(aoi_geom):
            n_outside += 1
            continue
        g = g.intersection(aoi_geom)
        if blocks_union is not None:
            g = g.difference(blocks_union)
        for part in _line_parts(g):
            if part.length() >= tol:
                pieces.append((part, None))

    if n_outside:
        feedback.pushInfo(f"  {n_outside} elemento(s) de la capa de drenaje quedan totalmente fuera del AOI y se ignoran.")
    empty = {"segments": [], "purged": [], "outlets": [], "unlinked": [], "to_sink": []}
    if not pieces:
        feedback.pushWarning(
            "La red de drenaje natural no tiene tramos dentro del area de interes "
            "(o quedaron todos dentro de manzanas). Se omite."
        )
        return empty

    # --- 2. Nodar (cruces + confluencias en T) ---
    pieces = _node_lines(pieces, tol)

    # --- 3. Orientar segun DEM (aguas arriba -> aguas abajo) ---
    segments = []
    for seg_id, (g, _t) in enumerate(pieces, start=1):
        p0, p1 = _endpoints(g)
        e0, e1 = elev_fn(p0.x(), p0.y()), elev_fn(p1.x(), p1.y())
        if e0 is not None and e1 is not None and e0 < e1:
            g = _reverse(g)
        segments.append({"geom": g, "seg_id": seg_id, "length": g.length()})

    starts = [_endpoints(s["geom"])[0] for s in segments]
    ends = [_endpoints(s["geom"])[1] for s in segments]
    starts_idx = _point_index(starts)
    ends_idx = _point_index(ends)

    street_pts = _street_end_points(streets_layer)
    street_idx = _point_index(street_pts) if street_pts else None

    def near_any(idx, pts, p, dist):
        if idx is None:
            return False
        for i in idx.intersects(_box(p, dist)):
            if pts[i].distance(p) <= dist:
                return True
        return False

    def upstream_set(group):
        seen, stack = set(), list(group)
        while stack:
            k = stack.pop()
            if k in seen:
                continue
            seen.add(k)
            st = starts[k]
            for j in ends_idx.intersects(_box(st, tol)):
                if j != k and j not in seen and ends[j].distance(st) <= tol:
                    stack.append(j)
        return seen

    # --- 4. Outlets (extremo aguas abajo que ningun tramo continua) y vinculo ---
    outlets = []
    connected_seed = []          # indices de tramos que terminan en un outlet conectado
    used_pts = []
    n_unlinked = 0
    unlinked = []
    sink_seeds = set()           # tramos que terminan en una salida (rio/storage)
    n_sink = 0
    for k in range(len(segments)):
        outlet = ends[k]
        continues = any(
            i != k and starts[i].distance(outlet) <= tol
            for i in starts_idx.intersects(_box(outlet, tol)))
        if continues:
            continue
        if any(u.distance(outlet) <= tol for u in used_pts):
            continue   # varios tramos terminan en el mismo punto: un solo outlet
        used_pts.append(outlet)
        group = [i for i in ends_idx.intersects(_box(outlet, tol))
                 if ends[i].distance(outlet) <= tol]

        # salida (rio/storage) a <= connect_buffer y mas cerca que la calle mas cercana
        d_sink = sink_ref.distance(outlet, connect_buffer) if sink_ref is not None else None
        if d_sink is not None:
            d_street = None
            if street_idx is not None:
                nn0 = street_idx.nearestNeighbor(outlet, 1)
                if nn0:
                    d_street = street_pts[nn0[0]].distance(outlet)
            if d_street is None or d_sink < d_street:
                n_sink += 1
                sink_seeds.update(group)
                continue

        target = None
        if near_any(street_idx, street_pts, outlet, tol):
            target = QgsPointXY(outlet)          # ya coincide con un nodo de calle
        elif street_idx is not None:
            nn = street_idx.nearestNeighbor(outlet, 1)
            if nn and street_pts[nn[0]].distance(outlet) <= connect_buffer:
                target = QgsPointXY(street_pts[nn[0]])
        if target is None:
            n_unlinked += 1
            near = None
            if street_idx is not None:
                nn2 = street_idx.nearestNeighbor(outlet, 1)
                if nn2:
                    near = (QgsPointXY(street_pts[nn2[0]]), street_pts[nn2[0]].distance(outlet))
            unlinked.append({"point": QgsPointXY(outlet),
                             "nearest": near[0] if near else None,
                             "dist": near[1] if near else None})
            continue
        link = None
        if target.distance(outlet) > tol:
            link = QgsGeometry.fromPolylineXY([outlet, target])
        outlets.append({
            "point": QgsPointXY(outlet), "target": target,
            "seg_ids": [segments[i]["seg_id"] for i in group], "link": link,
            "_group": group,
        })
        connected_seed.extend(group)

    # --- 5. PURGA: conservar solo lo que drena (aguas abajo) a un outlet conectado ---
    kept = set()
    for og in outlets:
        up = upstream_set(og.pop("_group"))
        og["net_seg_ids"] = [segments[i]["seg_id"] for i in up]   # toda la red que drena a este outlet
        kept |= up

    to_sink_set = set()
    for k0 in sink_seeds:
        to_sink_set |= upstream_set([k0])
    to_sink_set -= kept
    kept_segments = [s for i, s in enumerate(segments) if i in kept]
    to_sink = [s for i, s in enumerate(segments) if i in to_sink_set]
    purged = [s for i, s in enumerate(segments) if i not in kept and i not in to_sink_set]

    feedback.pushInfo(
        f"  Drenaje natural: {len(segments)} tramos; {len(outlets)} outlet(s) conectado(s) "
        f"a calles (buffer {connect_buffer} m); {n_sink} outlet(s) caen a salida (rio/storage, "
        f"sin aporte; {len(to_sink)} tramo(s) 'a_salida'); {n_unlinked} outlet(s) sin calle cercana; "
        f"{len(purged)} tramo(s) purgado(s), {len(kept_segments)} conservado(s)."
    )
    ordered = sorted((u for u in unlinked if u["dist"] is not None), key=lambda u: u["dist"])
    if ordered:
        feedback.pushInfo(
            "  Outlets sin conexion (mas cercanos a una calle primero): "
            + "; ".join(f"({u['point'].x():.0f}, {u['point'].y():.0f}) a {u['dist']:.0f} m" for u in ordered[:15])
            + (" ..." if len(ordered) > 15 else "")
            + ". Se dibujan como 'sin_conexion' en la capa de red.")
    return {"segments": kept_segments, "purged": purged, "outlets": outlets, "unlinked": unlinked,
            "to_sink": to_sink, "n_sink_outlets": n_sink}


# ----------------------------------------------------------------------
# 2) Hietograma -> volumen por paso (mm)
# ----------------------------------------------------------------------
def _clock_hours(s):
    s = str(s).strip()
    if ":" in s:
        parts = [float(p) for p in s.split(":")]
        h = parts[0]
        m = parts[1] if len(parts) > 1 else 0.0
        sec = parts[2] if len(parts) > 2 else 0.0
        return h + m / 60.0 + sec / 3600.0
    return float(s)


def _parse_date(s):
    s = str(s).strip()
    for fmt in ("%m/%d/%Y", "%Y-%m-%d"):
        try:
            return datetime.strptime(s, fmt)
        except ValueError:
            continue
    raise ValueError(f"Fecha no reconocida: {s}")


def parse_hyetograph(values, rain_format):
    """values: [(fecha|'', hora, valor), ...] (como en `timeseries` del algoritmo).
    rain_format: 'VOLUME' (mm por paso), 'INTENSITY' (mm/h) o 'CUMULATIVE' (mm).
    Devuelve dict {dt_h, t0_h, dated, step_mm (np.array), total_mm, irregular} o None."""
    if values is None or len(values) < 2:
        return None
    epoch = datetime(1970, 1, 1)
    dated = all(bool(v[0]) for v in values)
    t_h = None
    if dated:
        try:
            t_h = [(_parse_date(v[0]) - epoch).total_seconds() / 3600.0 + _clock_hours(v[1])
                   for v in values]
        except ValueError:
            dated = False
    if not dated:
        t_h = [_clock_hours(v[1]) for v in values]
    vals = np.array([float(v[2]) for v in values], dtype=float)
    t_h = np.array(t_h, dtype=float)
    order = np.argsort(t_h, kind="stable")
    t_h, vals = t_h[order], vals[order]
    diffs = np.diff(t_h)
    diffs = diffs[diffs > 1e-9]
    if diffs.size == 0:
        return None
    dt_h = float(np.median(diffs))
    irregular = bool(np.max(np.abs(diffs - dt_h)) > 0.01 * dt_h)

    fmt = str(rain_format).upper()
    if fmt == "INTENSITY":
        vol = vals * dt_h
    elif fmt == "CUMULATIVE":
        vol = np.clip(np.diff(vals, prepend=0.0), 0.0, None)
    else:  # VOLUME
        vol = vals
    idx = np.rint((t_h - t_h[0]) / dt_h).astype(int)
    step = np.zeros(int(idx.max()) + 1)
    np.add.at(step, idx, vol)
    return {"dt_h": dt_h, "t0_h": float(t_h[0]), "dated": dated,
            "step_mm": step, "total_mm": float(step.sum()), "irregular": irregular}


# ----------------------------------------------------------------------
# 3) Formulas hidrologicas
# ----------------------------------------------------------------------
def tc_california_culverts(L_m, H_m):
    """Tc [min] = 60 * (11.9 * L^3 / H)^0.385 con L en MILLAS y H en PIES.
    Entradas en metros; la conversion a unidades imperiales se hace aca."""
    L_mi = L_m / M_PER_MILE
    H_ft = H_m / M_PER_FOOT
    if L_mi <= 0.0 or H_ft <= 0.0:
        raise ValueError("L y H deben ser positivos para calcular Tc.")
    return 60.0 * ((11.9 * L_mi ** 3) / H_ft) ** 0.385


def scs_triangular_uh(area_km2, tc_min, dt_min):
    """HU triangular SCS para lluvia efectiva de 1 mm en `dt_min`.
    t_lag = 0.6 Tc; Tp = dt/2 + t_lag; Tb = 2.67 Tp; Qp = 0.208 A / Tp [m3/s/mm].
    Ordenadas cada `dt_min`, normalizadas para conservar el volumen (1000*A m3/mm).
    Devuelve (U [m3/s/mm], Qp, Tp_min, Tb_min)."""
    tlag = SCS_LAG_FACTOR * tc_min
    tp = dt_min / 2.0 + tlag
    tb = SCS_BASE_FACTOR * tp
    qp = SCS_QP_COEF * area_km2 / (tp / 60.0)
    n = int(math.ceil(tb / dt_min)) + 1
    t = np.arange(n + 1) * dt_min
    u = np.where(t <= tp, qp * t / tp,
                 np.where(t < tb, qp * (tb - t) / (tb - tp), 0.0))
    vol = float(u.sum()) * dt_min * 60.0
    if vol > 0.0:
        u = u * (1000.0 * area_km2 / vol)
    return u, qp, tp, tb


def scs_effective_rain(step_mm, cn, ia_ratio=0.2):
    """Lluvia efectiva incremental (mm por paso) por SCS-CN sobre la lluvia
    acumulada. Ia = ia_ratio * S (0.2 estandar), S = 25400/CN - 254 (mm)."""
    cn = min(max(float(cn), 1.0), 100.0)
    S = 25400.0 / cn - 254.0
    ia = ia_ratio * S
    cum = np.cumsum(step_mm)
    if S <= 0.0:
        q = cum
    else:
        excess = np.clip(cum - ia, 0.0, None)
        q = excess ** 2 / (excess + S)
    return np.diff(q, prepend=0.0)


# ----------------------------------------------------------------------
# 4) Hidrologia raster (numpy): relleno, D8, acumulacion, cuenca
# ----------------------------------------------------------------------
def _step_lengths(dx, dy):
    diag = math.hypot(dx, dy)
    out = []
    for dr, dc in _NBR:
        if dr != 0 and dc != 0:
            out.append(diag)
        elif dr != 0:
            out.append(dy)
        else:
            out.append(dx)
    return np.array(out, dtype=float)


def _rings_mask(rings, x0, y0, dx, dy, nrows, ncols):
    """Mascara (nrows x ncols) de celdas cuyo centro cae dentro de los anillos
    (regla par-impar). (x0, y0) = esquina superior izquierda de la ventana."""
    cross = np.zeros((nrows, ncols + 1), dtype=np.int32)
    for ring in rings:
        pts = np.asarray(ring, dtype=float)
        if len(pts) < 3:
            continue
        xa, ya = pts[:, 0], pts[:, 1]
        xb, yb = np.roll(xa, -1), np.roll(ya, -1)
        ymin, ymax = np.minimum(ya, yb), np.maximum(ya, yb)
        r_lo = np.floor((y0 - ymax) / dy - 0.5).astype(np.int64) + 1
        r_hi = np.floor((y0 - ymin) / dy - 0.5).astype(np.int64)
        r_lo = np.clip(r_lo, 0, nrows)
        r_hi = np.clip(r_hi, -1, nrows - 1)
        n_e = np.where((ymax > ymin), np.maximum(r_hi - r_lo + 1, 0), 0)
        total = int(n_e.sum())
        if total == 0:
            continue
        e_idx = np.repeat(np.arange(len(pts)), n_e)
        first = np.repeat(np.cumsum(n_e) - n_e, n_e)
        rows = np.repeat(r_lo, n_e) + (np.arange(total) - first)
        yc = y0 - (rows + 0.5) * dy
        x1, y1, x2, y2 = xa[e_idx], ya[e_idx], xb[e_idx], yb[e_idx]
        xint = x1 + (yc - y1) * (x2 - x1) / (y2 - y1)
        cols = np.clip(np.ceil((xint - x0) / dx - 0.5).astype(np.int64), 0, ncols)
        np.add.at(cross, (rows, cols), 1)
    return (np.cumsum(cross[:, :ncols], axis=1) % 2) == 1


def polygon_mask(geom, x0, y0, dx, dy, nrows, ncols):
    """Mascara raster de una geometria poligonal (multipartes y huecos incluidos)."""
    rings = []
    if geom is not None and not geom.isEmpty():
        for g in geom.asGeometryCollection():
            if g.type() != QgsWkbTypes.PolygonGeometry:
                continue
            polys = g.asMultiPolygon() if g.isMultipart() else [g.asPolygon()]
            for poly in polys:
                for ring in poly:
                    rings.append([(p.x(), p.y()) for p in ring])
    if not rings:
        return np.zeros((nrows, ncols), dtype=bool)
    return _rings_mask(rings, x0, y0, dx, dy, nrows, ncols)


def _fill_depressions(zp, vp):
    """Priority-flood con gradiente epsilon (resuelve tambien las zonas planas).
    zp/vp: arrays con borde de 1 celda invalido. Las celdas validas vecinas de
    una celda invalida (fuera del AOI, manzana, sin dato) son salidas del dominio."""
    nr, nc = zp.shape
    W = nc
    zl = np.where(vp, zp, 0.0).ravel().tolist()
    vl = vp.ravel().tolist()
    vis = bytearray(nr * nc)
    inner = vp[1:-1, 1:-1]
    all_nb = np.ones_like(inner)
    for dr, dc in _NBR:
        all_nb &= vp[1 + dr:nr - 1 + dr, 1 + dc:nc - 1 + dc]
    seeds = np.argwhere(inner & ~all_nb)
    flat = ((seeds[:, 0] + 1) * W + (seeds[:, 1] + 1)).tolist()
    heap = [(zl[i], i) for i in flat]
    for _z, i in heap:
        vis[i] = 1
    heapq.heapify(heap)
    offs = [dr * W + dc for dr, dc in _NBR]
    pop, push = heapq.heappop, heapq.heappush
    while heap:
        zc, i = pop(heap)
        for o in offs:
            j = i + o
            if vl[j] and not vis[j]:
                vis[j] = 1
                zn = zl[j]
                if zn <= zc:
                    zn = zc + FILL_EPS
                    zl[j] = zn
                push(heap, (zn, j))
    out = np.array(zl, dtype=float).reshape(nr, nc)
    out[~vp] = np.nan
    return out


def _d8(zp, vp, dx, dy):
    """Direccion de flujo D8 (indice 0..7 de _NBR, -1 = sin receptor / sale del dominio)."""
    nr, nc = zp.shape
    inner_z = zp[1:-1, 1:-1]
    steps = _step_lengths(dx, dy)
    best = np.zeros(inner_z.shape)
    dirs = np.full(inner_z.shape, -1, dtype=np.int8)
    with np.errstate(invalid="ignore"):
        for d, (dr, dc) in enumerate(_NBR):
            nb = zp[1 + dr:nr - 1 + dr, 1 + dc:nc - 1 + dc]
            slope = (inner_z - nb) / steps[d]
            better = slope > best
            best = np.where(better, slope, best)
            dirs = np.where(better, np.int8(d), dirs).astype(np.int8)
    out = np.full(zp.shape, -1, dtype=np.int8)
    out[1:-1, 1:-1] = dirs
    out[~vp] = -1
    return out


def _receivers(dirs, W):
    flat = dirs.ravel()
    recv = np.full(flat.size, -1, dtype=np.int64)
    offs = np.array([dr * W + dc for dr, dc in _NBR], dtype=np.int64)
    idx = np.flatnonzero(flat >= 0)
    recv[idx] = idx + offs[flat[idx]]
    return recv


def _flow_accumulation(recv, valid_flat):
    """Acumulacion de flujo (numero de celdas) por olas (algoritmo de Kahn vectorizado)."""
    N = recv.size
    acc = valid_flat.astype(float)
    m = recv >= 0
    indeg = np.bincount(recv[m], minlength=N)
    frontier = np.flatnonzero(valid_flat & (indeg == 0))
    while frontier.size:
        r = recv[frontier]
        keep = r >= 0
        src, dst = frontier[keep], r[keep]
        if dst.size == 0:
            break
        np.add.at(acc, dst, acc[src])
        u, c = np.unique(dst, return_counts=True)
        indeg[u] -= c
        frontier = u[indeg[u] == 0]
    return acc


def _donor_csr(recv):
    N = recv.size
    donors = np.flatnonzero(recv >= 0)
    tgt = recv[donors]
    order = np.argsort(tgt, kind="stable")
    donors_sorted = donors[order]
    counts = np.bincount(tgt, minlength=N)
    starts = np.cumsum(counts) - counts
    return donors_sorted, starts, counts


def _basin_cells(cell, donors_sorted, starts, counts, dirs_flat, step_len):
    """Celdas aguas arriba de `cell` y su distancia de flujo (m) hasta ella."""
    frontier = np.array([cell], dtype=np.int64)
    fdist = np.array([0.0])
    members, mdist = [frontier], [fdist]
    while True:
        cnt = counts[frontier]
        tot = int(cnt.sum())
        if tot == 0:
            break
        off = np.repeat(starts[frontier] - (np.cumsum(cnt) - cnt), cnt) + np.arange(tot)
        nxt = donors_sorted[off]
        nd_ = np.repeat(fdist, cnt) + step_len[dirs_flat[nxt]]
        members.append(nxt)
        mdist.append(nd_)
        frontier, fdist = nxt, nd_
    return np.concatenate(members), np.concatenate(mdist)


def _cells_to_polygon(rows, cols, x0, y0, dx, dy):
    """Poligono (union de rectangulos por tramos horizontales) de un conjunto de celdas."""
    order = np.lexsort((cols, rows))
    r, c = rows[order], cols[order]
    brk = np.ones(len(r), dtype=bool)
    brk[1:] = (r[1:] != r[:-1]) | (c[1:] != c[:-1] + 1)
    s_idx = np.flatnonzero(brk)
    e_idx = np.append(s_idx[1:], len(r)) - 1
    rects = []
    for s, e in zip(s_idx.tolist(), e_idx.tolist()):
        xa = x0 + c[s] * dx
        xb = x0 + (c[e] + 1) * dx
        yt = y0 - r[s] * dy
        rects.append(QgsGeometry.fromRect(QgsRectangle(xa, yt - dy, xb, yt)))
    return QgsGeometry.unaryUnion(rects)


# ----------------------------------------------------------------------
# 5) Lectura de rasters
# ----------------------------------------------------------------------
def _dem_window(dem_layer, aoi_geom):
    """Ventana del DEM alineada a su grilla que cubre el bbox del AOI."""
    ext = dem_layer.extent()
    dx = dem_layer.rasterUnitsPerPixelX()
    dy = dem_layer.rasterUnitsPerPixelY()
    bb = aoi_geom.boundingBox().intersect(ext)
    if bb.isEmpty() or bb.width() <= 0 or bb.height() <= 0:
        return None
    W, H = dem_layer.width(), dem_layer.height()
    c0 = max(int(math.floor((bb.xMinimum() - ext.xMinimum()) / dx + 1e-9)), 0)
    c1 = min(int(math.ceil((bb.xMaximum() - ext.xMinimum()) / dx - 1e-9)), W)
    r0 = max(int(math.floor((ext.yMaximum() - bb.yMaximum()) / dy + 1e-9)), 0)
    r1 = min(int(math.ceil((ext.yMaximum() - bb.yMinimum()) / dy - 1e-9)), H)
    ncols, nrows = c1 - c0, r1 - r0
    if ncols <= 0 or nrows <= 0:
        return None
    x0 = ext.xMinimum() + c0 * dx
    y0 = ext.yMaximum() - r0 * dy
    return {"x0": x0, "y0": y0, "dx": dx, "dy": dy, "ncols": ncols, "nrows": nrows,
            "rect": QgsRectangle(x0, y0 - nrows * dy, x0 + ncols * dx, y0)}


def _read_raster(provider, win, band=1):
    """Lee la ventana como array float64 (nan = sin dato). None si no se puede."""
    block = provider.block(band, win["rect"], win["ncols"], win["nrows"])
    if block is None or not block.isValid():
        return None
    dt = block.dataType()
    name = getattr(dt, "name", None) or str(dt).split(".")[-1]
    npdt = _NP_BY_NAME.get(name)
    if npdt is None:
        return None
    arr = np.frombuffer(bytes(block.data()), dtype=npdt)
    if arr.size != win["ncols"] * win["nrows"]:
        return None
    arr = arr.reshape(win["nrows"], win["ncols"]).astype(np.float64)
    if block.hasNoDataValue():
        ndv = block.noDataValue()
        arr[np.abs(arr - ndv) <= 1e-6 * max(1.0, abs(ndv))] = np.nan
    arr[~np.isfinite(arr)] = np.nan
    arr[np.abs(arr) > 1e30] = np.nan
    return arr


# ----------------------------------------------------------------------
# 6) Delineado de cuencas de aporte + hidrograma por outlet
# ----------------------------------------------------------------------
def delineate_natural_basins(dem_layer, landuse_layer, cn_table, kept_segments, outlets,
                             aoi_geom, blocks_union, hyeto_fn, feedback, ia_ratio=0.2):
    """
    outlets: lista de dicts con "point" (QgsPointXY) y "node" (id SWMM destino).
    hyeto_fn(xs, ys) -> (texto de pluviometros/pesos, dict de parse_hyetograph | None);
    recibe los centros de celda (arrays) de la cuenca para ponderar por Thiessen
    Devuelve lista de registros (uno por outlet con cuenca valida):
        node, outlet_pt, geom, area_ha, L_m, H_m, tc_min, cn, p_mm, pe_mm, qp_unit,
        qmax, gage, dt_h, dated, t_h (np), series (np, m3/s por mm efectivo total),
        pe_total (Scale Factor), q (np, m3/s)
    """
    if not outlets:
        return []
    win = _dem_window(dem_layer, aoi_geom)
    if win is None:
        feedback.pushWarning("El AOI no intersecta el DEM: no se pueden delinear cuencas naturales.")
        return []
    x0, y0, dx, dy = win["x0"], win["y0"], win["dx"], win["dy"]
    nrows, ncols = win["nrows"], win["ncols"]
    if nrows * ncols > MAX_WINDOW_CELLS:
        feedback.pushWarning(
            f"Ventana de DEM muy grande ({nrows * ncols / 1e6:.1f} M celdas): el delineado puede tardar mucho.")
    feedback.pushInfo(f"  Delineado D8 sobre ventana DEM de {ncols} x {nrows} celdas...")

    z0 = _read_raster(dem_layer.dataProvider(), win)
    if z0 is None:
        feedback.pushWarning("No se pudo leer el DEM como matriz (tipo de dato no soportado).")
        return []

    inside = polygon_mask(aoi_geom, x0, y0, dx, dy, nrows, ncols)
    valid = np.isfinite(z0) & inside
    if blocks_union is not None:
        valid &= ~polygon_mask(blocks_union, x0, y0, dx, dy, nrows, ncols)
    if not valid.any():
        feedback.pushWarning("El dominio natural (AOI menos manzanas) no tiene celdas de DEM validas.")
        return []

    Wp = ncols + 2
    z0p = np.pad(np.where(valid, z0, np.nan), 1, constant_values=np.nan)
    vp = np.isfinite(z0p)
    outside_p = np.pad(~inside, 1, constant_values=True)

    # burn: bajar el DEM bajo los ejes conservados para que el flujo los siga
    burn = np.zeros(z0p.shape, dtype=bool)
    seg_cells = {}
    step = 0.5 * min(dx, dy)
    for seg in kept_segments:
        dens = seg["geom"].densifyByDistance(step)
        pl = dens.asPolyline() or [p for part in dens.asMultiPolyline() for p in part]
        if not pl:
            continue
        xs = np.array([p.x() for p in pl])
        ys = np.array([p.y() for p in pl])
        cc = np.floor((xs - x0) / dx).astype(int)
        rr = np.floor((y0 - ys) / dy).astype(int)
        ok = (cc >= 0) & (cc < ncols) & (rr >= 0) & (rr < nrows)
        burn[rr[ok] + 1, cc[ok] + 1] = True
        seg_cells[seg["seg_id"]] = ((rr[ok] + 1) * Wp + (cc[ok] + 1)).astype(np.int64)
    burn &= vp
    zb = z0p.copy()
    zb[burn] -= BURN_DEPTH_M

    filled = _fill_depressions(zb, vp)
    dirs = _d8(filled, vp, dx, dy)
    recv = _receivers(dirs, Wp)
    valid_flat = vp.ravel()
    acc = _flow_accumulation(recv, valid_flat)
    donors_sorted, starts, counts = _donor_csr(recv)
    dirs_flat = dirs.ravel()
    step_len = _step_lengths(dx, dy)
    z0_flat = z0p.ravel()
    burn_flat = burn.ravel()
    outside_flat = outside_p.ravel()
    nb_offs = np.array([dr * Wp + dc for dr, dc in _NBR], dtype=np.int64)
    claimed = np.zeros(z0_flat.size, dtype=bool)

    if not cn_table:
        feedback.pushWarning(
            f"Sin columna 'cn' en la tabla de cobertura: se usa CN {DEFAULT_CN} en todas las cuencas.")
    lu_flat = None
    if landuse_layer is not None:
        lu = _read_raster(landuse_layer.dataProvider(), win)
        if lu is not None:
            lu_flat = np.pad(lu, 1, constant_values=np.nan).ravel()

    # pre-paso: ajustar cada outlet a su celda y procesar de menor a mayor acumulacion
    # (las cuencas de aguas arriba se descuentan de las que las contienen)
    snapped = []
    R = SNAP_RADIUS_CELLS
    for og in outlets:
        pt = og["point"]
        rc = int(math.floor((y0 - pt.y()) / dy)) + 1
        cc_ = int(math.floor((pt.x() - x0) / dx)) + 1
        best, best_key = None, None
        # 1) preferente: la celda de mayor acumulacion entre las celdas de TODA la red que
        #    drena a este outlet (recoge el flujo aunque el agua se deslice por el borde
        #    de las manzanas antes de salir del dominio)
        net_ids = [i for i in og.get("net_seg_ids", []) if i in seg_cells]
        if net_ids:
            cand = np.unique(np.concatenate([seg_cells[i] for i in net_ids]))
            cand = cand[valid_flat[cand]]
            if cand.size:
                best = int(cand[np.argmax(acc[cand])])
                snapped.append((og, best))
                continue
        # 2) respaldo: celda de mayor acumulacion en un radio de pocas celdas
        for r in range(rc - R, rc + R + 1):
            for c in range(cc_ - R, cc_ + R + 1):
                if 1 <= r <= nrows and 1 <= c <= ncols and vp[r, c]:
                    i = r * Wp + c
                    key = (acc[i], burn_flat[i], -((r - rc) ** 2 + (c - cc_) ** 2))
                    if best_key is None or key > best_key:
                        best, best_key = i, key
        snapped.append((og, best))
    snapped.sort(key=lambda t: acc[t[1]] if t[1] is not None else -1.0)

    records = []
    n_no_runoff = 0
    for og, best in snapped:
        if feedback.isCanceled():
            break
        pt = og["point"]
        label = f"({pt.x():.1f}, {pt.y():.1f})"
        if best is None:
            feedback.pushWarning(f"Outlet {label}: sin celdas de DEM validas cerca; se omite.")
            continue

        idx, dist = _basin_cells(best, donors_sorted, starts, counts, dirs_flat, step_len)
        if idx.size < 4:
            feedback.pushWarning(f"Outlet {label}: cuenca de aporte demasiado chica ({idx.size} celdas); se omite.")
            continue
        overlap = claimed[idx]
        if overlap.any():
            feedback.pushInfo(
                f"  Outlet {label}: se descuentan {int(overlap.sum())} celdas ya asignadas a cuencas de aguas arriba.")
            idx, dist = idx[~overlap], dist[~overlap]
            if idx.size < 4:
                feedback.pushWarning(f"Outlet {label}: la cuenca propia queda vacia tras el descuento; se omite.")
                continue
        claimed[idx] = True

        rows_i = idx // Wp - 1
        cols_i = idx % Wp - 1
        bx = x0 + ((best % Wp) - 1 + 0.5) * dx
        by = y0 - ((best // Wp) - 1 + 0.5) * dy
        snap_m = math.hypot(bx - pt.x(), by - pt.y())
        area_m2 = idx.size * dx * dy
        L = float(dist.max())
        elev = z0_flat[idx]
        H = float(np.nanmax(elev) - z0_flat[best])
        touches_aoi = bool(outside_flat[(idx[:, None] + nb_offs[None, :])].any())
        if touches_aoi:
            feedback.pushInfo(f"  Outlet {label}: la divisoria alcanza el borde del AOI (cuenca recortada por el AOI).")
        if L < min(dx, dy):
            L = min(dx, dy)
        if H <= 0.0:
            feedback.pushWarning(f"Outlet {label}: H <= 0 (cuenca plana); se usa H = 0.1 m.")
            H = 0.1

        # CN ponderado por area segun cobertura
        cn = DEFAULT_CN
        if lu_flat is not None and cn_table:
            cls = lu_flat[idx]
            cls = cls[np.isfinite(cls)]
            if cls.size:
                vals, cnts = np.unique(np.rint(cls).astype(int), return_counts=True)
                w = np.array([cn_table.get(int(v), np.nan) for v in vals])
                ok = np.isfinite(w)
                if ok.any():
                    cn = float((w[ok] * cnts[ok]).sum() / cnts[ok].sum())
                else:
                    feedback.pushWarning(f"Outlet {label}: ninguna clase de cobertura tiene CN; se usa CN {DEFAULT_CN}.")
            else:
                feedback.pushWarning(f"Outlet {label}: sin cobertura en la cuenca; se usa CN {DEFAULT_CN}.")

        gage, hy = hyeto_fn(x0 + (cols_i + 0.5) * dx, y0 - (rows_i + 0.5) * dy)
        if hy is None:
            feedback.pushWarning(f"Outlet {label}: hietograma no utilizable; se omite el aporte.")
            continue

        try:
            tc = tc_california_culverts(L, H)
        except ValueError as e:
            feedback.pushWarning(f"Outlet {label}: {e}")
            continue
        dt_min = hy["dt_h"] * 60.0
        area_km2 = area_m2 / 1e6
        uh, qp, tp, tb = scs_triangular_uh(area_km2, tc, dt_min)
        pe_inc = scs_effective_rain(hy["step_mm"], cn, ia_ratio)
        pe_total = float(pe_inc.sum())
        has_inflow = pe_total > 1e-9
        if not has_inflow:
            n_no_runoff += 1
            pe_inc = np.zeros_like(pe_inc)
        q = np.convolve(pe_inc, uh)                       # m3/s (todo 0 si no hay escorrentia)
        n = q.size
        t_h = hy["t0_h"] + np.arange(n) * hy["dt_h"]
        geom = _cells_to_polygon(rows_i, cols_i, x0, y0, dx, dy)
        feedback.pushInfo(
            f"  Cuenca -> {og['node']}: A={area_m2 / 1e4:.2f} ha, L={L:.0f} m, H={H:.1f} m, "
            f"Tc={tc:.1f} min, CN={cn:.1f}, P={hy['total_mm']:.1f} mm, Pe={pe_total:.2f} mm, "
            f"Qmax={q.max():.3f} m3/s" + (f", outlet ajustado {snap_m:.0f} m" if snap_m > 3 * max(dx, dy) else "") + ("" if has_inflow else " (SIN escorrentia: P <= Ia; no se crea INFLOW)") + ".")
        records.append({
            "node": og["node"], "outlet_pt": pt, "geom": geom,
            "area_ha": area_m2 / 1e4, "L_m": L, "H_m": H, "tc_min": tc, "cn": cn,
            "p_mm": hy["total_mm"], "pe_mm": pe_total, "qp_unit": qp, "qmax": float(q.max()),
            "gage": gage, "dt_h": hy["dt_h"], "dated": hy["dated"],
            "t_h": t_h, "series": (q / pe_total) if has_inflow else q, "pe_total": pe_total, "q": q,
            "inflow": has_inflow, "ia_mm": ia_ratio * (25400.0 / min(max(cn, 1.0), 100.0) - 254.0),
        })
    if n_no_runoff:
        feedback.pushWarning(
            f"{n_no_runoff} de {len(records)} cuenca(s) no generan escorrentia (lluvia total <= Ia): "
            f"no tienen INFLOW, pero se dibujan en la capa de cuencas. Revisa el hietograma, el CN "
            f"o la relacion Ia/S ({ia_ratio}).")
    return records


# ----------------------------------------------------------------------
# 7) Ensamblado de INFLOWS y series para el .inp
# ----------------------------------------------------------------------
def _label_time(t_h, dated):
    if dated:
        d = datetime(1970, 1, 1) + timedelta(seconds=int(round(t_h * 3600.0)))
        tm = d.strftime("%H:%M") if d.second == 0 else d.strftime("%H:%M:%S")
        return d.strftime("%m/%d/%Y"), tm
    total = int(round(t_h * 3600.0))
    h, rem = divmod(total, 3600)
    m, s = divmod(rem, 60)
    return "", (f"{h}:{m:02d}" if s == 0 else f"{h}:{m:02d}:{s:02d}")


def assemble_node_inflows(records, feedback=None):
    """Agrupa los registros por nodo. Un nodo con un solo outlet usa la serie
    unitaria (m3/s por mm) con Scale Factor = lluvia efectiva total. Un nodo con
    varios outlets recibe una unica serie (suma de los caudales, interpolada) con
    Scale Factor 1.0, porque SWMM admite un solo INFLOW FLOW por nodo.
    Devuelve (inflows, timeseries):
        inflows    = [{"node", "series", "scale"}]
        timeseries = {nombre: [(fecha, hora, valor), ...]}"""
    by_node = {}
    for r in records:
        if not r.get("inflow", True):
            continue
        by_node.setdefault(r["node"], []).append(r)
    inflows, tsd = [], {}
    for node, recs in by_node.items():
        name = f"UH_{node}"
        if len(recs) == 1:
            r = recs[0]
            t_h, vals, scale, dated = r["t_h"], r["series"], r["pe_total"], r["dated"]
        else:
            dated = all(r["dated"] for r in recs)
            ts = []
            for r in recs:
                t = r["t_h"] if (dated or not r["dated"]) else r["t_h"] - r["t_h"][0]
                ts.append(t)
            step = min(r["dt_h"] for r in recs)
            t_a = min(t[0] for t in ts)
            t_b = max(t[-1] for t in ts)
            t_h = t_a + np.arange(int(round((t_b - t_a) / step)) + 1) * step
            vals = np.zeros(t_h.size)
            for r, t in zip(recs, ts):
                vals += np.interp(t_h, t, r["q"], left=0.0, right=0.0)
            scale = 1.0
            if feedback is not None:
                feedback.pushInfo(f"  Nodo {node}: {len(recs)} outlets confluyen; se suman sus hidrogramas (Scale Factor 1.0).")
        rows = []
        for t, v in zip(t_h, vals):
            d, tm = _label_time(float(t), dated)
            rows.append((d, tm, float(v)))
        tsd[name] = rows
        inflows.append({"node": node, "series": name, "scale": float(scale)})
    return inflows, tsd
