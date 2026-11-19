# -*- coding: utf-8 -*-
"""
rain_thiessen.py
-----------------
Asignacion de lluvia por poligonos de Thiessen (Voronoi) para subcuencas urbanas
y cuencas de drenaje natural.

SWMM admite UN solo pluviometro por subcuenca, asi que la ponderacion se resuelve
antes de escribir el .inp:

  1. Se arma el diagrama de Voronoi de las estaciones validas.
  2. Para cada subcuenca (poligono) o cuenca natural (celdas) se calcula la
     fraccion de area que cae en cada estacion (los pesos).
  3. Los pesos se cuantizan (por defecto en pasos de 5 %) para limitar la
     cantidad de combinaciones distintas.
  4. Si queda una sola estacion, se usa el pluviometro real. Si quedan varias,
     se genera un hietograma ponderado  sum_i w_i * P_i(t)  (volumen por paso, mm)
     sobre una grilla de tiempo comun (el paso mas fino) conservando la masa.
     Las subcuencas urbanas reciben un pluviometro sintetico RG_TH<n> (formato
     VOLUME) con su serie TS_TH<n>; las cuencas naturales usan el hietograma
     ponderado directamente para el HU y la lluvia efectiva (SCS-CN).
"""

import math

import numpy as np

from qgis.core import QgsGeometry, QgsPointXY, QgsRectangle, QgsSpatialIndex, QgsFeature

from .natural_drainage import _label_time


# ----------------------------------------------------------------------
# Utilidades de series
# ----------------------------------------------------------------------
def interval_str(dt_h):
    """Intervalo de registro para [RAINGAGES]: H:MM si son minutos enteros, si no horas decimales."""
    minutes = dt_h * 60.0
    if abs(minutes - round(minutes)) < 1e-6:
        m = int(round(minutes))
        return f"{m // 60}:{m % 60:02d}"
    return f"{dt_h:.5f}"


def quantize_weights(weights, step=0.05):
    """Cuantiza pesos {nombre: w} a multiplos de `step` (metodo del mayor resto)
    sumando exactamente 1. Devuelve tupla ordenada ((nombre, w), ...) sin ceros."""
    units = max(int(round(1.0 / step)), 1)
    names = sorted(weights)
    raw = np.array([weights[n] for n in names], dtype=float)
    tot = raw.sum()
    if tot <= 0:
        return ()
    raw = raw / tot * units
    base = np.floor(raw + 1e-9).astype(int)
    left = units - int(base.sum())
    if left > 0:
        order = np.argsort(-(raw - base), kind="stable")
        for i in order[:left]:
            base[i] += 1
    return tuple((n, b / units) for n, b in zip(names, base) if b > 0)


def combine_hyetographs(items):
    """items: [(hy, peso)] con hy de natural_drainage.parse_hyetograph.
    Devuelve un hietograma ponderado (misma estructura) sobre la grilla comun."""
    if len(items) == 1:
        return dict(items[0][0])
    dated = all(h["dated"] for h, _ in items)
    mixed = (not dated) and any(h["dated"] for h, _ in items)
    starts = [0.0 if mixed else h["t0_h"] for h, _ in items]
    dts = [h["dt_h"] for h, _ in items]
    dtg = min(dts)
    ends = [s + len(h["step_mm"]) * h["dt_h"] for s, (h, _) in zip(starts, items)]
    g0, g1 = min(starts), max(ends)
    n = max(int(math.ceil((g1 - g0) / dtg - 1e-6)), 1)
    edges = g0 + np.arange(n + 1) * dtg
    vol = np.zeros(n)
    for s, (h, w) in zip(starts, items):
        m = len(h["step_mm"])
        e = s + np.arange(m + 1) * h["dt_h"]
        c = np.concatenate([[0.0], np.cumsum(h["step_mm"])])
        cg = np.interp(edges, e, c, left=0.0, right=c[-1])
        vol += w * np.diff(cg)
    return {"dt_h": dtg, "t0_h": g0, "dated": bool(dated), "step_mm": vol,
            "total_mm": float(vol.sum()), "irregular": False, "mixed": mixed}


def hyetograph_to_rows(hy):
    """Filas (fecha, hora, mm por paso) para [TIMESERIES] (formato VOLUME)."""
    rows = []
    for k, v in enumerate(hy["step_mm"]):
        d, tm = _label_time(hy["t0_h"] + k * hy["dt_h"], hy["dated"])
        rows.append((d, tm, float(v)))
    return rows


# ----------------------------------------------------------------------
# Asignador Thiessen
# ----------------------------------------------------------------------
class RainAssigner:
    """stations: [{"name": nombre del pluviometro SWMM, "pt": QgsPointXY, "hy": dict parse_hyetograph}]"""

    def __init__(self, stations, aoi_geom, step_pct, feedback):
        self.feedback = feedback
        self.step = min(max(step_pct, 1.0), 50.0) / 100.0
        # duplicados exactos: se conserva el primero
        seen, st = set(), []
        for s in stations:
            key = (round(s["pt"].x(), 6), round(s["pt"].y(), 6))
            if key in seen:
                feedback.pushWarning(
                    f"Pluviometro {s['name']}: misma ubicacion que otro; se omite en Thiessen.")
                continue
            seen.add(key)
            st.append(s)
        self.stations = st
        self.sx = np.array([s["pt"].x() for s in st])
        self.sy = np.array([s["pt"].y() for s in st])
        self.hy_by_name = {s["name"]: s["hy"] for s in st}
        self._hy_cache = {}
        self._urban = {}       # key -> nombre RG_TH<n>
        self.cells = None
        self._cell_idx = None
        self._build_cells(aoi_geom)
        irregular = [s["name"] for s in st if s["hy"].get("irregular")]
        if irregular:
            feedback.pushWarning(
                "Hietogramas con paso no uniforme (se usa el paso mediano): " + ", ".join(irregular))

    # -- poligonos de Thiessen --------------------------------------------
    def _build_cells(self, aoi_geom):
        if len(self.stations) < 2:
            return
        bb = aoi_geom.boundingBox()
        bb.combineExtentWith(QgsRectangle(self.sx.min(), self.sy.min(), self.sx.max(), self.sy.max()))
        pad = max(bb.width(), bb.height(), 1000.0)
        bb.grow(pad)
        ext = QgsGeometry.fromRect(bb)
        mp = QgsGeometry.fromMultiPointXY([s["pt"] for s in self.stations])
        try:
            vor = mp.voronoiDiagram(ext, 0.0, False)
            parts = vor.asGeometryCollection() if vor is not None and not vor.isEmpty() else []
        except Exception:
            parts = []
        cells = []
        for s in self.stations:
            pg = QgsGeometry.fromPointXY(s["pt"])
            match = next((p for p in parts if p.contains(pg)), None)
            if match is None:
                cells = []
                break
            cells.append((s["name"], match))
        if len(cells) != len(self.stations):
            self.feedback.pushWarning(
                "No se pudo construir el diagrama de Thiessen; se asigna el pluviometro mas cercano.")
            return
        self.cells = cells
        self._cell_idx = QgsSpatialIndex()
        for i, (_n, g) in enumerate(cells):
            f = QgsFeature()
            f.setId(i)
            f.setGeometry(g)
            self._cell_idx.addFeature(f)
        self.feedback.pushInfo(f"  Thiessen: {len(cells)} pluviometro(s); pesos en pasos de {self.step * 100:.0f} %.")

    def cell_polygons(self, aoi_geom):
        """[(nombre, geometria recortada al AOI, area_ha)] para la capa opcional."""
        out = []
        if self.cells is None:
            return out
        for name, g in self.cells:
            c = g.intersection(aoi_geom)
            if c is not None and not c.isEmpty():
                out.append((name, c, c.area() / 1e4))
        return out

    # -- pesos --------------------------------------------------------------
    def _nearest_name(self, x, y):
        i = int(np.argmin((self.sx - x) ** 2 + (self.sy - y) ** 2))
        return self.stations[i]["name"]

    def _weights_geom(self, geom):
        if len(self.stations) == 1:
            return {self.stations[0]["name"]: 1.0}
        c = geom.centroid().asPoint()
        if self.cells is None:
            return {self._nearest_name(c.x(), c.y()): 1.0}
        cand = self._cell_idx.intersects(geom.boundingBox())
        if not cand:
            return {self._nearest_name(c.x(), c.y()): 1.0}
        if len(cand) == 1:
            return {self.cells[cand[0]][0]: 1.0}
        w = {}
        try:
            g = geom if geom.isGeosValid() else geom.makeValid()
            for i in cand:
                a = g.intersection(self.cells[i][1]).area()
                if a > 0.0:
                    w[self.cells[i][0]] = w.get(self.cells[i][0], 0.0) + a
        except Exception:
            w = {}
        tot = sum(w.values())
        if tot <= 0.0:
            return {self._nearest_name(c.x(), c.y()): 1.0}
        return {k: v / tot for k, v in w.items()}

    def _weights_cells(self, xs, ys):
        if len(self.stations) == 1:
            return {self.stations[0]["name"]: 1.0}
        best_d = np.full(xs.shape, np.inf)
        best_i = np.zeros(xs.shape, dtype=np.int64)
        for i in range(len(self.stations)):
            d = (xs - self.sx[i]) ** 2 + (ys - self.sy[i]) ** 2
            m = d < best_d
            best_d[m] = d[m]
            best_i[m] = i
        cnt = np.bincount(best_i, minlength=len(self.stations)).astype(float)
        cnt /= cnt.sum()
        return {self.stations[i]["name"]: float(cnt[i]) for i in range(len(self.stations)) if cnt[i] > 0}

    @staticmethod
    def label(key):
        return "|".join(f"{n}:{w:.2f}" for n, w in key)

    def _hy_for_key(self, key):
        if key not in self._hy_cache:
            items = [(self.hy_by_name[n], w) for n, w in key]
            hy = combine_hyetographs(items)
            if hy.get("mixed"):
                self.feedback.pushWarning(
                    f"Thiessen {self.label(key)}: hietogramas con y sin fecha; se alinean los primeros registros.")
            self._hy_cache[key] = hy
        return self._hy_cache[key]

    # -- API -----------------------------------------------------------------
    def gage_for_geometry(self, geom):
        """Subcuenca urbana -> (nombre de pluviometro SWMM, texto de pesos)."""
        key = quantize_weights(self._weights_geom(geom), self.step)
        if len(key) == 1:
            return key[0][0], f"{key[0][0]}:1.00"
        if key not in self._urban:
            self._urban[key] = f"RG_TH{len(self._urban) + 1}"
        return self._urban[key], self.label(key)

    def hyetograph_for_cells(self, xs, ys):
        """Cuenca natural (centros de celda) -> (texto de pesos, hietograma ponderado)."""
        key = quantize_weights(self._weights_cells(np.asarray(xs), np.asarray(ys)), self.step)
        return self.label(key), self._hy_for_key(key)

    def synthetic_gages(self):
        """Pluviometros sinteticos usados por subcuencas urbanas -> (raingages, timeseries)."""
        gages, ts = [], {}
        for key, name in self._urban.items():
            hy = self._hy_for_key(key)
            ts_name = "TS_" + name[3:]
            gages.append({"name": name, "format": "VOLUME", "interval": interval_str(hy["dt_h"]),
                          "scf": 1.0, "source": ts_name})
            ts[ts_name] = hyetograph_to_rows(hy)
        return gages, ts
