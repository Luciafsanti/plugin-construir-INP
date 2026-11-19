# -*- coding: utf-8 -*-
"""
inp_writer.py
-------------
Construye el texto de un archivo .inp de EPA SWMM 5.2 para un modelo de
escurrimiento superficial en calles (conductos con seccion STREET), con
soporte opcional para:
  - Multiples pluviometros (RAINGAGES) con su propia serie TIMESERIES
  - Conductos de rio con seccion irregular (TRANSECTS, formato HEC-2)
  - Nodos de almacenamiento (STORAGE) para embalses/lagunas

Este modulo NO depende de PyQGIS: recibe estructuras de datos ya calculadas
(listas de dicts) para que pueda probarse de forma aislada o reutilizarse
fuera del plugin (p.ej. en un script batch).

Unidades asumidas: sistema metrico (FLOW_UNITS = CMS), coordenadas en metros
(CRS proyectado), tiempo en el formato HH:MM.
"""

from datetime import datetime


def _fmt(value, decimals=3):
    if value is None:
        return "*"
    return f"{value:.{decimals}f}"


def build_inp(
    title,
    junctions,
    outfalls,
    conduits,
    street_categories,
    subcatchments,
    raingages,
    timeseries,
    coordinates,
    polygons,
    storages=None,
    transects=None,
    options=None,
    inflows=None,
    inflow_timeseries=None,
    extra_timeseries=None,
):
    """
    Parametros
    ----------
    title : str
    junctions : list[dict]  -> id, elev, maxdepth, initdepth, surdepth, aponded
    outfalls : list[dict]   -> id, elev, type ('FREE'|'NORMAL'), stage_data
    conduits : list[dict]   -> id, from_node, to_node, length, roughness,
                                xsection_type ('STREET'|'IRREGULAR'),
                                category (nombre en street_categories, si STREET)
                                transect (nombre en transects, si IRREGULAR)
    street_categories : dict[str, dict] -> tcrown, hcurb, sx, n, a, w, sides,
                                            tback, sback, nback
    subcatchments : list[dict] -> id, raingage, outlet, area_ha, pct_imperv,
                                   width_m, pct_slope, n_imp, n_perv,
                                   dstore_imp_mm, dstore_perv_mm, pct_zero,
                                   infil_params (dict)
    raingages : list[dict] -> name, format, interval, scf, source (=nombre
                               de una serie en `timeseries`)
    timeseries : dict[str, list[tuple(date_str, time_str, value)]]
                 Una entrada por cada serie referenciada desde `raingages`.
    coordinates : dict[node_id] -> (x, y)
    polygons : dict[subcatchment_id] -> list[(x, y), ...]
    storages : list[dict] opcional -> id, elev, ymax, y0, a1, a2, a0,
                                       ysur, fevap
    transects : dict[str, dict] opcional -> nombre -> {
                    "n_left", "n_right", "n_channel",
                    "x_left", "x_right", "station_elev": [(sta, elev), ...]
                }
    options : dict opcional para sobreescribir OPTIONS por defecto
    inflows : list[dict] opcional -> node, series (nombre en `inflow_timeseries`),
              scale (Scale Factor). Se escribe en [INFLOWS] como Direct Inflow
              (constituyente FLOW, tipo FLOW, Mfactor 1.0).
    extra_timeseries : dict[str, list[tuple(date_str, time_str, value)]] opcional.
              Series adicionales (p. ej. hietogramas ponderados por Thiessen de pluviometros
              sinteticos); se escriben en [TIMESERIES] con 6 decimales.
    inflow_timeseries : dict[str, list[tuple(date_str, time_str, value)]] opcional.
              Series de los INFLOW (hidrograma en m3/s por mm de lluvia efectiva);
              se escriben en [TIMESERIES] con 6 decimales.

    Retorna
    -------
    str : contenido completo del archivo .inp
    """
    storages = storages or []
    transects = transects or {}
    inflows = inflows or []
    inflow_timeseries = inflow_timeseries or {}
    extra_timeseries = extra_timeseries or {}

    default_options = {
        "FLOW_UNITS": "CMS",
        "INFILTRATION": "HORTON",
        "FLOW_ROUTING": "DYNWAVE",
        "LINK_OFFSETS": "DEPTH",
        "MIN_SLOPE": "0",
        "ALLOW_PONDING": "YES",
        "SKIP_STEADY_STATE": "NO",
        "START_DATE": "01/01/2024",
        "START_TIME": "00:00:00",
        "REPORT_START_DATE": "01/01/2024",
        "REPORT_START_TIME": "00:00:00",
        "END_DATE": "01/01/2024",
        "END_TIME": "06:00:00",
        "SWEEP_START": "01/01",
        "SWEEP_END": "12/31",
        "DRY_DAYS": "0",
        "REPORT_STEP": "00:05:00",
        "WET_STEP": "00:05:00",
        "DRY_STEP": "01:00:00",
        "ROUTING_STEP": "0:00:05",
        "LENGTHENING_STEP": "0",
        "VARIABLE_STEP": "0.75",
        "MINIMUM_STEP": "0.5",
        "INERTIAL_DAMPING": "PARTIAL",
        "NORMAL_FLOW_LIMITED": "BOTH",
        "FORCE_MAIN_EQUATION": "H-W",
        "MIN_SURFAREA": "1.14",
        "MAX_TRIALS": "8",
        "HEAD_TOLERANCE": "0.0015",
        "SYS_FLOW_TOL": "5",
        "LAT_FLOW_TOL": "5",
    }
    if options:
        default_options.update(options)

    lines = []
    lines.append("[TITLE]")
    lines.append(f";;Generado automaticamente - SWMM Street Builder ({datetime.now().isoformat(timespec='seconds')})")
    lines.append(title)
    lines.append("")

    # ---------------- OPTIONS ----------------
    lines.append("[OPTIONS]")
    lines.append(";;Key            Value")
    for k, v in default_options.items():
        # ancho fijo + espacio explicito: garantiza separacion aun con
        # claves largas (p.ej. NORMAL_FLOW_LIMITED, FORCE_MAIN_EQUATION)
        lines.append(f"{k:<20} {v}")
    lines.append("")

    # ---------------- RAINGAGES ----------------
    lines.append("[RAINGAGES]")
    lines.append(";;Name           Format    Interval SCF      Source")
    for rg in raingages:
        src = rg.get("source", rg["name"])
        lines.append(
            f"{rg['name']:<16}{rg.get('format','INTENSITY'):<10}"
            f"{rg.get('interval','0:05'):<9}{rg.get('scf',1.0):<9}TIMESERIES {src}"
        )
    lines.append("")

    # ---------------- JUNCTIONS ----------------
    lines.append("[JUNCTIONS]")
    lines.append(";;Name           Elevation  MaxDepth   InitDepth  SurDepth   Aponded")
    for j in junctions:
        lines.append(
            f"{j['id']:<16}{_fmt(j['elev'])}   {_fmt(j.get('maxdepth', 1.0))}   "
            f"{_fmt(j.get('initdepth', 0.0))}   {_fmt(j.get('surdepth', 0.0))}   "
            f"{_fmt(j.get('aponded', 0.0))}"
        )
    lines.append("")

    # ---------------- OUTFALLS ----------------
    lines.append("[OUTFALLS]")
    lines.append(";;Name           Elevation  Type       StageData  Gated    RouteTo")
    for o in outfalls:
        lines.append(
            f"{o['id']:<16}{_fmt(o['elev'])}   {o.get('type','FREE'):<10} "
            f"{o.get('stage_data',''):<10} NO"
        )
    lines.append("")

    # ---------------- STORAGE ----------------
    # Formato FUNCTIONAL: Name Elev Ymax Y0 A1 A2 A0 (Ysur Fevap Psi Ksat IMD)
    # Area(d) = A0 + A1*d^A2  -> con A1=0, A2=0 se modela area de espejo
    # de agua constante (aproximacion prismatica, sin batimetria).
    if storages:
        lines.append("[STORAGE]")
        lines.append(";;Name           Elev       Ymax       Y0         Type      A1         A2         A0         Ysur      Fevap")
        for s in storages:
            lines.append(
                f"{s['id']:<16}{_fmt(s['elev'])}   {_fmt(s.get('ymax',2.0))}   {_fmt(s.get('y0',0.0))}   "
                f"{'FUNCTIONAL':<9} {_fmt(s.get('a1',0.0))}   {_fmt(s.get('a2',0.0))}   {_fmt(s.get('a0',1000.0),2)}   "
                f"{_fmt(s.get('ysur',0.0))}  {_fmt(s.get('fevap',0.0),2)}"
            )
        lines.append("")

    # ---------------- CONDUITS ----------------
    lines.append("[CONDUITS]")
    lines.append(";;Name           From Node        To Node          Length     Roughness  InOffset   OutOffset  InitFlow   MaxFlow")
    for c in conduits:
        lines.append(
            f"{c['id']:<16}{c['from_node']:<17}{c['to_node']:<17}"
            f"{_fmt(c['length'])}   {_fmt(c.get('roughness',0.016),3)}      0          0          "
            f"{_fmt(c.get('init_flow',0.0))}   {_fmt(c.get('max_flow',0.0))}"
        )
    lines.append("")

    # ---------------- XSECTIONS ----------------
    # Dos formatos posibles por conducto (Seccion D del manual SWMM 5.2):
    #   Link  STREET     NombreCategoria   (calles, geometria en [STREETS])
    #   Link  IRREGULAR  NombreTransecto   (rios, geometria en [TRANSECTS])
    lines.append("[XSECTIONS]")
    lines.append(";;Link           Shape        Tsect/Street")
    for c in conduits:
        xtype = c.get("xsection_type", "STREET")
        if xtype == "IRREGULAR":
            lines.append(f"{c['id']:<16}{'IRREGULAR':<13}{c['transect']}")
        else:
            lines.append(f"{c['id']:<16}{'STREET':<13}{c['category']}")
    lines.append("")

    # ---------------- STREETS ----------------
    if street_categories:
        lines.append("[STREETS]")
        lines.append(";;Name           Tcrown     Hcurb      Sx         nRoad      a          W          Sides      Tback      Sback      nBack")
        for name, p in street_categories.items():
            lines.append(
                f"{name:<16}{_fmt(p['tcrown'])}   {_fmt(p['hcurb'])}   {_fmt(p['sx'],4)}   "
                f"{_fmt(p['n'],4)}   {_fmt(p.get('a',1.5))}   {_fmt(p['w'])}   "
                f"{p.get('sides',2)}          {_fmt(p.get('tback',0.0))}   {_fmt(p.get('sback',0.0),4)}   {_fmt(p.get('nback',p['n']),4)}"
            )
        lines.append("")

    # ---------------- TRANSECTS (rios / canales de seccion irregular) ------
    # Formato HEC-2:
    #   NC  Nleft Nright Nchanl
    #   X1  Name Nsta Xleft Xright Lfactor Wfactor Eoffset  0 0 0
    #   GR  Elev Station  Elev Station  ... (hasta 4 pares por linea)
    if transects:
        lines.append("[TRANSECTS]")
        lines.append(";;Formato HEC-2: NC (rugosidades) / X1 (cabecera) / GR (estacion-elevacion)")
        for name, t in transects.items():
            lines.append(
                f"NC  {_fmt(t.get('n_left',0.05),3)}  {_fmt(t.get('n_right',0.05),3)}  {_fmt(t.get('n_channel',0.035),3)}"
            )
            sta_elev = t["station_elev"]
            nsta = len(sta_elev)
            x_left = t.get("x_left", sta_elev[0][0])
            x_right = t.get("x_right", sta_elev[-1][0])
            lines.append(
                f"X1  {name}  {nsta}  {_fmt(x_left,2)}  {_fmt(x_right,2)}  "
                f"{_fmt(t.get('lfactor',0.0),2)}  {_fmt(t.get('wfactor',0.0),2)}  {_fmt(t.get('eoffset',0.0),2)}  0  0  0"
            )
            # GR: hasta 4 pares Elev-Station por linea (orden Elev, Station)
            pairs = [f"{_fmt(elev,3)} {_fmt(sta,2)}" for sta, elev in sta_elev]
            for i in range(0, len(pairs), 4):
                lines.append("GR  " + "  ".join(pairs[i:i + 4]))
        lines.append("")

    # ---------------- SUBCATCHMENTS ----------------
    lines.append("[SUBCATCHMENTS]")
    lines.append(";;Name           Raingage         Outlet           Area       %Imperv    Width      %Slope     CurbLen  SnowPack")
    for s in subcatchments:
        lines.append(
            f"{s['id']:<16}{s.get('raingage','RG1'):<17}{s['outlet']:<17}"
            f"{_fmt(s['area_ha'])}   {_fmt(s['pct_imperv'],1)}      {_fmt(s['width_m'])}   "
            f"{_fmt(s['pct_slope'],2)}      0"
        )
    lines.append("")

    # ---------------- SUBAREAS ----------------
    lines.append("[SUBAREAS]")
    lines.append(";;Subcatchment   N-Imperv   N-Perv     S-Imperv   S-Perv     PctZero    RouteTo")
    for s in subcatchments:
        lines.append(
            f"{s['id']:<16}{_fmt(s.get('n_imp',0.013),3)}      {_fmt(s.get('n_perv',0.24),3)}      "
            f"{_fmt(s.get('dstore_imp_mm',2.0))}   {_fmt(s.get('dstore_perv_mm',5.0))}   "
            f"{_fmt(s.get('pct_zero',25.0),1)}      OUTLET"
        )
    lines.append("")

    # ---------------- INFILTRATION (Horton por defecto) ----------------
    lines.append("[INFILTRATION]")
    lines.append(";;Subcatchment   Param1     Param2     Param3     Param4     Param5")
    for s in subcatchments:
        ip = s.get("infil_params", {"f0": 75, "fmin": 13, "decay": 4, "drytime": 7, "maxinfil": 0})
        lines.append(
            f"{s['id']:<16}{_fmt(ip.get('f0',75))}   {_fmt(ip.get('fmin',13))}   "
            f"{_fmt(ip.get('decay',4))}   {_fmt(ip.get('drytime',7))}   {_fmt(ip.get('maxinfil',0))}"
        )
    lines.append("")

    # ---------------- INFLOWS (Direct Inflow, p. ej. drenaje natural) ----------------
    if inflows:
        lines.append("[INFLOWS]")
        lines.append(";;Node           Constituent      Time Series      Type     Mfactor  Sfactor  Baseline Pattern")
        for inf in inflows:
            lines.append(
                f"{inf['node']:<16}{'FLOW':<17}{inf['series']:<17}{'FLOW':<9}"
                f"{'1.0':<9}{_fmt(inf.get('scale', 1.0), 4):<9}0"
            )
        lines.append("")

    # ---------------- TIMESERIES ----------------
    lines.append("[TIMESERIES]")
    lines.append(";;Name           Date       Time       Value")
    for ts_name, values in timeseries.items():
        for date_str, time_str, value in values:
            if date_str:
                lines.append(f"{ts_name:<16}{date_str:<11}{time_str:<11}{_fmt(value,2)}")
            else:
                lines.append(f"{ts_name:<16}{time_str:<11}{_fmt(value,2)}")
    for ts_name, values in [*extra_timeseries.items(), *inflow_timeseries.items()]:
        for date_str, time_str, value in values:
            if date_str:
                lines.append(f"{ts_name:<16}{date_str:<11}{time_str:<11}{_fmt(value,6)}")
            else:
                lines.append(f"{ts_name:<16}{time_str:<11}{_fmt(value,6)}")
    lines.append("")

    # ---------------- REPORT ----------------
    lines.append("[REPORT]")
    lines.append("INPUT      NO")
    lines.append("CONTROLS   NO")
    lines.append("SUBCATCHMENTS ALL")
    lines.append("NODES      ALL")
    lines.append("LINKS      ALL")
    lines.append("")

    # ---------------- COORDINATES ----------------
    lines.append("[COORDINATES]")
    lines.append(";;Node           X-Coord            Y-Coord")
    for node_id, (x, y) in coordinates.items():
        lines.append(f"{node_id:<16}{_fmt(x,3)}         {_fmt(y,3)}")
    lines.append("")

    # ---------------- POLYGONS (subcuencas) ----------------
    lines.append("[Polygons]")
    lines.append(";;Subcatchment   X-Coord            Y-Coord")
    for sub_id, coords in polygons.items():
        for (x, y) in coords:
            lines.append(f"{sub_id:<16}{_fmt(x,3)}         {_fmt(y,3)}")
    lines.append("")

    lines.append("[TAGS]")
    lines.append("")

    return "\n".join(lines)
