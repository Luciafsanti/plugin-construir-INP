# -*- coding: utf-8 -*-
"""
swmm_algorithm.py
------------------
Algoritmo de QGIS Processing: "SWMM Street Builder"

Construye un modelo EPA SWMM 5.2 de escurrimiento superficial en calles
(sistema mayor) a partir de:
    DEM, poligono de area de interes, ejes de calle, manzanas y uso de suelo.

Opcionalmente incorpora:
    - Multiples pluviometros (shp de puntos + tabla CSV de asociacion)
    - Rios (shp lineal) y embalses/lagunas (shp de poligonos) como SALIDAS
      del modelo (OUTFALL en el extremo de calle que llega a ellos; no son
      conductos ni Storage Units)
    - Red de drenaje natural (shp lineal): sus outlets aportan a la red de
      calles como INFLOW con hidrograma unitario sintetico SCS (sin modelar
      cauces ni cuencas naturales como elementos SWMM)

SUPUESTOS / LIMITACIONES (documentados explicitamente, revisar antes de usar
en produccion):

 1. La capa de ejes de calle NO necesita venir pre-noded: el algoritmo la
    parte automaticamente en cada cruce (native:splitwithlines contra si
    misma) y luego agrupa (clusteriza) los extremos cercanos en nodos segun
    NODE_SNAP_TOLERANCE. Los rios NO entran en la red (ver 3).
 2. La cota de cada nodo se obtiene del DEM promediando el punto central y
    4 puntos vecinos (patron en cruz) para reducir el ruido de un DEM de
    resolucion gruesa (p.ej. 5 m) en cruces de calle.
 3. Los OUTFALLS son nodos de grado 1 (extremos de calle) a una distancia
    <= OUTFALL_BUFFER de un rio o de un embalse/laguna (recortados al AOI).
    Los poligonos de embalse se restan a las calles antes del noding, para
    que terminen en la orilla. Un OUTFALL admite un solo conducto: una calle
    que cruza un rio sin terminar ahi no genera outfall. Si no hay ninguno,
    se usa el nodo de grado 1 de menor cota. Revisar tras generar el modelo.
 4. La subdivision de manzanas grandes (> MAX_BLOCK_AREA_HA) se hace con una
    grilla rectangular regular recortada por la manzana. Es una aproximacion
    geometrica simple, NO sigue necesariamente las lineas naturales de
    escurrimiento interno del predio. Revisar/ajustar manualmente en casos
    criticos.
 5. El "Width" (ancho de subcuenca) se aproxima como Area / (coef * Perimetro).
    Es la aproximacion clasica cuando no se cuenta con longitud de
    escurrimiento medida en campo. coef es configurable (por defecto 0.5).
 6. La categoria de calle (para asignar Ts/Sx/W/n de la seccion STREET) se
    toma de un campo del usuario en la capa de ejes de calle. Si el campo no
    existe o el valor no matchea la tabla, se usa la categoria por defecto
    (la primera fila de la tabla de categorias).
 7. El %Impermeable se calcula reclasificando el raster de uso de suelo a
    un raster de %imperviousness (segun tabla de lookup) y promediando por
    zona (subcuenca).
 8. MaxDepth de cada JUNCTION conectado a calles = minimo Hcurb entre las
    categorias de calle de los conductos que llegan a ese nodo (el agua
    desborda por el cordon mas bajo). Nodos que solo conectan a rios usan
    un MaxDepth generico mas amplio (no hay cordon fisico).
 9. Rios y embalses/lagunas no se modelan (ni conductos, ni transectos, ni
    Storage Units): son solo referencia de salida (nd.SinkRef).
10. Un drenaje natural cuyo outlet cae a <= DRAIN_CONNECT_BUFFER de un rio o
    embalse (y mas cerca que el extremo de calle mas cercano) desemboca en
    una salida: se descarta (sin INFLOW) y sus tramos se marcan "a_salida".
11. (Opcional) Red de drenaje natural: ver natural_drainage.py. Los ejes se
    recortan al AOI, se les descuenta lo que cae sobre manzanas, se nodan y se
    orientan segun el DEM. Solo el extremo aguas abajo (outlet) de cada drenaje
    se vincula al nodo de calle mas cercano dentro de DRAIN_CONNECT_BUFFER; los
    drenajes que no desembocan en ninguna calle se PURGAN. Los cauces NO son
    conductos y sus cuencas NO son subcuencas SWMM: por cada nodo vinculado se
    delimita la cuenca de aporte con el DEM (D8, numpy), se calcula Tc por
    California Culverts Practice (L en millas, H en pies), se genera el HU
    triangular SCS con el paso de tiempo del hietograma, se convoluciona con la
    lluvia efectiva (SCS-CN por cobertura) y se escribe como INFLOW (FLOW) con
    una serie en [TIMESERIES] (m3/s por mm efectivo) y Scale Factor = Pe total.
"""

import os
import csv
import math
from collections import defaultdict

from qgis.PyQt.QtCore import QVariant
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingException,
    QgsProcessingParameterRasterLayer,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterField,
    QgsProcessingParameterFile,
    QgsProcessingParameterNumber,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterFileDestination,
    QgsProcessingParameterString,
    QgsFields,
    QgsField,
    QgsFeature,
    QgsGeometry,
    QgsPointXY,
    QgsWkbTypes,
    QgsFeatureSink,
    QgsSpatialIndex,
    QgsVectorLayer,
)
import processing

from .inp_writer import build_inp
from . import natural_drainage as nd
from . import street_orientation as so
from . import rain_thiessen as th


PLUGIN_DIR = os.path.dirname(__file__)
DEFAULT_STREET_CATEGORIES_CSV = os.path.join(PLUGIN_DIR, "lookup_tables", "street_categories.csv")
DEFAULT_LANDUSE_IMPERV_CSV = os.path.join(PLUGIN_DIR, "lookup_tables", "landuse_imperviousness.csv")
DEFAULT_RAINFALL_CSV = os.path.join(PLUGIN_DIR, "lookup_tables", "hietograma_ejemplo.csv")

FALLBACK_JUNCTION_MAXDEPTH = 0.15    # si un nodo de calle no logra vincularse a ninguna categoria
RIVER_ONLY_NODE_MAXDEPTH = 5.0       # nodos que solo conectan a rio (sin cordon fisico)


class SwmmStreetBuilderAlgorithm(QgsProcessingAlgorithm):

    INPUT_DEM = "INPUT_DEM"
    INPUT_AOI = "INPUT_AOI"
    INPUT_STREETS = "INPUT_STREETS"
    INPUT_BLOCKS = "INPUT_BLOCKS"
    INPUT_LANDUSE = "INPUT_LANDUSE"

    FIELD_STREET_CATEGORY = "FIELD_STREET_CATEGORY"
    STREET_CATEGORIES_CSV = "STREET_CATEGORIES_CSV"
    LANDUSE_IMPERV_CSV = "LANDUSE_IMPERV_CSV"
    RAINFALL_CSV = "RAINFALL_CSV"

    INPUT_RAINGAGES = "INPUT_RAINGAGES"
    FIELD_RAINGAGE_CODE = "FIELD_RAINGAGE_CODE"
    RAINGAGES_CSV = "RAINGAGES_CSV"

    INPUT_RIVERS = "INPUT_RIVERS"
    INPUT_STORAGE = "INPUT_STORAGE"

    INPUT_DRAINS = "INPUT_DRAINS"
    DRAIN_CONNECT_BUFFER = "DRAIN_CONNECT_BUFFER"
    IA_RATIO = "IA_RATIO"
    THIESSEN_STEP = "THIESSEN_STEP"
    OUTPUT_THIESSEN = "OUTPUT_THIESSEN"

    NODE_SNAP_TOLERANCE = "NODE_SNAP_TOLERANCE"
    MAX_BLOCK_AREA_HA = "MAX_BLOCK_AREA_HA"
    WIDTH_COEF = "WIDTH_COEF"
    OUTFALL_BUFFER = "OUTFALL_BUFFER"
    MODEL_TITLE = "MODEL_TITLE"

    OUTPUT_INP = "OUTPUT_INP"
    OUTPUT_NODES = "OUTPUT_NODES"
    OUTPUT_CONDUITS = "OUTPUT_CONDUITS"
    OUTPUT_SUBCATCHMENTS = "OUTPUT_SUBCATCHMENTS"
    OUTPUT_NATURAL_BASINS = "OUTPUT_NATURAL_BASINS"
    OUTPUT_NATURAL_NETWORK = "OUTPUT_NATURAL_NETWORK"

    def createInstance(self):
        return SwmmStreetBuilderAlgorithm()

    def name(self):
        return "build_swmm_street_model"

    def displayName(self):
        return "Construir modelo SWMM de escurrimiento en calles"

    def group(self):
        return "SWMM Street Builder"

    def groupId(self):
        return "swmm_street_builder"

    def shortHelpString(self):
        return (
            "Genera un modelo EPA SWMM 5.2 (.inp) de escurrimiento superficial en "
            "calles (conductos con seccion STREET) a partir de un DEM, un poligono "
            "de area de interes, los ejes de calle, las manzanas y un raster de uso "
            "de suelo clasificado.\n\n"
            "Opcional: pluviometros multiples (shp puntos + CSV codigo/rain_format/"
            "file), rios (shp lineal) y embalses/lagunas (shp poligonos). Rios y "
            "embalses son SALIDAS del modelo: no se modelan como conductos ni Storage "
            "Units; se crea un OUTFALL en el extremo de calle (grado 1) que quede a "
            "<= la distancia indicada de ellos, y los embalses se restan a las calles.\n\n"
            "Opcional: RED DE DRENAJE NATURAL (shp lineal). Sus tramos se orientan \n"
            "segun el DEM; solo el extremo aguas abajo (outlet) se conecta al nodo \n"
            "de calle mas cercano (dentro de la distancia indicada); si desemboca en un \n"
            "rio/embalse (mas cerca que cualquier calle) no aporta (a_salida) y los \n"
            "drenajes que no desembocan en ninguna calle se eliminan. Cada outlet aporta como \n"
            "INFLOW con un hidrograma unitario triangular SCS (Tc por California \n"
            "Culverts) de la cuenca delimitada con el DEM; requiere la columna 'cn' \n"
            "en la tabla de cobertura. Los cauces y cuencas naturales no se modelan \n"
            "como elementos SWMM (capas opcionales de salida para verificarlos)."
        )

    # ------------------------------------------------------------------
    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterRasterLayer(self.INPUT_DEM, "DEM"))
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.INPUT_AOI, "Poligono de area de interes",
            [QgsProcessing.TypeVectorPolygon]))
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.INPUT_STREETS, "Ejes de calle",
            [QgsProcessing.TypeVectorLine]))
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.INPUT_BLOCKS, "Manzanas",
            [QgsProcessing.TypeVectorPolygon]))
        self.addParameter(QgsProcessingParameterRasterLayer(
            self.INPUT_LANDUSE, "Uso de suelo (raster clasificado)"))

        self.addParameter(QgsProcessingParameterField(
            self.FIELD_STREET_CATEGORY, "Campo de categoria de calle (en Ejes de calle)",
            parentLayerParameterName=self.INPUT_STREETS, optional=True))
        self.addParameter(QgsProcessingParameterFile(
            self.STREET_CATEGORIES_CSV, "Tabla de categorias de calle (CSV)",
            defaultValue=DEFAULT_STREET_CATEGORIES_CSV))
        self.addParameter(QgsProcessingParameterFile(
            self.LANDUSE_IMPERV_CSV, "Tabla uso de suelo -> %impermeable (CSV)",
            defaultValue=DEFAULT_LANDUSE_IMPERV_CSV))
        self.addParameter(QgsProcessingParameterFile(
            self.RAINFALL_CSV,
            "Hietograma por defecto (CSV: time,value o date,time,value) -- "
            "se usa solo si NO se cargan pluviometros abajo",
            defaultValue=DEFAULT_RAINFALL_CSV, optional=True))

        # ---- Pluviometros (opcional) ----
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.INPUT_RAINGAGES, "Pluviometros / Rain Gauges (shp puntos, opcional)",
            [QgsProcessing.TypeVectorPoint], optional=True))
        self.addParameter(QgsProcessingParameterField(
            self.FIELD_RAINGAGE_CODE, "Campo de codigo de pluviometro",
            parentLayerParameterName=self.INPUT_RAINGAGES, optional=True))
        self.addParameter(QgsProcessingParameterFile(
            self.RAINGAGES_CSV,
            "Tabla de pluviometros (CSV: codigo, rain_format, file) -- rain_format "
            "en {intensity, volume, cumulative}",
            optional=True))

        # ---- Rios (opcional) ----
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.INPUT_RIVERS, "Rios / cauces receptores = SALIDAS del modelo (shp lineal, opcional)",
            [QgsProcessing.TypeVectorLine], optional=True))

        # ---- Storage units (opcional) ----
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.INPUT_STORAGE, "Embalses / lagunas receptores = SALIDAS del modelo (shp poligonos, opcional)",
            [QgsProcessing.TypeVectorPolygon], optional=True))

        # ---- Drenaje natural (opcional) ----
        self.addParameter(QgsProcessingParameterFeatureSource(
            self.INPUT_DRAINS, "Red de drenaje natural (shp lineal, opcional)",
            [QgsProcessing.TypeVectorLine], optional=True))
        self.addParameter(QgsProcessingParameterNumber(
            self.NODE_SNAP_TOLERANCE, "Tolerancia de union de nodos (m)",
            type=QgsProcessingParameterNumber.Double, defaultValue=1.0, minValue=0.01))
        self.addParameter(QgsProcessingParameterNumber(
            self.MAX_BLOCK_AREA_HA, "Area maxima de subcuenca antes de subdividir (ha)",
            type=QgsProcessingParameterNumber.Double, defaultValue=1.5, minValue=0.01))
        self.addParameter(QgsProcessingParameterNumber(
            self.WIDTH_COEF, "Coeficiente para Width = Area / (coef * Perimetro)",
            type=QgsProcessingParameterNumber.Double, defaultValue=0.5, minValue=0.05, maxValue=1.0))
        self.addParameter(QgsProcessingParameterNumber(
            self.OUTFALL_BUFFER, "Distancia max. calle->rio/storage para crear outfall (m)",
            type=QgsProcessingParameterNumber.Double, defaultValue=20.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber(
            self.DRAIN_CONNECT_BUFFER,
            "Distancia maxima (buffer) para vincular un drenaje natural con un nodo de calle (m)",
            type=QgsProcessingParameterNumber.Double, defaultValue=50.0, minValue=0.0))
        self.addParameter(QgsProcessingParameterNumber(
            self.IA_RATIO,
            "Drenaje natural: relacion de abstraccion inicial Ia/S del metodo SCS-CN (estandar 0.2)",
            type=QgsProcessingParameterNumber.Double, defaultValue=0.2, minValue=0.0, maxValue=0.5))
        self.addParameter(QgsProcessingParameterNumber(
            self.THIESSEN_STEP,
            "Pluviometros: resolucion de los pesos de Thiessen (%)",
            type=QgsProcessingParameterNumber.Double, defaultValue=5.0, minValue=1.0, maxValue=50.0))
        self.addParameter(QgsProcessingParameterString(
            self.MODEL_TITLE, "Titulo del modelo", defaultValue="Modelo SWMM - escurrimiento en calles"))

        self.addParameter(QgsProcessingParameterFileDestination(
            self.OUTPUT_INP, "Archivo .inp de salida", fileFilter="SWMM input (*.inp)"))
        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT_NODES, "Capa de salida: Nodos"))
        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT_CONDUITS, "Capa de salida: Conductos (calles)"))
        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT_SUBCATCHMENTS, "Capa de salida: Subcuencas"))
        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT_THIESSEN, "Capa de salida: Poligonos de Thiessen de pluviometros (opcional)",
            QgsProcessing.TypeVectorPolygon, optional=True, createByDefault=False))
        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT_NATURAL_BASINS, "Capa de salida: Cuencas de aporte de drenaje natural (opcional, solo verificacion)",
            QgsProcessing.TypeVectorPolygon, optional=True, createByDefault=True))
        self.addParameter(QgsProcessingParameterFeatureSink(
            self.OUTPUT_NATURAL_NETWORK, "Capa de salida: Red de drenaje natural conservada/purgada + conectores (opcional)",
            QgsProcessing.TypeVectorLine, optional=True, createByDefault=True))

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------
    @staticmethod
    def _load_csv_dict(path, key_field):
        rows = {}
        with open(path, newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            for row in reader:
                rows[row[key_field]] = row
        return rows

    @staticmethod
    def _snap_key(x, y, tol):
        return (round(x / tol), round(y / tol))

    @staticmethod
    def _sample_dem_avg(dem_provider, x, y, cell_size):
        offsets = [(0, 0), (cell_size, 0), (-cell_size, 0), (0, cell_size), (0, -cell_size)]
        values = []
        for dx, dy in offsets:
            val, ok = dem_provider.sample(QgsPointXY(x + dx, y + dy), 1)
            if ok and val is not None:
                values.append(val)
        if not values:
            return None
        return sum(values) / len(values)

    @staticmethod
    def _attr(feat, fields_lower_map, name, default=None, cast=float):
        """Lectura case-insensitive de un atributo, con default si falta/vacio/invalido."""
        key = fields_lower_map.get(name.lower())
        if key is None:
            return default
        try:
            val = feat[key]
        except KeyError:
            return default
        if val is None or (isinstance(val, str) and val.strip() == ""):
            return default
        if cast is None:
            return val
        try:
            return cast(val)
        except (TypeError, ValueError):
            return default

    # ------------------------------------------------------------------
    def processAlgorithm(self, parameters, context, feedback):
        dem_layer = self.parameterAsRasterLayer(parameters, self.INPUT_DEM, context)
        streets_source = self.parameterAsSource(parameters, self.INPUT_STREETS, context)
        landuse_layer = self.parameterAsRasterLayer(parameters, self.INPUT_LANDUSE, context)

        category_field = self.parameterAsString(parameters, self.FIELD_STREET_CATEGORY, context)
        street_cat_csv = self.parameterAsFile(parameters, self.STREET_CATEGORIES_CSV, context)
        landuse_csv = self.parameterAsFile(parameters, self.LANDUSE_IMPERV_CSV, context)
        rainfall_csv = self.parameterAsFile(parameters, self.RAINFALL_CSV, context)

        raingages_source = self.parameterAsSource(parameters, self.INPUT_RAINGAGES, context)
        raingage_code_field = self.parameterAsString(parameters, self.FIELD_RAINGAGE_CODE, context)
        raingages_csv = self.parameterAsFile(parameters, self.RAINGAGES_CSV, context)

        rivers_source = self.parameterAsSource(parameters, self.INPUT_RIVERS, context)
        storage_source = self.parameterAsSource(parameters, self.INPUT_STORAGE, context)

        drains_source = self.parameterAsSource(parameters, self.INPUT_DRAINS, context)
        drain_connect_buffer = self.parameterAsDouble(parameters, self.DRAIN_CONNECT_BUFFER, context)
        ia_ratio = self.parameterAsDouble(parameters, self.IA_RATIO, context)
        thiessen_step = self.parameterAsDouble(parameters, self.THIESSEN_STEP, context)

        tol = self.parameterAsDouble(parameters, self.NODE_SNAP_TOLERANCE, context)
        max_block_area_ha = self.parameterAsDouble(parameters, self.MAX_BLOCK_AREA_HA, context)
        width_coef = self.parameterAsDouble(parameters, self.WIDTH_COEF, context)
        outfall_buffer = self.parameterAsDouble(parameters, self.OUTFALL_BUFFER, context)
        model_title = self.parameterAsString(parameters, self.MODEL_TITLE, context)

        inp_path = self.parameterAsFileOutput(parameters, self.OUTPUT_INP, context)

        if dem_layer is None or landuse_layer is None:
            raise QgsProcessingException("DEM y Uso de Suelo deben ser rasters validos.")

        has_rivers = rivers_source is not None
        has_storage = storage_source is not None
        has_drains = drains_source is not None
        has_raingages = raingages_source is not None
        if has_raingages and not raingages_csv:
            raise QgsProcessingException(
                "Se cargo una capa de Pluviometros pero no la tabla CSV asociada "
                "(codigo, rain_format, file)."
            )
        if has_raingages and not raingage_code_field:
            raise QgsProcessingException(
                "Se cargo una capa de Pluviometros pero no se indico el campo de codigo."
            )

        crs = streets_source.sourceCrs()
        dem_provider = dem_layer.dataProvider()
        cell_size = dem_layer.rasterUnitsPerPixelX()

        # ==============================================================
        # 1. AOI: materializar poligono de area de interes (disuelto)
        # ==============================================================
        feedback.pushInfo("1/10 Preparando area de interes...")
        dissolve_result = processing.run("native:dissolve", {
            "INPUT": parameters[self.INPUT_AOI], "FIELD": [], "OUTPUT": "memory:aoi_dissolved"
        }, context=context, feedback=feedback)["OUTPUT"]
        if isinstance(dissolve_result, str):
            aoi_layer_tmp = context.temporaryLayerStore().mapLayer(dissolve_result)
        else:
            aoi_layer_tmp = dissolve_result
        aoi_geom = None
        for f in aoi_layer_tmp.getFeatures():
            aoi_geom = f.geometry()
            break
        if aoi_geom is None:
            raise QgsProcessingException("No se pudo leer el poligono de area de interes.")

        # ==============================================================
        # 2. Clip de calles, manzanas, rios y storage al AOI (storage se resta a las calles)
        # ==============================================================
        feedback.pushInfo("2/10 Recortando capas al area de interes...")
        streets_clipped = processing.run("native:clip", {
            "INPUT": parameters[self.INPUT_STREETS],
            "OVERLAY": parameters[self.INPUT_AOI],
            "OUTPUT": "memory:streets_clip"
        }, context=context, feedback=feedback)["OUTPUT"]

        blocks_clipped = processing.run("native:clip", {
            "INPUT": parameters[self.INPUT_BLOCKS],
            "OVERLAY": parameters[self.INPUT_AOI],
            "OUTPUT": "memory:blocks_clip"
        }, context=context, feedback=feedback)["OUTPUT"]

        rivers_clipped = None
        if has_rivers:
            rivers_clipped = processing.run("native:clip", {
                "INPUT": parameters[self.INPUT_RIVERS],
                "OVERLAY": parameters[self.INPUT_AOI],
                "OUTPUT": "memory:rivers_clip"
            }, context=context, feedback=feedback)["OUTPUT"]

        storage_clipped = None
        if has_storage:
            storage_clipped = processing.run("native:clip", {
                "INPUT": parameters[self.INPUT_STORAGE],
                "OVERLAY": parameters[self.INPUT_AOI],
                "OUTPUT": "memory:storage_clip"
            }, context=context, feedback=feedback)["OUTPUT"]
            # Lectura B: el storage es una SALIDA; las calles terminan en su orilla.
            storage_fixed = processing.run("native:fixgeometries", {
                "INPUT": storage_clipped, "OUTPUT": "memory:storage_fixed"
            }, context=context, feedback=feedback)["OUTPUT"]
            streets_clipped = processing.run("native:difference", {
                "INPUT": streets_clipped, "OVERLAY": storage_fixed,
                "OUTPUT": "memory:streets_minus_storage"
            }, context=context, feedback=feedback)["OUTPUT"]

        # ==============================================================
        # 3. Topologia: fix, explode, auto-noding en cruces (solo calles)
        # ==============================================================
        feedback.pushInfo("3/10 Limpiando topologia y particionando en cruces...")
        streets_fixed = processing.run("native:fixgeometries", {
            "INPUT": streets_clipped, "OUTPUT": "memory:streets_fixed"
        }, context=context, feedback=feedback)["OUTPUT"]
        streets_single = processing.run("native:multiparttosingleparts", {
            "INPUT": streets_fixed, "OUTPUT": "memory:streets_single"
        }, context=context, feedback=feedback)["OUTPUT"]
        streets_noded = processing.run("native:splitwithlines", {
            "INPUT": streets_single, "LINES": streets_single, "OUTPUT": "memory:streets_noded"
        }, context=context, feedback=feedback)["OUTPUT"]
        streets_snapped_raw = processing.run("native:snapgeometries", {
            "INPUT": streets_noded, "REFERENCE_LAYER": streets_noded,
            "TOLERANCE": tol, "BEHAVIOR": 0, "OUTPUT": "memory:streets_snapped_raw"
        }, context=context, feedback=feedback)["OUTPUT"]
        # native:snapgeometries puede introducir auto-intersecciones al mover
        # vertices (sobre todo con calles muy cortas o casi paralelas dentro
        # de la tolerancia) -- se corrige con un fixgeometries adicional,
        # y una geometria invalida abortaria pasos posteriores.
        streets_snapped = processing.run("native:fixgeometries", {
            "INPUT": streets_snapped_raw, "OUTPUT": "memory:streets_snapped"
        }, context=context, feedback=feedback)["OUTPUT"]

        # ==============================================================
        # 4. Construir nodos (cluster de extremos de calles) y conductos
        # ==============================================================
        feedback.pushInfo("4/10 Generando nodos y conductos (calles)...")

        raw_segments = []  # (geom, kind, extra_dict)
        for feat in streets_snapped.getFeatures():
            geom = feat.geometry()
            if geom is None or geom.isEmpty():
                continue
            cat_val = feat[category_field] if category_field else None
            raw_segments.append((geom, "street", {"category": cat_val}))

        # ---- Salidas del modelo: rios + storages recortados al AOI (no son conductos) ----
        sink_geoms = []
        for lyr_ in (rivers_clipped, storage_clipped):
            if lyr_ is None:
                continue
            for f_ in lyr_.getFeatures():
                g_ = f_.geometry()
                if g_ is not None and not g_.isEmpty():
                    sink_geoms.append(g_)
        sink_ref = nd.SinkRef(sink_geoms)
        feedback.pushInfo(f"   Salidas de referencia (rios/storages en el AOI): {len(sink_ref)} geometria(s).")

        # ---- Drenaje natural (opcional): tramos + conectores a calles ----
        drain_data = None
        if has_drains:
            feedback.pushInfo("   Procesando red de drenaje natural...")
            blocks_fix_tmp = processing.run("native:fixgeometries", {
                "INPUT": blocks_clipped, "OUTPUT": "memory:blocks_fix_drain"
            }, context=context, feedback=feedback)["OUTPUT"]
            blocks_union = nd.union_geometries(blocks_fix_tmp)
            drain_data = nd.prepare_drainage_segments(
                drains_source, crs, context.transformContext(), aoi_geom, blocks_union,
                streets_snapped, sink_ref,
                lambda x, y: self._sample_dem_avg(dem_provider, x, y, cell_size),
                tol, drain_connect_buffer, feedback,
            )
            drain_data["blocks_union"] = blocks_union
            # Los cauces naturales NO entran como conductos: solo aportan por INFLOW (ver mas abajo).

        # Cluster de extremos por proximidad (snap grid) -- compartido entre
        # calles para que las confluencias se fusionen en un nodo.
        key_to_cluster = {}
        cluster_points = []

        def get_cluster(x, y):
            key = self._snap_key(x, y, tol)
            neighbors = [key]
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    neighbors.append((key[0] + dx, key[1] + dy))
            for nk in neighbors:
                if nk in key_to_cluster:
                    cid = key_to_cluster[nk]
                    cluster_points[cid].append((x, y))
                    key_to_cluster[key] = cid
                    return cid
            cid = len(cluster_points)
            cluster_points.append([(x, y)])
            key_to_cluster[key] = cid
            return cid

        conduits_raw = []
        for geom, kind, extra in raw_segments:
            pl = geom.asPolyline() if not geom.isMultipart() else geom.asMultiPolyline()[0]
            if len(pl) < 2:
                continue
            x0, y0 = pl[0].x(), pl[0].y()
            x1, y1 = pl[-1].x(), pl[-1].y()
            c0 = get_cluster(x0, y0)
            c1 = get_cluster(x1, y1)
            if c0 == c1:
                continue  # segmento degenerado / micro-loop, se descarta
            conduits_raw.append({
                "from_cluster": c0, "to_cluster": c1,
                "geom": geom, "kind": kind, "extra": extra, "length": geom.length()
            })

        if not conduits_raw:
            raise QgsProcessingException(
                "No se generaron conductos. Revisar que la capa de ejes de calle "
                "intersecte el area de interes y tenga "
                "geometria valida."
            )

        cluster_coord = []
        for pts in cluster_points:
            mx = sum(p[0] for p in pts) / len(pts)
            my = sum(p[1] for p in pts) / len(pts)
            cluster_coord.append((mx, my))

        degree = defaultdict(int)
        for c in conduits_raw:
            degree[c["from_cluster"]] += 1
            degree[c["to_cluster"]] += 1

        # cabeceras de drenaje natural (extremo aguas arriba sin nada mas): nunca son outfall
        drain_headwater_clusters = {
            c["from_cluster"] for c in conduits_raw
            if c["kind"] == "drain" and degree[c["from_cluster"]] == 1
        }

        # ==============================================================
        # 5. Cotas de nodo (DEM, promedio en cruz) + outfalls (extremos de calle junto a rio/storage)
        # ==============================================================
        feedback.pushInfo("5/10 Muestreando cotas de terreno en nodos...")
        node_ids = {}
        junctions = []
        outfalls = []
        coordinates = {}
        cluster_z = {}        # cluster -> cota de nodo (para orientar conductos por DEM)

        for cid, (x, y) in enumerate(cluster_coord):
            node_id = f"N{cid+1}"
            node_ids[cid] = node_id
            elev = self._sample_dem_avg(dem_provider, x, y, cell_size)
            if elev is None:
                feedback.pushWarning(f"Nodo {node_id}: sin valor de DEM valido, se asigna 0.0")
                elev = 0.0
            coordinates[node_id] = (x, y)
            cluster_z[cid] = elev

            # Lectura B: outfall = extremo de calle (grado 1) a <= outfall_buffer de un rio/storage
            is_outfall = (degree[cid] == 1
                                and cid not in drain_headwater_clusters
                                and sink_ref.distance(QgsPointXY(x, y), outfall_buffer) is not None)

            if is_outfall:
                outfalls.append({"id": node_id, "elev": elev, "type": "FREE", "stage_data": ""})
            else:
                junctions.append({
                    "id": node_id, "elev": elev, "maxdepth": FALLBACK_JUNCTION_MAXDEPTH,
                    "initdepth": 0.0, "surdepth": 0.0, "aponded": 0.0,
                })

        id_to_cid = {nid: cid for cid, nid in node_ids.items()}

        # ==============================================================
        # 6. (Lectura B) Los storages ya NO son nodos/Storage Units: son salidas
        # ==============================================================
        storages = []
        remap = {}

        def resolve_node(nid):
            seen = set()
            while nid in remap and nid not in seen:
                seen.add(nid)
                nid = remap[nid]
            return nid

        # ==============================================================
        # 6b. Asegurar al menos un OUTFALL
        # ==============================================================
        if not outfalls:
            feedback.pushWarning(
                "Ningun nodo de grado 1 queda a <= "
                f"{outfall_buffer:g} m de un rio/storage "
                + ("(no se cargo ningun rio ni storage dentro del AOI). " if len(sink_ref) == 0 else ". ")
                + "Se usa como ultimo recurso el nodo de grado 1 de menor cota de toda la red. "
                "Revisa el buffer de outfalls o la geometria de rios/storages."
            )
            degree1 = [j for j in junctions if degree[id_to_cid.get(j["id"], -1)] == 1]
            if degree1:
                fallback = min(degree1, key=lambda j: j["elev"])
                junctions = [j for j in junctions if j["id"] != fallback["id"]]
                outfalls.append({"id": fallback["id"], "elev": fallback["elev"], "type": "FREE", "stage_data": ""})
            else:
                raise QgsProcessingException(
                    "No se pudo determinar ningun punto de salida valido (OUTFALL) para "
                    "la red: no hay nodos de grado 1 en toda la red (probablemente el area "
                    "de interes forma una malla cerrada sin extremos). SWMM requiere que la "
                    "red drene a algun lado -- opciones: (1) ampliar el area de interes o el "
                    "buffer de outfalls para que una calle termine en un rio/storage, (2) "
                    "revisar que las calles terminen en la orilla (un OUTFALL admite un solo "
                    "conducto: una calle que cruza el rio sin terminar ahi no genera outfall), "
                    "o (3) editar manualmente el .inp para designar un OUTFALL."
                )

        outfall_ids = {o["id"] for o in outfalls}

        # ==============================================================
        # 6c. Orientar los conductos de calle de mayor a menor cota (coherente a nivel de red)
        # ==============================================================
        st_idx = [i for i, c in enumerate(conduits_raw) if c["kind"] == "street"]
        if st_idx:
            outfall_cids = {id_to_cid[nid] for nid in outfall_ids if nid in id_to_cid}
            flips, ost = so.orient_edges(
                [(conduits_raw[i]["from_cluster"], conduits_raw[i]["to_cluster"],
                  conduits_raw[i]["length"]) for i in st_idx],
                cluster_z, outfall_cids)
            n_applied = 0
            for i, do_flip in zip(st_idx, flips):
                if not do_flip:
                    continue
                c = conduits_raw[i]
                g = c["geom"]
                if g.isMultipart():
                    mp = g.asMultiPolyline()
                    if len(mp) != 1:
                        continue      # multiparte real: se deja como esta
                    g = QgsGeometry.fromPolylineXY(list(reversed(mp[0])))
                else:
                    g = QgsGeometry.fromPolylineXY(list(reversed(g.asPolyline())))
                c["geom"] = g
                c["from_cluster"], c["to_cluster"] = c["to_cluster"], c["from_cluster"]
                n_applied += 1
            feedback.pushInfo(
                f"   Orientacion por DEM (red): {n_applied}/{ost['edges']} conducto(s) de calle invertido(s); "
                f"{ost['flat']} tramo(s) plano(s) resueltos por posicion en la red; "
                f"{ost['comp_seeds']} componente(s) sin outfall (sembrados en su nodo mas bajo); "
                f"{ost['adverse_left']} con pendiente adversa residual por ruido/depresion del DEM.")

        # ==============================================================
        # 7. Categorias de calle -> tabla STREETS
        # ==============================================================
        feedback.pushInfo("7/10 Asignando categorias de calle...")
        cat_table_raw = self._load_csv_dict(street_cat_csv, "category")
        street_categories = {}
        for name, row in cat_table_raw.items():
            street_categories[name] = {
                "tcrown": float(row["tcrown_m"]),
                "hcurb": float(row["hcurb_m"]),
                "sx": float(row["sx"]),
                "n": float(row["n"]),
                "a": float(row.get("a", 1.5)),
                "w": float(row["w_m"]),
                "sides": int(row.get("sides", 2)),
                "tback": float(row.get("tback_m", 0.0)),
                "sback": float(row.get("sback", 0.0)),
                "nback": float(row.get("nback", row["n"])),
                "roughness": float(row["n"]),
            }
        default_category = list(street_categories.keys())[0]

        def build_trapezoid_profile(bottom_w, depth, side_slope, bank_w):
            """Perfil sintetico relativo (0 = fondo del canal): fondo plano,
            taludes laterales side_slope (H:V) hasta profundidad `depth`, y
            luego una banca plana `bank_w` a cada lado para contener
            desbordes sin que SWMM falle por falta de recorrido lateral."""
            half_bottom = bottom_w / 2.0
            top_half = half_bottom + side_slope * depth
            pts = [
                (0.0, depth + 1.0),
                (top_half + bank_w, depth + 1.0),
                (top_half, depth),
                (half_bottom, 0.0),
                (-half_bottom, 0.0),
                (-top_half, depth),
                (-(top_half + bank_w), depth + 1.0),
            ]
            # reordenar como estaciones crecientes de izquierda a derecha
            pts = sorted(pts, key=lambda p: p[0])
            # normalizar estaciones a partir de 0
            min_sta = pts[0][0]
            return [(round(sta - min_sta, 3), round(elev, 3)) for sta, elev in pts]

        def load_profile_csv(path):
            pts = []
            with open(path, newline="", encoding="utf-8-sig") as pf:
                reader = csv.reader(pf)
                for row in reader:
                    if not row or row[0].strip().lower() in ("station", "sta", ""):
                        continue
                    try:
                        sta = float(row[0]); elev = float(row[1])
                    except (ValueError, IndexError):
                        continue
                    pts.append((sta, elev))
            return sorted(pts, key=lambda p: p[0])

        transects = {}
        conduits = []
        conduit_categories_by_node = defaultdict(list)  # node_id -> [hcurb, ...] para MaxDepth
        conduit_kinds_by_node = defaultdict(set)
        drain_depths_by_node = defaultdict(list)   # drenaje natural: profundidad de cauce por nodo
        drain_outlet_by_seg = {}                   # seg_id de drenaje -> nodo aguas abajo
        drain_street_node_by_seg = {}              # seg_id de drenaje -> nodo de calle vinculado por buffer

        for i, c in enumerate(conduits_raw):
            from_id = resolve_node(node_ids[c["from_cluster"]])
            to_id = resolve_node(node_ids[c["to_cluster"]])
            if from_id == to_id:
                continue  # colapsado por fusion con storage

            cid_name = f"C{i+1}"
            if c["kind"] == "street":
                cat_val = str(c["extra"]["category"]) if c["extra"]["category"] not in (None, "") else default_category
                if cat_val not in street_categories:
                    feedback.pushWarning(
                        f"Categoria de calle '{cat_val}' no encontrada en la tabla, "
                        f"se usa '{default_category}' por defecto."
                    )
                    cat_val = default_category
                conduits.append({
                    "id": cid_name, "from_node": from_id, "to_node": to_id,
                    "length": round(c["length"], 3),
                    "roughness": street_categories[cat_val]["roughness"],
                    "xsection_type": "STREET", "category": cat_val,
                    "geom": c["geom"], "kind": "street",
                })
                conduit_categories_by_node[from_id].append(street_categories[cat_val]["hcurb"])
                conduit_categories_by_node[to_id].append(street_categories[cat_val]["hcurb"])
                conduit_kinds_by_node[from_id].add("street")
                conduit_kinds_by_node[to_id].add("street")
            elif c["kind"] == "drain":
                extra = c["extra"]
                tp = extra["type_params"]
                transect_name = "TRD_" + str(extra["type_name"]).replace(" ", "_")
                if transect_name not in transects:
                    transects[transect_name] = {
                        "n_left": tp["manning_n"], "n_right": tp["manning_n"],
                        "n_channel": tp["manning_n"],
                        "station_elev": build_trapezoid_profile(
                            tp["bottom_w"], tp["depth"], tp["side_slope"], tp["bank_w"]),
                    }
                conduits.append({
                    "id": cid_name, "from_node": from_id, "to_node": to_id,
                    "length": round(c["length"], 3),
                    "roughness": tp["manning_n"],
                    "xsection_type": "IRREGULAR", "transect": transect_name,
                    "geom": c["geom"], "kind": "drain",
                })
                for nid in (from_id, to_id):
                    conduit_kinds_by_node[nid].add("drain")
                    drain_depths_by_node[nid].append(tp["depth"])
                if extra.get("seg_id") is not None:
                    drain_outlet_by_seg[extra["seg_id"]] = to_id
                for sid in (extra.get("link_of") or []):
                    drain_street_node_by_seg[sid] = to_id
            else:
                extra = c["extra"]
                profile_path = extra.get("profile_csv")
                if profile_path and os.path.isfile(profile_path):
                    station_elev = load_profile_csv(profile_path)
                else:
                    station_elev = build_trapezoid_profile(
                        extra["bottom_w"], extra["depth_m"], extra["side_slope"], extra["bank_w"]
                    )
                transect_name = extra.get("name") or f"TR_{cid_name}"
                transect_name = str(transect_name)
                transects[transect_name] = {
                    "n_left": extra["manning_n"], "n_right": extra["manning_n"],
                    "n_channel": extra["manning_n"], "station_elev": station_elev,
                }
                conduits.append({
                    "id": cid_name, "from_node": from_id, "to_node": to_id,
                    "length": round(c["length"], 3),
                    "roughness": extra["manning_n"],
                    "xsection_type": "IRREGULAR", "transect": transect_name,
                    "geom": c["geom"], "kind": "river",
                })
                conduit_kinds_by_node[from_id].add("river")
                conduit_kinds_by_node[to_id].add("river")

        # MaxDepth: minimo Hcurb entre las calles conectadas a cada nodo. Si
        # el nodo solo conecta con rio (sin cordon fisico), se usa un valor
        # generico mas amplio para no limitar artificialmente la profundidad.
        for j in junctions:
            hcurbs = conduit_categories_by_node.get(j["id"])
            if hcurbs:
                j["maxdepth"] = min(hcurbs)
            elif "river" in conduit_kinds_by_node.get(j["id"], set()):
                j["maxdepth"] = RIVER_ONLY_NODE_MAXDEPTH
            elif "drain" in conduit_kinds_by_node.get(j["id"], set()):
                j["maxdepth"] = max(drain_depths_by_node[j["id"]]) + 1.0

        # ==============================================================
        # 8. Subcuencas a partir de manzanas (subdivision si son grandes)
        # ==============================================================
        feedback.pushInfo("8/10 Generando subcuencas desde manzanas...")
        blocks_final = self._prepare_subcatchment_polygons(
            blocks_clipped, max_block_area_ha, context, feedback
        )

        imperv_raster = self._reclassify_landuse(landuse_layer, landuse_csv, context, feedback)

        slope_raster_path = processing.run("native:slope", {
            "INPUT": dem_layer, "Z_FACTOR": 1, "OUTPUT": "TEMPORARY_OUTPUT"
        }, context=context, feedback=feedback)["OUTPUT"]

        blocks_with_imperv = processing.run("native:zonalstatisticsfb", {
            "INPUT": blocks_final, "INPUT_RASTER": imperv_raster, "RASTER_BAND": 1,
            "COLUMN_PREFIX": "imp_", "STATISTICS": [2], "OUTPUT": "TEMPORARY_OUTPUT"
        }, context=context, feedback=feedback)["OUTPUT"]

        blocks_with_slope = processing.run("native:zonalstatisticsfb", {
            "INPUT": blocks_with_imperv, "INPUT_RASTER": slope_raster_path, "RASTER_BAND": 1,
            "COLUMN_PREFIX": "slp_", "STATISTICS": [2], "OUTPUT": "TEMPORARY_OUTPUT"
        }, context=context, feedback=feedback)["OUTPUT"]
        # Salidas escritas a disco (TEMPORARY_OUTPUT) a proposito: encadenar
        # varias capas "memory:" como INPUT/OUTPUT de sucesivos
        # processing.run() puede hacer que QGIS destruya el objeto C++
        # intermedio ("wrapped C/C++ object ... has been deleted") con
        # algoritmos in-place como zonalstatisticsfb.
        if isinstance(blocks_with_slope, str):
            blocks_with_slope = QgsVectorLayer(blocks_with_slope, "blocks_with_slope", "ogr")

        # indice espacial de nodos (para Outlet de subcuenca por vecino mas cercano)
        node_lookup = {}
        idx = QgsSpatialIndex()
        tmp_id = 1
        for nid, (x, y) in coordinates.items():
            if nid in remap:
                continue  # nodos absorbidos por un storage no son destino valido
            if conduit_kinds_by_node.get(nid) == {"drain"}:
                continue  # nodos solo de drenaje natural: las manzanas no descargan ahi
            f = QgsFeature()
            f.setId(tmp_id)
            f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x, y)))
            idx.insertFeature(f)
            node_lookup[tmp_id] = nid
            tmp_id += 1

        # indice espacial de pluviometros (si hay) para asignar Gage por subcuenca
        raingage_idx = None
        raingage_lookup = {}
        raingages = []
        stations = []          # pluviometros validos para Thiessen
        timeseries = {}

        if has_raingages:
            feedback.pushInfo("   Cargando tabla de pluviometros...")
            rg_table = self._load_csv_dict(raingages_csv, "codigo")
            rg_csv_dir = os.path.dirname(raingages_csv)
            raingage_idx = QgsSpatialIndex()
            rg_tmp_id = 1
            fmt_map = {"intensity": "INTENSITY", "volume": "VOLUME", "cumulative": "CUMULATIVE"}
            for feat in raingages_source.getFeatures():
                code = str(feat[raingage_code_field])
                if code not in rg_table:
                    feedback.pushWarning(
                        f"Pluviometro con codigo '{code}' no encontrado en la tabla CSV, se omite."
                    )
                    continue
                row = rg_table[code]
                rain_format = fmt_map.get(str(row.get("rain_format", "intensity")).strip().lower(), "INTENSITY")
                file_path = row["file"]
                if not os.path.isabs(file_path):
                    file_path = os.path.join(rg_csv_dir, file_path)
                if not os.path.isfile(file_path):
                    feedback.pushWarning(
                        f"Pluviometro '{code}': no se encontro el archivo de hietograma "
                        f"'{file_path}', se omite."
                    )
                    continue
                gage_name = f"RG_{code}"
                ts_name = f"TS_{code}"
                raingages.append({
                    "name": gage_name, "format": rain_format, "interval": "0:05",
                    "scf": 1.0, "source": ts_name,
                })
                timeseries[ts_name] = self._load_rainfall_csv(file_path)
                hy_rg = nd.parse_hyetograph(timeseries[ts_name], rain_format)
                if hy_rg is not None:
                    raingages[-1]["interval"] = th.interval_str(hy_rg["dt_h"])

                geom = feat.geometry()
                if geom is None or geom.isEmpty():
                    continue
                pt = geom.asPoint()
                gf = QgsFeature()
                gf.setId(rg_tmp_id)
                gf.setGeometry(QgsGeometry.fromPointXY(pt))
                raingage_idx.insertFeature(gf)
                raingage_lookup[rg_tmp_id] = gage_name
                rg_tmp_id += 1
                if hy_rg is not None:
                    stations.append({"name": gage_name, "pt": QgsPointXY(pt), "hy": hy_rg})

            if not raingages:
                raise QgsProcessingException(
                    "Se cargo la capa/tabla de pluviometros pero ningun registro pudo "
                    "asociarse correctamente (revisar codigos y rutas de archivo)."
                )

        if not raingages:
            # fallback: pluviometro unico, comportamiento previo del plugin
            if not rainfall_csv:
                raise QgsProcessingException(
                    "No se cargaron pluviometros y tampoco se indico un hietograma por "
                    "defecto (Hietograma por defecto)."
                )
            raingages = [{
                "name": "RG1", "format": "INTENSITY", "interval": "0:05",
                "scf": 1.0, "source": "TS_RAIN",
            }]
            timeseries["TS_RAIN"] = self._load_rainfall_csv(rainfall_csv)
            hy_def = nd.parse_hyetograph(timeseries["TS_RAIN"], "INTENSITY")
            if hy_def is not None:
                raingages[0]["interval"] = th.interval_str(hy_def["dt_h"])

        default_gage_name = raingages[0]["name"]
        assigner = None
        if has_raingages and stations:
            assigner = th.RainAssigner(stations, aoi_geom, thiessen_step, feedback)

        subcatchments = []
        polygons_coords = {}
        sub_count = 0
        for feat in blocks_with_slope.getFeatures():
            geom = feat.geometry()
            if geom is None or geom.isEmpty():
                continue
            area_m2 = geom.area()
            if area_m2 <= 1.0:
                continue
            perimeter_m = geom.length()
            sub_count += 1
            sub_id = f"S{sub_count}"
            centroid = geom.centroid().asPoint()
            nearest_ids = idx.nearestNeighbor(centroid, 1)
            if not nearest_ids:
                feedback.pushWarning(f"{sub_id}: no se encontro nodo cercano, se omite.")
                sub_count -= 1
                continue
            outlet_id = node_lookup[nearest_ids[0]]

            if assigner is not None:
                gage_name, pesos = assigner.gage_for_geometry(geom)     # Thiessen
            else:
                gage_name, pesos = default_gage_name, ""

            pct_imperv = feat["imp_mean"] if feat["imp_mean"] is not None else 50.0
            pct_slope = feat["slp_mean"] if feat["slp_mean"] is not None else 1.0
            pct_slope = max(pct_slope, 0.1)

            width_m = area_m2 / max(width_coef * perimeter_m, 1.0)

            subcatchments.append({
                "id": sub_id,
                "raingage": gage_name,
                "pesos": pesos,
                "outlet": outlet_id,
                "area_ha": area_m2 / 10000.0,
                "pct_imperv": round(pct_imperv, 1),
                "width_m": round(width_m, 1),
                "pct_slope": round(pct_slope, 2),
                "n_imp": 0.013,
                "n_perv": 0.24,
                "dstore_imp_mm": 2.0,
                "dstore_perv_mm": 5.0,
                "pct_zero": 25.0,
                "infil_params": {"f0": 75, "fmin": 13, "decay": 4, "drytime": 7, "maxinfil": 0},
                "geom": geom,
            })

            if geom.isMultipart():
                ring = geom.asMultiPolygon()[0][0]
            else:
                ring = geom.asPolygon()[0]
            polygons_coords[sub_id] = [(p.x(), p.y()) for p in ring]

        # ---- Pluviometros sinteticos (Thiessen) para subcuencas urbanas con varias estaciones ----
        extra_ts = {}
        if assigner is not None:
            syn_gages, extra_ts = assigner.synthetic_gages()
            raingages.extend(syn_gages)
            if syn_gages:
                feedback.pushInfo(f"  -> {len(syn_gages)} pluviometro(s) sintetico(s) Thiessen (RG_TH#) para subcuencas entre estaciones.")

        # ---- Aporte de drenaje natural (opcional): INFLOW + HU sintetico SCS ----
        natural_inflows, natural_ts = [], {}
        basin_records, network_rows = [], []
        if drain_data is not None:
            feedback.pushInfo("   Delineando cuencas de aporte y calculando hidrogramas SCS...")
            # nodos destino validos: que no sean OUTFALL
            in_idx = QgsSpatialIndex()
            in_lookup = {}
            k_in = 1
            for nid, (x, y) in coordinates.items():
                if nid in remap or nid in outfall_ids:
                    continue
                f = QgsFeature()
                f.setId(k_in)
                f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x, y)))
                in_idx.insertFeature(f)
                in_lookup[k_in] = nid
                k_in += 1

            usable = []
            for og in drain_data["outlets"]:
                node_id = None
                for cand in (in_idx.nearestNeighbor(og["point"], 3) if in_lookup else []):
                    nid = in_lookup[cand]
                    cx, cy = coordinates[nid]
                    if QgsPointXY(cx, cy).distance(og["point"]) <= drain_connect_buffer + tol:
                        node_id = nid
                        break
                if node_id is None:
                    feedback.pushWarning(
                        f"Outlet de drenaje natural en ({og['point'].x():.1f}, {og['point'].y():.1f}): "
                        f"no hay ningun nodo (no OUTFALL) a <= {drain_connect_buffer} m; se omite.")
                    continue
                usable.append({**og, "node": node_id})

            cn_table = nd.load_cn_table(landuse_csv)
            if not cn_table:
                feedback.pushWarning(
                    "La tabla de cobertura no tiene la columna 'cn': se usa CN "
                    f"{nd.DEFAULT_CN} en todas las cuencas naturales.")

            hy_default = nd.parse_hyetograph(
                timeseries.get(raingages[0]["source"]), raingages[0].get("format", "INTENSITY"))

            def hyeto_for(xs, ys):
                if assigner is not None:
                    return assigner.hyetograph_for_cells(xs, ys)         # Thiessen
                return default_gage_name, hy_default

            basin_records = nd.delineate_natural_basins(
                dem_layer, landuse_layer, cn_table, drain_data["segments"], usable,
                aoi_geom, drain_data["blocks_union"], hyeto_for, feedback, ia_ratio)
            natural_inflows, natural_ts = nd.assemble_node_inflows(basin_records, feedback)
            feedback.pushInfo(f"  -> {len(natural_inflows)} INFLOW(s) de drenaje natural en nodos de calle.")
            for r in basin_records:
                dur_h = float(r["t_h"][-1] - r["t_h"][0])
                if dur_h > 6.0:
                    feedback.pushWarning(
                        f"El hidrograma del nodo {r['node']} dura {dur_h:.1f} h; el END_TIME por defecto del "
                        f".inp es 06:00:00 -- ajustalo en Simulation Options de SWMM.")
                    break

            for sg in drain_data["segments"]:
                network_rows.append((sg["geom"], "drenaje", sg["seg_id"], "", sg["length"]))
            for sg in drain_data["purged"]:
                network_rows.append((sg["geom"], "purgado", sg["seg_id"], "", sg["length"]))
            for sg in drain_data.get("to_sink", []):
                network_rows.append((sg["geom"], "a_salida", sg["seg_id"], "", sg["length"]))
            for u in drain_data.get("unlinked", []):
                if u["nearest"] is not None:
                    network_rows.append((QgsGeometry.fromPolylineXY([u["point"], u["nearest"]]),
                                         "sin_conexion", None, "", u["dist"]))
            for og in usable:
                if og["link"] is not None:
                    network_rows.append((og["link"], "conector", None, og["node"], og["link"].length()))

        # ==============================================================
        # 9. Ensamblar .inp y capas de salida
        # ==============================================================
        feedback.pushInfo("9/10 Escribiendo archivo .inp...")
        inp_text = build_inp(
            title=model_title,
            junctions=junctions,
            outfalls=outfalls,
            conduits=[{k: v for k, v in c.items() if k != "geom"} for c in conduits],
            street_categories=street_categories,
            subcatchments=[{k: v for k, v in s.items() if k != "geom"} for s in subcatchments],
            raingages=raingages,
            timeseries=timeseries,
            coordinates=coordinates,
            polygons=polygons_coords,
            storages=storages,
            transects=transects,
            inflows=natural_inflows,
            inflow_timeseries=natural_ts,
            extra_timeseries=extra_ts,
        )
        with open(inp_path, "w", encoding="utf-8") as f:
            f.write(inp_text)

        feedback.pushInfo("10/10 Escribiendo capas de salida...")
        nodes_dest = self._write_nodes_layer(parameters, context, junctions, outfalls, storages, coordinates, crs)
        conduits_dest = self._write_conduits_layer(parameters, context, conduits, crs)
        subs_dest = self._write_subcatchments_layer(parameters, context, subcatchments, crs)
        thiessen_dest = self._write_thiessen_layer(parameters, context, assigner, aoi_geom, crs)
        basins_dest = self._write_natural_basins_layer(parameters, context, basin_records, crs)
        network_dest = self._write_natural_network_layer(parameters, context, network_rows, crs)

        feedback.pushInfo(
            f"Listo: {len(junctions)} junctions, {len(outfalls)} outfalls, "
            f"{len(conduits)} conductos, "
            f"{len(subcatchments)} subcuencas, {len(raingages)} pluviometro(s)."
        )

        return {
            self.OUTPUT_INP: inp_path,
            self.OUTPUT_NODES: nodes_dest,
            self.OUTPUT_CONDUITS: conduits_dest,
            self.OUTPUT_SUBCATCHMENTS: subs_dest,
            **({self.OUTPUT_THIESSEN: thiessen_dest} if thiessen_dest else {}),
            **({self.OUTPUT_NATURAL_BASINS: basins_dest} if basins_dest else {}),
            **({self.OUTPUT_NATURAL_NETWORK: network_dest} if network_dest else {}),
        }

    # ------------------------------------------------------------------
    @staticmethod
    def _generate_grid_cells(bbox, hspacing, vspacing):
        """Genera geometrias de celdas rectangulares cubriendo bbox, sin
        pasar por processing.run (ver nota en _prepare_subcatchment_polygons)."""
        cells = []
        x = bbox.xMinimum()
        while x < bbox.xMaximum():
            x2 = min(x + hspacing, bbox.xMaximum())
            y = bbox.yMinimum()
            while y < bbox.yMaximum():
                y2 = min(y + vspacing, bbox.yMaximum())
                ring = [
                    QgsPointXY(x, y), QgsPointXY(x2, y),
                    QgsPointXY(x2, y2), QgsPointXY(x, y2), QgsPointXY(x, y),
                ]
                cells.append(QgsGeometry.fromPolygonXY([ring]))
                y += vspacing
            x += hspacing
        return cells

    # ------------------------------------------------------------------
    def _prepare_subcatchment_polygons(self, blocks_layer, max_area_ha, context, feedback):
        """Subdivide manzanas mayores a max_area_ha con una grilla regular."""
        max_area_m2 = max_area_ha * 10000.0

        big_ok_layer = processing.run("native:fixgeometries", {
            "INPUT": blocks_layer, "OUTPUT": "memory:blocks_fix"
        }, context=context, feedback=feedback)["OUTPUT"]

        small_feats = []
        big_feats = []
        for feat in big_ok_layer.getFeatures():
            geom = feat.geometry()
            if geom is None or geom.isEmpty():
                continue
            if geom.area() > max_area_m2:
                big_feats.append(feat)
            else:
                small_feats.append(feat)

        if not big_feats:
            return big_ok_layer

        feedback.pushInfo(
            f"  -> {len(big_feats)} manzana(s) superan {max_area_ha} ha, subdividiendo con grilla..."
        )

        crs_auth = blocks_layer.sourceCrs().authid()
        out_layer = QgsVectorLayer(f"Polygon?crs={crs_auth}", "blocks_subdivided", "memory")
        out_layer.dataProvider().addAttributes(blocks_layer.fields())
        out_layer.updateFields()

        feats_out = []
        for feat in small_feats:
            nf = QgsFeature(out_layer.fields())
            nf.setGeometry(feat.geometry())
            nf.setAttributes(feat.attributes())
            feats_out.append(nf)

        for feat in big_feats:
            geom = feat.geometry()
            bbox = geom.boundingBox()
            side = math.sqrt(max_area_m2)
            hspacing = min(side, max(bbox.width() / 2.0, 0.01))
            vspacing = min(side, max(bbox.height() / 2.0, 0.01))
            # Grilla generada directamente con geometria de PyQGIS (sin pasar
            # por processing.run/native:creategrid): encadenar cientos de
            # capas "memory:" en un bucle de 375 iteraciones dispara el
            # mismo problema de estabilidad ya visto ("wrapped C/C++ object
            # ... has been deleted"), asi que se evita por completo.
            for cell_geom in self._generate_grid_cells(bbox, hspacing, vspacing):
                clipped = geom.intersection(cell_geom)
                if clipped is None or clipped.isEmpty() or clipped.area() < 1.0:
                    continue
                nf = QgsFeature(out_layer.fields())
                nf.setGeometry(clipped)
                nf.setAttributes(feat.attributes())
                feats_out.append(nf)

        out_layer.dataProvider().addFeatures(feats_out)
        out_layer.updateExtents()
        return out_layer

    # ------------------------------------------------------------------
    def _reclassify_landuse(self, landuse_layer, csv_path, context, feedback):
        table = []
        with open(csv_path, newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            for row in reader:
                v = float(row["class_value"])
                pct = float(row["pct_imperv"])
                table.append([v - 0.5, v + 0.5, pct])

        flat_table = [x for triplet in table for x in triplet]
        result = processing.run("native:reclassifybytable", {
            "INPUT_RASTER": landuse_layer, "RASTER_BAND": 1,
            "TABLE": flat_table, "NO_DATA": -9999, "RANGE_BOUNDARIES": 0,
            "NODATA_FOR_MISSING": True, "DATA_TYPE": 5,
            "OUTPUT": "TEMPORARY_OUTPUT"
        }, context=context, feedback=feedback)["OUTPUT"]
        return result

    # ------------------------------------------------------------------
    @staticmethod
    def _load_rainfall_csv(path):
        values = []
        with open(path, newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            fieldnames = [fn.lower() for fn in reader.fieldnames]
            has_date = "date" in fieldnames
            for row in reader:
                row = {k.lower(): v for k, v in row.items()}
                date_str = row.get("date", "") if has_date else ""
                time_str = row["time"]
                value = float(row["value"])
                values.append((date_str, time_str, value))
        return values

    # ------------------------------------------------------------------
    def _write_nodes_layer(self, parameters, context, junctions, outfalls, storages, coordinates, crs):
        fields = QgsFields()
        fields.append(QgsField("id", QVariant.String))
        fields.append(QgsField("node_type", QVariant.String))
        fields.append(QgsField("elevation", QVariant.Double))

        sink, dest_id = self.parameterAsSink(
            parameters, self.OUTPUT_NODES, context, fields, QgsWkbTypes.Point, crs
        )
        for j in junctions:
            f = QgsFeature(fields)
            x, y = coordinates[j["id"]]
            f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x, y)))
            f.setAttributes([j["id"], "JUNCTION", j["elev"]])
            sink.addFeature(f, QgsFeatureSink.FastInsert)
        for o in outfalls:
            f = QgsFeature(fields)
            x, y = coordinates[o["id"]]
            f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x, y)))
            f.setAttributes([o["id"], "OUTFALL", o["elev"]])
            sink.addFeature(f, QgsFeatureSink.FastInsert)
        for s in storages:
            f = QgsFeature(fields)
            x, y = coordinates[s["id"]]
            f.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x, y)))
            f.setAttributes([s["id"], "STORAGE", s["elev"]])
            sink.addFeature(f, QgsFeatureSink.FastInsert)
        return dest_id

    # ------------------------------------------------------------------
    def _write_conduits_layer(self, parameters, context, conduits, crs):
        fields = QgsFields()
        fields.append(QgsField("id", QVariant.String))
        fields.append(QgsField("from_node", QVariant.String))
        fields.append(QgsField("to_node", QVariant.String))
        fields.append(QgsField("length_m", QVariant.Double))
        fields.append(QgsField("kind", QVariant.String))
        fields.append(QgsField("category", QVariant.String))
        fields.append(QgsField("roughness", QVariant.Double))

        sink, dest_id = self.parameterAsSink(
            parameters, self.OUTPUT_CONDUITS, context, fields, QgsWkbTypes.LineString, crs
        )
        for c in conduits:
            f = QgsFeature(fields)
            f.setGeometry(c["geom"])
            cat_display = c["category"] if c["kind"] == "street" else c["transect"]
            f.setAttributes([c["id"], c["from_node"], c["to_node"], c["length"], c["kind"], cat_display, c["roughness"]])
            sink.addFeature(f, QgsFeatureSink.FastInsert)
        return dest_id

    # ------------------------------------------------------------------
    def _write_thiessen_layer(self, parameters, context, assigner, aoi_geom, crs):
        """Capa opcional: poligonos de Thiessen de los pluviometros (recortados al AOI)."""
        if assigner is None:
            return None
        fields = QgsFields()
        fields.append(QgsField("raingage", QVariant.String))
        fields.append(QgsField("area_ha", QVariant.Double))
        sink, dest_id = self.parameterAsSink(
            parameters, self.OUTPUT_THIESSEN, context, fields, QgsWkbTypes.MultiPolygon, crs
        )
        if sink is None:
            return None
        for name, geom, area_ha in assigner.cell_polygons(aoi_geom):
            geom.convertToMultiType()
            f = QgsFeature(fields)
            f.setGeometry(geom)
            f.setAttributes([name, round(area_ha, 3)])
            sink.addFeature(f, QgsFeatureSink.FastInsert)
        return dest_id

    def _write_natural_basins_layer(self, parameters, context, records, crs):
        """Capa (shp/gpkg) opcional de cuencas de aporte de drenaje natural (solo verificacion)."""
        fields = QgsFields()
        for name, qtype in (("node", QVariant.String), ("gage", QVariant.String),
                            ("area_ha", QVariant.Double), ("L_m", QVariant.Double),
                            ("H_m", QVariant.Double), ("Tc_min", QVariant.Double),
                            ("CN", QVariant.Double), ("P_mm", QVariant.Double),
                            ("Pe_mm", QVariant.Double), ("dt_min", QVariant.Double),
                            ("Qp_unit", QVariant.Double), ("Qmax_cms", QVariant.Double),
                            ("inflow", QVariant.Int)):
            fields.append(QgsField(name, qtype))
        sink, dest_id = self.parameterAsSink(
            parameters, self.OUTPUT_NATURAL_BASINS, context, fields, QgsWkbTypes.Polygon, crs
        )
        if sink is None:
            return None
        for r in records:
            f = QgsFeature(fields)
            f.setGeometry(r["geom"])
            f.setAttributes([r["node"], r["gage"], round(r["area_ha"], 3), round(r["L_m"], 1),
                             round(r["H_m"], 2), round(r["tc_min"], 2), round(r["cn"], 1),
                             round(r["p_mm"], 2), round(r["pe_mm"], 3), round(r["dt_h"] * 60.0, 3),
                             round(r["qp_unit"], 4), round(r["qmax"], 4), 1 if r.get("inflow", True) else 0])
            sink.addFeature(f, QgsFeatureSink.FastInsert)
        return dest_id

    def _write_natural_network_layer(self, parameters, context, rows, crs):
        """Capa (shp/gpkg) opcional: drenaje conservado, purgado y conectores outlet -> nodo."""
        fields = QgsFields()
        for name, qtype in (("tipo", QVariant.String), ("seg_id", QVariant.Int),
                            ("node", QVariant.String), ("length_m", QVariant.Double)):
            fields.append(QgsField(name, qtype))
        sink, dest_id = self.parameterAsSink(
            parameters, self.OUTPUT_NATURAL_NETWORK, context, fields, QgsWkbTypes.LineString, crs
        )
        if sink is None:
            return None
        for geom, tipo, seg_id, node, length in rows:
            f = QgsFeature(fields)
            f.setGeometry(geom)
            f.setAttributes([tipo, seg_id, node, round(length, 2)])
            sink.addFeature(f, QgsFeatureSink.FastInsert)
        return dest_id

    # ------------------------------------------------------------------
    def _write_subcatchments_layer(self, parameters, context, subcatchments, crs):
        fields = QgsFields()
        fields.append(QgsField("id", QVariant.String))
        fields.append(QgsField("outlet", QVariant.String))
        fields.append(QgsField("raingage", QVariant.String))
        fields.append(QgsField("pesos", QVariant.String))
        fields.append(QgsField("area_ha", QVariant.Double))
        fields.append(QgsField("pct_imperv", QVariant.Double))
        fields.append(QgsField("width_m", QVariant.Double))
        fields.append(QgsField("pct_slope", QVariant.Double))

        sink, dest_id = self.parameterAsSink(
            parameters, self.OUTPUT_SUBCATCHMENTS, context, fields, QgsWkbTypes.Polygon, crs
        )
        for s in subcatchments:
            f = QgsFeature(fields)
            f.setGeometry(s["geom"])
            f.setAttributes([s["id"], s["outlet"], s["raingage"], s.get("pesos", ""), round(s["area_ha"], 4),
                              s["pct_imperv"], s["width_m"], s["pct_slope"]])
            sink.addFeature(f, QgsFeatureSink.FastInsert)
        return dest_id
