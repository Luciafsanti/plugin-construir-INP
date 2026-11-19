# -*- coding: utf-8 -*-
"""
SWMM Street Builder
--------------------
Plugin de QGIS (Processing) que genera un modelo EPA SWMM 5.2 de escurrimiento
superficial en calles (sistema mayor) a partir de:
    - DEM (raster)
    - Poligono de extent / area de interes
    - Ejes de calle (vector de lineas)
    - Manzanas (vector de poligonos)
    - Uso de suelo clasificado (raster)

Salidas:
    - Archivo .inp de SWMM
    - Capas GIS equivalentes (nodos, conductos, subcuencas) en GeoPackage
"""


def classFactory(iface):
    from .swmm_street_builder_plugin import SwmmStreetBuilderPlugin
    return SwmmStreetBuilderPlugin(iface)
