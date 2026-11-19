# -*- coding: utf-8 -*-
from qgis.core import QgsApplication
from .swmm_provider import SwmmStreetBuilderProvider


class SwmmStreetBuilderPlugin:
    """Registra el proveedor de Processing con el algoritmo SWMM Street Builder."""

    def __init__(self, iface):
        self.iface = iface
        self.provider = None

    def initGui(self):
        self.provider = SwmmStreetBuilderProvider()
        QgsApplication.processingRegistry().addProvider(self.provider)

    def unload(self):
        QgsApplication.processingRegistry().removeProvider(self.provider)
