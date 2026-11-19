# -*- coding: utf-8 -*-
from qgis.core import QgsProcessingProvider
from qgis.PyQt.QtGui import QIcon
import os

from .swmm_algorithm import SwmmStreetBuilderAlgorithm


class SwmmStreetBuilderProvider(QgsProcessingProvider):

    def id(self):
        return "swmm_street_builder"

    def name(self):
        return "SWMM Street Builder"

    def icon(self):
        path = os.path.join(os.path.dirname(__file__), "icons", "icon.png")
        if os.path.exists(path):
            return QIcon(path)
        return super().icon()

    def loadAlgorithms(self):
        self.addAlgorithm(SwmmStreetBuilderAlgorithm())
