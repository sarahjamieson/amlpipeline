from illuminate import InteropTileMetrics, InteropControlMetrics, InteropErrorMetrics, InteropExtractionMetrics, \
    InteropIndexMetrics, InteropQualityMetrics, InteropCorrectedIntensityMetrics
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import ScalarFormatter
import plotly.plotly as py
import plotly.tools as tls
import os
import datetime
from pandas import ExcelWriter
InterOp = '/media/sf_sarah_share/InterOp/'

tilemetrics = InteropTileMetrics('%sTileMetricsOut.bin' % InterOp)
print tilemetrics

read_1_phas = tilemetrics.mean_phasing[0] * 100
read_2_phas = tilemetrics.mean_phasing[1] * 100
read_1_prephas = tilemetrics.mean_prephasing[0] * 100
read_2_prephas = tilemetrics.mean_prephasing[1] * 100
