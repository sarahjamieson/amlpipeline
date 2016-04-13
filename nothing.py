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

tilemetrics = InteropTileMetrics('/media/sf_sarah_share/AML_data/InterOp/TileMetricsOut.bin')
indexmetrics = InteropIndexMetrics('/media/sf_sarah_share/AML_data/InterOp/IndexMetricsOut.bin')

total_aligned_clusters = float(tilemetrics.num_clusters * tilemetrics.aligned)

pf_aligned_clusters = float(indexmetrics.df['clusters'].sum())

percent_pf = format(float(pf_aligned_clusters / total_aligned_clusters) * 100, '.4f')

print percent_pf

print indexmetrics




