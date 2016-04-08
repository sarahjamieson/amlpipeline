from illuminate import InteropTileMetrics, InteropControlMetrics, InteropErrorMetrics, InteropExtractionMetrics, \
    InteropIndexMetrics, InteropQualityMetrics, InteropCorrectedIntensityMetrics
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import ScalarFormatter
import plotly.plotly as py
import plotly.tools as tls

indexmetrics = InteropIndexMetrics('/media/sf_sarah_share/AML_data/InterOp/IndexMetricsOut.bin')
index_df = indexmetrics.df


cols_to_drop = ['lane', 'project_str', 'tile', 'read']
index_df = index_df.drop(cols_to_drop, axis=1)
grouped_by_index = index_df.groupby(['name_str', 'index_str'])
x = grouped_by_index.aggregate(np.sum)
print x
x = x.sort_values('name_str', ascending=True)

# sorting dataframe by index (patient) to find number of reads/clusters per index to see if balanced.


