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

tile_df = tilemetrics.df[tilemetrics.df['code'] == 100]
cols_to_drop = ['code', 'lane']
tile_df = tile_df.drop(cols_to_drop, axis=1)
tile_df = tile_df.sort_values('tile')

fig = plt.figure()
ax = fig.add_subplot(111)
my_cmap = plt.cm.get_cmap('RdYlGn')

num_tiles = tilemetrics.num_tiles
tiles_11 = tile_df['value'][0:num_tiles/2].tolist()
tiles_21 = tile_df['value'][num_tiles/2:num_tiles].tolist()
myint = 1000
tiles_11_div = [x / myint for x in tiles_11]
tiles_21_div = [x / myint for x in tiles_21]
max_clusters = max(tile_df['value']) / 1000
min_clusters = min(tile_df['value']) / 1000

data = []
item = 0
while item < (num_tiles/2):
    data.append([tiles_11_div[item], tiles_21_div[item]])
    item += 1

heatmap = ax.pcolor(data, cmap=my_cmap, vmin=min_clusters, vmax=max_clusters)

plt.locator_params(axis='x', nbins=4)  # check these
plt.locator_params(axis='y', nbins=28)

x_ticks = [1, 2]
x_labels = ['11', '21']
y_labels = []
y_ticks = []
y = 1
while y <= (num_tiles/2):
    if y < 10:
        y_labels.append('0%s' % str(y))
    else:
        y_labels.append(str(y))
    y_ticks.append(float(y))
    y += 1

y_ticks_cen = [x - 0.5 for x in y_ticks]
x_ticks_cen = [x - 0.5 for x in x_ticks]

ax.set_xticks(x_ticks_cen, minor=False)
ax.set_xticklabels(x_labels, ha='center')

ax.set_yticks(y_ticks_cen, minor=False)
ax.set_yticklabels(y_labels, minor=False)

plt.colorbar(heatmap)
plt.show()




