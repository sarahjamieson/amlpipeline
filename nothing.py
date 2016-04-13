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


def calculate_percentage(value):
    tilemetrics = InteropTileMetrics('/media/sf_sarah_share/AML_data/InterOp/TileMetricsOut.bin')

    return (value / tilemetrics.num_clusters_pf) * 100


def new_function():
    qualitymetrics = InteropQualityMetrics('/media/sf_sarah_share/AML_data/InterOp/QMetricsOut.bin')

    quality_df = qualitymetrics.df
    cols_to_drop = ['lane', 'tile']
    quality_df = quality_df.drop(cols_to_drop, axis=1)
    quality_df.columns = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
                          26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
                          48, 49, 50, 'cycle']
    max_cycle = quality_df['cycle'].max()
    grouped_by_cycle = quality_df.groupby(['cycle'])
    x = grouped_by_cycle.aggregate(np.sum)
    x = x.transpose()
    x = x.apply(calculate_percentage)
    print x
    os.system('mv /home/cuser/PycharmProjects/AMLpipeline/quality.csv /media/sf_sarah_share/')
    fig = plt.figure()
    ax = fig.add_subplot(111)
    my_cmap = plt.cm.get_cmap('GnBu')
    my_cmap.set_under(color='white')
    heatmap = ax.pcolor(x, cmap=my_cmap)

    major_xticks = np.arange(0, max_cycle + 20, 20)
    major_yticks = np.arange(0, 50, 10)
    ax.set_xticks(major_xticks)
    ax.set_yticks(major_yticks)
    ax.grid(b=True, which='both', color='0.85', linestyle='-')
    ax.spines['right'].set_color('0.85')
    ax.spines['top'].set_color('0.85')
    ax.spines['bottom'].set_color('0.85')
    ax.set_xlabel('Cycle')
    ax.set_ylabel('Q Score')
    plt.colorbar(heatmap)
    heatmap.set_clim(vmin=0, vmax=100)
    plt.show()

new_function()
