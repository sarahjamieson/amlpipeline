from illuminate import InteropTileMetrics, InteropControlMetrics, InteropErrorMetrics, InteropExtractionMetrics, \
    InteropIndexMetrics, InteropQualityMetrics, InteropCorrectedIntensityMetrics
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import plotly.plotly as py
import plotly.tools as tls
qualitymetrics = InteropQualityMetrics('/media/sf_sarah_share/AML_data/InterOp/QMetricsOut.bin')

quality_df = qualitymetrics.df

cols_to_drop = ['cycle', 'lane', 'tile']
quality_df = quality_df.drop(cols_to_drop, axis=1)
qual_scores = quality_df.sum(axis=0).values
qual_list = map(int, qual_scores.tolist())
less_than_q30 = qual_list[0: 28]
more_than_q30 = qual_list[29: 50]


myint = 1000000
y1 = [x / myint for x in less_than_q30]
y2 = [x / myint for x in more_than_q30]
x1 = range(1, 29)
x2 = range(29, 50)

plt.bar(x1, y1, color='b')
plt.bar(x2, y2, color='g')
plt.ylabel('Clusters (millions)')
plt.xlabel('Q score')
plt.show()

'''
x = [range(0, 50)]
y = something.tolist()
pyplot.bar(x, y, width=3.0)
pyplot.show()
'''



# remove lane and tile as well


'''
# Method for grouping by cycle
grouped = shortened_df.groupby(['cycle'])
for name, group in grouped:
    print name
    print group
print grouped.aggregate(np.sum)  # to sum all per cycle
'''


