from illuminate import InteropTileMetrics, InteropControlMetrics, InteropErrorMetrics, InteropExtractionMetrics, \
    InteropIndexMetrics, InteropQualityMetrics, InteropCorrectedIntensityMetrics

tilemetrics = InteropTileMetrics('/media/sf_sarah_share/AML_data/InterOp/TileMetricsOut.bin')
indexmetrics = InteropIndexMetrics('/media/sf_sarah_share/AML_data/InterOp/IndexMetricsOut.bin')

print indexmetrics.df

