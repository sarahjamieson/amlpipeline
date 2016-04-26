from ruffus import *
import os
import glob
from illuminate import InteropTileMetrics, InteropControlMetrics, InteropErrorMetrics, InteropExtractionMetrics, \
    InteropIndexMetrics, InteropQualityMetrics, InteropCorrectedIntensityMetrics
from create_pdf import CreatePDF

# paths to relevant tools
Trimmomatic = '/home/cuser/programs/Trimmomatic-0.36/trimmomatic-0.36.jar'
InterOp = '/media/sf_sarah_share/AML_data/InterOp/'


def assess_quality():
    """Obtains quality statistics from MiSeq InterOp files and produces PDF summary report.

    Notes:
        Uses Python module Illuminate 0.6.2. to read InterOp file data (https://pypi.python.org/pypi/illuminate/).
        See http://rpackages.ianhowson.com/bioc/savR/ for details on column meanings in InterOp files.

    """
    tilemetrics = InteropTileMetrics('%sTileMetricsOut.bin' % InterOp)
    controlmetrics = InteropControlMetrics('%sControlMetricsOut.bin' % InterOp)
    errormetrics = InteropErrorMetrics('%sErrorMetricsOut.bin' % InterOp)
    extractionmetrics = InteropExtractionMetrics('%sExtractionMetricsOut.bin' % InterOp)
    indexmetrics = InteropIndexMetrics('%sIndexMetricsOut.bin' % InterOp)
    qualitymetrics = InteropQualityMetrics('%sQMetricsOut.bin' % InterOp)
    corintmetrics = InteropCorrectedIntensityMetrics('%sCorrectedIntMetricsOut.bin' % InterOp)

    pdf = CreatePDF(tilemetrics, controlmetrics, errormetrics, extractionmetrics, indexmetrics, qualitymetrics,
                    corintmetrics)
    pdf.create_pdf()  # creates PDF document using LaTeX of quality data.

    os.system('mv /home/cuser/PycharmProjects/AMLpipeline/output.pdf /media/sf_sarah_share/')  # move file to location.

    # Could then have something like:
    # if tilemetrics.mean_cluster_density in range(1100, 1300) or tilemetrics.percent_pf_clusters > 90:
    #       continue with pipeline
    # else:
    #       stop pipeline and return error message


# Input: all files ending in fastq.gz; formatter collates files with same name ending either R1 or R2; outputs into file
# with common name. Uses Trimmomatic 0.36.
@collate("*fastq", formatter("([^/]+)R[12].fastq$"), "{1[0]}.fastq.gz")  # (input, filter, output)
@follows(assess_quality)
def quality_trim(infile, outfile):
    fastq1 = infile[0]  # first input file given
    fastq2 = infile[1]  # second input file given

    os.system("java -Xmx4g -jar %s "  # ?-Xmx4g
              "PE "  # paired-end mode
              "-threads 2 "  # multi-threading
              "-phred33 "  # use phred33 encoding for base quality (>= Illumina 1.8)
              "%s %s "  # input files R1 R2
              "%s.qfilter.fastq.gz %s.unpaired.fastq.gz "  # paired and unpaired output R1
              "%s.qfilter.fastq.gz %s.unpaired.fastq.gz "  # paired and unpaired output R2
              "CROP:150 "  # crop reads to max 150 from end
              "SLIDINGWINDOW:4:25 "  # perform trim on every 4 bases for min quality of 25
              "MINLEN:50"  # all reads should be min 50
              % (Trimmomatic, fastq1, fastq2, fastq1[:-6], fastq1, fastq2[:-6], fastq2))  # [:-6] removes file extension

    # glob module: finds path names matching specified pattern (https://docs.python.org/2/library/glob.html).
    filelist = glob.glob("*unpaired*")
    for f in filelist:
        os.remove(f)


pipeline_run(verbose=3)  # verbose = 3 gives medium level of information about run
