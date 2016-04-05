from ruffus import *
import os
import glob
from illuminate import InteropTileMetrics, InteropControlMetrics, InteropErrorMetrics, InteropExtractionMetrics, \
    InteropIndexMetrics, InteropQualityMetrics, InteropCorrectedIntensityMetrics
from pylatex import Document, Section, Tabular, Package, Subsection
from latex import build_pdf

Trimmomatic = '/home/cuser/programs/Trimmomatic-0.36/trimmomatic-0.36.jar'


def assess_quality():
    tilemetrics = InteropTileMetrics('/media/sf_sarah_share/AML_data/InterOp/TileMetricsOut.bin')
    '''
    controlmetrics = InteropControlMetrics('/media/sf_sarah_share/AML_data/InterOp/ControlMetricsOut.bin')
    errormetrics = InteropErrorMetrics('/media/sf_sarah_share/AML_data/InterOp/ErrorMetricsOut.bin')
    extractionmetrics = InteropExtractionMetrics('/media/sf_sarah_share/AML_data/InterOp/ExtractionMetricsOut.bin')
    indexmetrics = InteropIndexMetrics('/media/sf_sarah_share/AML_data/InterOp/IndexMetricsOut.bin')
    qualitymetrics = InteropQualityMetrics('/media/sf_sarah_share/AML_data/InterOp/QMetricsOut.bin')
    corintmetrics = InteropCorrectedIntensityMetrics(
        '/media/sf_sarah_share/AML_data/InterOp/CorrectedIntMetricsOut.bin')
    '''

    doc = Document()
    doc.packages.append(Package('geometry', options=['tmargin=1cm', 'lmargin=10cm']))

    with doc.create(Section('Simple stuff')):
        doc.append('Identification of some kind')

        with doc.create(Subsection('Summary statistics')):
            with doc.create(Tabular('l|l')) as table:
                table.add_hline()
                table.add_row(('Mean Cluster Density', format(tilemetrics.mean_cluster_density/1000, '.2f')))
                table.add_row(('Mean Cluster Density PF', format(tilemetrics.mean_cluster_density_pf/1000, '.2f')))
                table.add_row(('Percentage of Clusters PF', '%s%%' % (format(tilemetrics.percent_pf_clusters, '.2f'))))

    doc.generate_pdf('output')

    # document = doc.generate_tex('output')
    # build_pdf(document)
    '''
    # making PDF file using ReportLab
    doc = SimpleDocTemplate("quality_output.pdf", pagesize=letter)
    elements = []
    data = [['Mean Cluster Density', format(tilemetrics.mean_cluster_density/1000, '.2f')],
            ['Mean PF Cluster Density', format(tilemetrics.mean_cluster_density_pf/1000, '.2f')],
            ['Percentage of Clusters PF', '%s%%' % (format(tilemetrics.percent_pf_clusters, '.2f'))]]
    summary_table = Table(data)
    summary_table.hAlign = "LEFT"  # aligns table to left of page
    elements.append(summary_table)
    doc.build(elements)
    '''

    # os.system("mv /home/cuser/PycharmProjects/AMLpipeline/quality_output.pdf /media/sf_sarah_share/")

    # Could then have something like:
    # if tilemetrics.mean_cluster_density in range(1100, 1300) or tilemetrics.percent_pf_clusters > 90:
    #       continue with pipeline
    # else:
    #       stop pipeline and return error message


assess_quality()

'''
# Input: all files ending in fastq.gz; formatter collates files with same name either R1 or R2; outputs into file with
# common name.
@collate("*fastq", formatter("([^/]+)R[12].fastq$"), "{1[0]}.fastq.gz")
@follows(assess_quality)
def quality_trim(infile, outfile):
    FASTQ1 = infile[0]  # first input file given
    FASTQ2 = infile[1]  # second input file given

    os.system("java -Xmx4g -jar %s PE -threads 2 -phred33 %s %s %s.qfilter.fastq.gz %s.unpaired.fastq.gz "
              "%s.qfilter.fastq.gz %s.unpaired.fastq.gz CROP:150 SLIDINGWINDOW:4:25 MINLEN:50"
              % (Trimmomatic, FASTQ1, FASTQ2, FASTQ1[:-6], FASTQ1, FASTQ2[:-6], FASTQ2))

    filelist = glob.glob("*unpaired*")
    for f in filelist:
        os.remove(f)


pipeline_run(verbose=3)
'''