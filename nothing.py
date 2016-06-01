from illuminate import InteropTileMetrics, InteropControlMetrics, InteropErrorMetrics, InteropExtractionMetrics, \
    InteropIndexMetrics, InteropQualityMetrics, InteropCorrectedIntensityMetrics
import numpy as np
import pandas as pd
import os
import datetime
from pandas import ExcelWriter
from ruffus import *
from ruffus.proxy_logger import *
from pybedtools import BedTool
import vcf
from zipfile import ZipFile
from pylatex import Document, Section, Tabular, Package, Command, Figure, SubFigure, Description
from pylatex.utils import NoEscape, bold

'''
samtools = '/home/cuser/programs/samtools/samtools/bin/samtools'
genome = '/media/sf_sarah_share/ucsc_hg19/ucsc.hg19.fasta'
varscan = '/home/cuser/programs/VarScan2/VarScan.v2.4.0.jar'
bamToBed = '/home/cuser/programs/bedtools2/bin/bamToBed'
pairDiscordants = '/home/cuser/programs/hydra/pairDiscordants.py'
dedupDiscordants = '/home/cuser/programs/hydra/dedupDiscordants.py'
hydra = '/home/cuser/programs/hydra/hydra'

(logger_proxy, logging_mutex) = make_shared_logger_and_proxy(
    setup_std_shared_logger,
    "aml_logger",
    {"file_name": "/home/cuser/PycharmProjects/amlpipeline/aml_pipeline.log"})


@transform(["*.csv"], suffix(".csv"), ".bed", logger_proxy, logging_mutex)
def first_task(infile, outfile, logger_proxy, logging_mutex):
    with logging_mutex:
        logger_proxy.info("First task completed")
    bedfile = BedTool("%s" % infile)
    bedfile.saveas("%s" % outfile)


# ------------------------------------------------------------------------------
# for Hydra-SV
@transform(["*.bwa.drm.sorted.bam"], suffix(".bwa.drm.sorted.bam"), ".bedpe")
def convert_bam_to_bed(infile, outfile):
    os.system("%s -i %s -tag NM | "
              "%s -i stdin -m hydra -z 800 > %s"
              % (bamToBed, infile, pairDiscordants, outfile))


@follows(convert_bam_to_bed)
@transform(convert_bam_to_bed, suffix(".bedpe"), ".deduped.bedpe")
def dedup_discordants(infile, outfile):
    os.system("%s -i %s -s 3 > %s" % (dedupDiscordants, infile, outfile))


@follows(dedup_discordants)
@transform(dedup_discordants, suffix(".deduped.bedpe"), ".breaks")
def run_hydra(infile, outfile):
    os.system("%s -in %s -out %s" % (hydra, infile, outfile))
# ------------------------------------------------------------------------------



def get_output():
    vcf_file = 'cataracts.annovar.vcf'
    sample = vcf_file[:-12]
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    for record in vcf_reader:
        for sample in record:
            # PyVCF reader
            chr = record.CHROM
            pos = record.POS
            ref = record.REF
            alt = ",".join(str(a) for a in record.ALT)
            type = record.var_type
            gt = sample['GT']
            if record.QUAL is None:
                gq = '.'
            else:
                gq = record.QUAL
            nv = sample['NV']

            # ANNOVAR info
            info_dict = record.INFO
            gene = ",".join(str(g) for g in info_dict.get("Gene.refGene"))
            depth = info_dict.get("TC")
            ab = "%.1f" % ((float(nv)/float(depth))*100)
            tcf = info_dict.get("TCF")
            tcr = info_dict.get("TCR")
            nf = ",".join(str(f) for f in info_dict.get("NF"))
            nr = ",".join(str(r) for r in info_dict.get("NR"))
            qd = info_dict.get("QD")
            mq = ",".join(str(m) for m in info_dict.get("MQ"))
            sbpval = ",".join(str(s) for s in info_dict.get("SbPval"))
            hp = info_dict.get("HP")
            func = ",".join(str(f) for f in info_dict.get("ExonicFunc.refGene"))
            # HGVS in AAChange.refGene for exonic and in GeneDetail.refGene for intronic
            hgvs = ",".join(str(h) for h in info_dict.get("AAChange.refGene"))
'''
'''

fastqc = '/home/cuser/programs/FastQC/fastqc'
Trimmomatic = '/home/cuser/programs/Trimmomatic-0.36/trimmomatic-0.36.jar'
trim_adapters = '/home/cuser/programs/Trimmomatic-0.36/adapters/NexteraPE-PE.fa'


@collate("*.fastq.gz", formatter("([^/]+)R[12].fastq.gz$"), "{1[0]}.qfilter.fastq.gz")  # (input, filter, output)
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
              "ILLUMINACLIP:%s:2:30:10 "  # adaptor trimming, seed matches (16bp), allow 2 mismatches,
              # PE reads Q30, SE reads Q10
              "CROP:150 "  # crop reads to max 150 from end
              "SLIDINGWINDOW:4:25 "  # perform trim on every 4 bases for min quality of 25
              "MINLEN:50"  # all reads should be min 50
              % (Trimmomatic, fastq1, fastq2, fastq1[:-9], fastq1, fastq2[:-9], fastq2,
                 trim_adapters))  # [:-6] removes file extension


@follows(quality_trim)
@collate("*qfilter.fastq.gz", formatter("([^/]+)R[12].qfilter.fastq.gz$"), "{1[0]}.fastq.gz")
def run_fastqc(infile, outfile):
    fastq1 = infile[0]
    fastq2 = infile[1]
    os.system("%s --extract %s" % (fastqc, fastq1))
    os.system("%s --extract %s" % (fastqc, fastq2))


pipeline_run()
'''
# with ZipFile('04-R1.qfilter_fastqc.zip', 'w') as myzip:
# myzip.extractall()

pdflatex = '/usr/local/texlive/2015/bin/x86_64-linux/pdflatex'


summary_df = pd.read_table('04-R1.qfilter_fastqc/summary.txt', header=None, names=['Score', 'Parameter'],
                           usecols=[0, 1])
score_list = summary_df['Score'].tolist()
parameter_list = summary_df['Parameter'].tolist()
summary_dict = dict(zip(parameter_list, score_list))

basic_stats_df = pd.read_table('04-R1.qfilter_fastqc/fastqc_data.txt', header=None, names=['Property', 'Value'],
                               usecols=[0, 1], skiprows=3, nrows=7)

doc = Document()
doc.packages.append(Package('geometry', options=['tmargin=0.75in', 'lmargin=0.75in', 'rmargin=0.75in']))
doc.append(Command('begin', 'center'))
doc.append(Command('Large', bold('FastQC Quality Results')))
doc.append(Command('end', 'center'))

with doc.create(Section('Basic Statistics')):
    with doc.create(Description()) as desc:
        for row_index, row in basic_stats_df.iterrows():
            desc.add_item("%s:" % row['Property'], "%s" % row['Value'])

with doc.create(Section('Per base sequence quality: %s' % summary_dict.get('Per base sequence quality'))):
    with doc.create(Figure(position='htbp')) as plot:
        plot.add_image('04-R1.qfilter_fastqc/Images/per_base_quality.png', width='60.0\textwidth')

with doc.create(Section('Per tile sequence quality: %s' % summary_dict.get('Per tile sequence quality'))):
    with doc.create(Figure(position='htbp', placement=NoEscape(r'\centering'))) as plot:
        plot.add_image('04-R1.qfilter_fastqc/Images/per_tile_quality.png')


doc.generate_pdf('fastqc', clean_tex=False, compiler=pdflatex)
os.system('mv /home/cuser/PycharmProjects/amlpipeline/fastqc.pdf /media/sf_sarah_share/MiSeq_quality_outputs/')
