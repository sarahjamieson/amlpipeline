from illuminate import InteropTileMetrics, InteropControlMetrics, InteropErrorMetrics, InteropExtractionMetrics, \
    InteropIndexMetrics, InteropQualityMetrics, InteropCorrectedIntensityMetrics
import numpy as np
import pandas as pd
import os
import datetime
from pandas import ExcelWriter
from ruffus import *

Trimmomatic = '/home/cuser/programs/Trimmomatic-0.36/trimmomatic-0.36.jar'
trim_adapters = '/home/cuser/programs/Trimmomatic-0.36/adapters/NexteraPE-PE.fa'
InterOp = '/media/sf_sarah_share/InterOp/'
bwa = '/home/cuser/programs/bwa/bwa'
samblaster = '/home/cuser/programs/samblaster/samblaster'
samtools = '/home/cuser/programs/samtools/samtools/bin/samtools'
plot_bamstats = '/home/cuser/programs/samtools/samtools/bin/plot-bamstats'
bam2cfg = '/home/cuser/programs/breakdancer/breakdancer-max1.4.5/bam2cfg.pl'
breakdancer = '/home/cuser/programs/breakdancer/breakdancer-max'
pindel = '/home/cuser/programs/pindel/pindel'
genome = '/media/sf_sarah_share/ucsc_hg19/ucsc.hg19.fasta'
platypus = '/home/cuser/programs/Platypus/bin/Platypus.py'
pindel2vcf = '/home/cuser/programs/pindel/pindel2vcf'
annovar = '/home/cuser/programs/ANNOVAR/annovar/table_annovar.pl'

curr_datetime = datetime.datetime.now().isoformat()


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

    os.system('mv /home/cuser/PycharmProjects/amlpipeline/output.pdf /media/sf_sarah_share/')  # move file to location.

    # Could then have something like:
    # if tilemetrics.mean_cluster_density in range(1100, 1300) or tilemetrics.percent_pf_clusters > 90:
    #       continue with pipeline
    # else:
    #       stop pipeline and return error message


assess_quality()


# Input: all files ending in fastq.gz; formatter collates files with same name ending either R1 or R2; outputs into file
# with common name. Uses Trimmomatic 0.36.
@collate("*.fastq", formatter("([^/]+)R[12].fastq$"), "{1[0]}.fastq.gz")  # (input, filter, output)
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
              "ILLUMINACLIP:%s:2:30:10 "  # adaptor trimming, seed matches (16bp), allow 2 mismatches,
              # PE reads Q30, SE reads Q10
              "SLIDINGWINDOW:4:25 "  # perform trim on every 4 bases for min quality of 25
              "MINLEN:50"  # all reads should be min 50
              % (Trimmomatic, fastq1, fastq2, fastq1[:-6], fastq1, fastq2[:-6], fastq2,
                 trim_adapters))  # [:-6] removes file extension

    # glob module: finds path names matching specified pattern (https://docs.python.org/2/library/glob.html).
    filelist = glob.glob("*unpaired*")
    for f in filelist:
        os.remove(f)


@follows(quality_trim)
@collate("*qfilter.fastq.gz", formatter("([^/]+)R[12].qfilter.fastq.gz$"), "{1[0]}.bwa.drm.bam")
def align_bwa(infile, outfile):
    fastq1 = infile[0]
    fastq2 = infile[1]
    sample = fastq1[:-20]
    id = 'bwa'
    rg_header = '@RG\tID:%s\tCN:WMRGL\tDS:TruSight_Myeloid_Nextera_v1.1\tDT:%s' % (sample, curr_datetime)

    os.system("%s mem -M -a "  # mark shorter split hits as secondary, output alignments for SE and unpaired PE
              "-t 2 "  # number of threads
              "-k 18 "  # min seed length
              "-R \'%s\' "  # read group header
              "%s "  # genome (bwa indexed hg19)
              "%s %s"  # input R1, input R2
              "| sed \'s/@PG\tID:%s/@PG\tID:%s/\' - "  # "-" requests standard output, useful when combining tools
              "| %s -M -r -u "
              "> %s.samblaster.log - "
              "| %s view -Sb - > %s"  # -Sb puts output through samtools sort
              % (bwa, rg_header, genome, fastq1, fastq2, id, sample, samblaster, sample, samtools, outfile))

@follows(align_bwa)
@transform(["*.bwa.drm.bam"], suffix(".bwa.drm.bam"), ".bwa.drm.sorted.bam")
def sort_bam(infile, outfile):
    os.system("%s sort "
              "-@ 2 "  # number of threads = 2
              "-m 2G "  # memory per thread = 2G
              "-O bam "  # output file format = .bam
              "-T %s "  # temporary file name
              "-o %s "  # output file
              "%s"  # input file
              % (samtools, infile, outfile, infile))


@follows(sort_bam)
@transform(sort_bam, suffix(".bwa.drm.sorted.bam"), ".bwa.drm.sorted.bam.bai")
def index_bam(infile, outfile):
    os.system("%s index %s" % (samtools, infile))


@follows(index_bam)
@transform(sort_bam, suffix(".bwa.drm.sorted.bam"), ".bwa.drm.sorted.bam.stats")
def run_samtools_stats(infile, outfile):
    os.system("%s stats %s > %s" % (samtools, infile, outfile))
    # os.system("%s -p %s %s" % (plot_bamstats, outfile, outfile))  # error: no such file or directory gnuplot


@follows(run_samtools_stats)
@transform(sort_bam, suffix(".bwa.drm.sorted.bam"), ".breakdancer_config.txt")
def create_breakdancer_config(infile, outfile):
    config_file = open("%s" % outfile, "w")  # write into output file
    config_file.write("map:%s\tmean:240\tstd:40\treadlen:36.00\tsample:%s\texe:bwa-0.7.12\n" % (infile, infile[:-19]))
    config_file.close()


@follows(create_breakdancer_config)
@transform(create_breakdancer_config, suffix(".breakdancer_config.txt"), ".breakdancer_output.sv")
def run_breakdancer(infile, outfile):
    os.system("%s -q 10 %s > %s" % (breakdancer, infile, outfile))


@follows(run_breakdancer)
@transform(sort_bam, suffix(".bwa.drm.sorted.bam"), ".pindel_config.txt")
def create_pindel_config(infile, outfile):
    config_file = open("%s" % outfile, "w")  # write into output file
    config_file.write("%s\t240\t%s\n" % (infile, infile[:-19]))  # name of BAM; insert size; sample_label
    config_file.close()


@follows(create_pindel_config)
@transform(create_pindel_config, formatter(), ["{path[0]}/{basename[0]}._BP",
                                               "{path[0]}/{basename[0]}._D",
                                               "{path[0]}/{basename[0]}._INV",
                                               "{path[0]}/{basename[0]}._LI",
                                               "{path[0]}/{basename[0]}._RP",
                                               "{path[0]}/{basename[0]}._SI",
                                               "{path[0]}/{basename[0]}._TD"])
def run_pindel(infile, outfile):
    os.system("%s "  # pindel program
              "-f %s "  # reference genome
              "-i %s "
              "-o %s "  # ??output_prefix
              "-c ALL "  # look at all chromosomes
              "-T 2 "  # number of threads
              "-w 50 "  # window size
              "-x 4 "  # maximum size of SVs to detect, 4 = 8092    why 4???
              "-v 6 "  # minimum inversion size in bases     why 6???
              "-e 0.01 "  # expected sequencing error rate
              "-b %s"  # file name with BreakDancer results     or -Q???
              % (pindel, genome, infile, infile[:-18], 'cataracts.breakdancer_output.sv'))


@follows(run_pindel)
@transform(run_pindel, formatter(), "{path[0]}/{basename[0]}.pindel.merged.vcf", "{path[0]}/{basename[0]}.")
def pindel_to_vcf(infile, outfile, pindel_prefix):
    os.system("%s "  # pindel2vcf program
              "-P %s "  # input file
              "-r %s "  # reference genome on computer
              "-R hg19 "  # reference genome
              "-d hg19 "  # ??
              "-v %s"  # output file
              % (pindel2vcf, pindel_prefix, genome, outfile))


@follows(run_pindel)
@transform(sort_bam, suffix(".bwa.drm.sorted.bam"), ".platypus_output.unsorted.vcf")
def run_platypus(infile, outfile):
    os.system("python %s callVariants "
              "--bamFiles=%s "
              "--refFile=%s "
              "--output=%s "
              "--minReads=10 "  # min reads required to support variant
              "--maxReads=2000000 "  # why specify a maximum??
              "--minPosterior=1 "  # Phred score
              "--minFlank=0 "  # multiple variants within this distance of each other are discarded (0 = keep all)
              "--badReadsThreshold=0 "  # median base quality threshold in surrounding 7bp (0=keep all)
              "--minVarFreq=0 "  # threshold for fraction of reads supporting variant (0=keep all)
              "--minBaseQual=15 "  # min base quality score threshold
              "--filterDuplicates=1 "  # remove likely PCR copies
              "--minMapQual=0 "  # min mapping quality threshold (0=keep all)
              "--minGoodQualBases=1 "  # min good quality bases in a read
              "--nCPU=?"  # number of processors/cores, parallel jobs
              % (platypus, infile, genome, outfile))


# the number of arguments for "-protocol", "-arg" and "-operation" should all be equal and match in order.
@follows(run_platypus)
@transform(run_platypus, suffix(".platypus_output.unsorted.vcf"), ".annovar.vcf")
def annotate_vcf(infile, outfile):
    os.system("%s "  # table_annovar.pl
              "%s "  # infile
              "/home/cuser/programs/ANNOVAR/annovar/humandb/ "  # directory database files are stored in
              "-buildver hg19 "  # genome build
              "-out %s "  # outfile
              "-remove "  # removes all temporary files
              "-protocol refGene,knownGene,ensgene,esp6500siv2_all, "  # databases
              "-arg '-exonicsplicing -hgvs -splicing 30',, "  # optional arguments, splicing threshold 30bp
              "-operation g,f,f "  # gene-based or filter-based for protocols
              "-nastring . "  # if variant doesn't have a score from tools (e.g. intronic variant & SIFT), position="."
              "-vcfinput"  # required if input is vcf file
              % (annovar, infile, outfile))

    # rename vcf file to match expected
    os.rename("%s.hg19.multianno.vcf" % outfile, '%s' % outfile)


pipeline_run(verbose=4, forcedtorun_tasks=pindel_to_vcf)

# use these databases for actually running:
# refGene,knownGene,ensgene,esp6500_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_amr,1000g2014oct_eas,
# 1000g2014oct_eur,1000g2014oct_sas,snp137,clinvar_20140929,cosmic70,exac03
