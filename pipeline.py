from ruffus import *
import os
import glob
from illuminate import InteropTileMetrics, InteropControlMetrics, InteropErrorMetrics, InteropExtractionMetrics, \
    InteropIndexMetrics, InteropQualityMetrics, InteropCorrectedIntensityMetrics
from create_pdf import CreatePDF
import datetime
import warnings
warnings.simplefilter(action="ignore", category=FutureWarning)  # removes specific warning already seen

# paths to relevant tools
Trimmomatic = '/home/cuser/programs/Trimmomatic-0.36/trimmomatic-0.36.jar'
trim_adapters = '/home/cuser/programs/Trimmomatic-0.36/adapters/NexteraPE-PE.fa'
InterOp = '/media/sf_sarah_share/InterOp/'
genome = '/media/genomedata/ucsc_hg19/ucsc.hg19.fasta'
bwa = '/home/cuser/programs/bwa/bwa'
samblaster = '/home/cuser/programs/samblaster/samblaster'
samtools = '/home/cuser/programs/samtools/samtools/bin/samtools'
plot_bamstats = '/home/cuser/programs/samtools/samtools/bin/plot-bamstats'
bam2cfg = '/home/cuser/programs/breakdancer/build/bin/bam2cfg.pl'
breakdancer = '/home/cuser/programs/breakdancer/build/bin/breakdancer-max'
pindel = '/home/cuser/programs/pindel/pindel'
pindel2vcf = '/home/cuser/programs/pindel/pindel2vcf'
platypus = '/home/cuser/programs/Platypus/bin/Platypus.py'

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
    os.system("rm *bwa.drm.bam")


@follows(sort_bam)
@transform(["*.bwa.drm.sorted.bam"], suffix(".bwa.drm.sorted.bam"), ".bwa.drm.sorted.bam.bai")
def index_bam(infile, outfile):
    params = {'samtools': samtools,
              'input': infile}
    os.system("%(samtools)s index %(input)s" % params)


@follows(index_bam)
@transform(["*.bwa.drm.sorted.bam"], suffix(".bwa.drm.sorted.bam"), ".bwa.drm.sorted.bam.stats")
def run_samtools_stats(infile, outfile):
    os.system("%s stats %s > %s" % (samtools, infile, outfile))
    os.system("%s -p %s %s" % (plot_bamstats, outfile, outfile))


@transform(["*.bwa.drm.sorted.bam"], suffix(".bwa.drm.sorted.bam"), ".breakdancer_config")
def create_breakdancer_config(infile, outfile):
    os.system("perl %s %s > %s" % (bam2cfg, infile, outfile))  # infile = sorted, indexed bam, outfile = config_file


@follows(create_breakdancer_config)
def run_breakdancer(infile, outfile):
    os.system("%s -q 10 %s > %s" % (breakdancer, infile, outfile))  # infile = config file


def create_pindel_config(infile, outfile):
    config_file = open("%s" % outfile, "w")  # write into output file
    config_file.write("%s\t240\t%s\n" % (infile, infile[:-19]))  # name of BAM; insert size; sample_label
    config_file.close()


@follows(create_pindel_config)
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
              % (pindel, genome, infile, 'output_prefix', 'breakdancer_filename'))


@follows(run_pindel)
def pindel_to_vcf(infile, outfile):
    params = {
        'pindel2vcf': pindel2vcf,
        'input': infile,
        'reference': genome,
        'output': outfile
    }
    os.system("%(pindel2vcf)s "  # pindel2vcf program
              "-p %(input)s "  # input file
              "-r %(reference)s "  # reference genome on computer
              "-R hg19 "  # reference genome
              "-d hg19 "  # ??
              "-v %(output)s"  # output file
              % params)


# is nCPU=1 needed if default is 1??
# badReadsThreshold twice with different values
# https://wiki.gacrc.uga.edu/wiki/Platypus-Sapelo - arguments for platypus
def run_platypus(infile, outfile):
    params = {
        'platypus': platypus,
        'input': infile,
        'reference': genome,
        'output': outfile
    }
    os.system("python %(platypus)s callVariants "
              "--bamFiles=%(input)s "
              "--refFile=%(reference)s "
              "--output=%(output)s "
              "--regions=ROI_converted_to_platypus "  # ??
              "--minReads=10 "  # min reads required to support variant
              "--maxReads=2000000"  # why specify a maximum??
              "--minPosterior=1 "  # Phred score
              "--minFlank=0 "  # multiple variants within this distance of each other are discarded (0 = keep all)
              "--badReadsThreshold=0 "  # median base quality threshold in surrounding 7bp (0=keep all)
              "--minVarFreq=0 "  # threshold for fraction of reads supporting variant (0=keep all)
              "--minBaseQual=15 "  # min base quality score threshold
              "--filterDuplicates= "  # remove likely PCR copies
              "--minMapQual=0 "  # min mapping quality threshold (0=keep all)
              "--minGoodQualBases=1"  # min good quality bases in a read
              % params)


# verbose = 3 gives medium level of information about run
pipeline_run(verbose=6, target_tasks=sort_bam)
