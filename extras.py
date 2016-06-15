import argparse
from argparse import RawTextHelpFormatter
import os
from ruffus.proxy_logger import *
from ruffus import *
from pybedtools import BedTool

bamToBed = '/home/cuser/programs/bedtools2/bin/bamToBed'
pairDiscordants = '/home/cuser/programs/hydra/pairDiscordants.py'
dedupDiscordants = '/home/cuser/programs/hydra/dedupDiscordants.py'
hydra = '/home/cuser/programs/hydra/hydra'

# use of parser to create arguments for script to run from command line
parser = argparse.ArgumentParser(description="Runs pipeline for AML data.",
                                 formatter_class=RawTextHelpFormatter)

parser.add_argument('-s', '--sheet', action="store", dest='samplesheet', help='An illumina sample sheet for this run.',
                    required=True)
parser.add_argument('-d', '--directory', action='store', dest='directory', help='')
args = parser.parse_args()


def perform_calculation(number):
    print number * 2


perform_calculation(args.number)


# work in progress for TIGRA-ext - can't "make" tool
def run_tigra():
    os.system("perl TIGRA-ext.pl "
              "-o %s "  # output as ./sv.vcf
              "-r %s "  # reference genome
              "-k \"15,25\" "  # kmer size. Add 35 if reads >100bp.
              "-a 50 "  # assembly 50bp into SV breakpoints (recommended value)
              "-l 500 "  # assembly 500bp from SV breakpoints (recommended value)
              "-b "  # specifies breakdancer format
              "-T %s "  # TIGRA binary directory
              "%s "  # breakdancer output file
              "%s")  # bwa-aligned bam files


# work in progress for ANNOVAR
def annotate_vcf(infile, outfile):
    os.system("%s "  # table_annovar.pl
              "%s "  # infile
              "humandb/ "  # directory database files are stored in
              "-buildver hg19 "  # genome build
              "-out %s "  # outfile
              "-remove "  # removes all temporary files
              "-protocol refGene,knownGene,ensgene,esp6500_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_amr,"
              "1000g2014oct_eas,1000g2014oct_eur,1000g2014oct_sas,snp137,clinvar_20140929,cosmic70,exac03 "  # databases
              "-arg '-exonicsplicing -hgvs -splicing 30',,,,,,,,,,, "  # optional arguments, splicing threshold 30bp
              "-operation g,g,g,f,f,f,f,f,f,f,f,f,f,f "  # gene-based or filter-based for protocols
              "-nastring . "  # if variant doesn't have a score from tools (e.g. intronic variant & SIFT), position="."
              "-vcfinput"  # required if input is vcf file
              % ('annovar', infile, outfile))

# ---------------------------------------------------
# Ruffus logging
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
# -------------------------------------------------------------------

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

# ------------------------------------------------------------------------------
# SNV detection: Platypus & VarScan2
samtools = '/home/cuser/programs/samtools/samtools/bin/samtools'
platypus = '/home/cuser/programs/Platypus/bin/Platypus.py'
varscan = '/home/cuser/programs/VarScan2/VarScan.v2.4.0.jar'
genome = '/media/genomicdata/ucsc_hg19/ucsc.hg19.fasta'


@transform(["*.bwa.drm.sorted.bam"], suffix(".bwa.drm.sorted.bam"), ".platypus_output.unsorted.vcf")
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
              "--nCPU=2"  # number of processors/cores, parallel jobs
              % (platypus, infile, genome, outfile))


# need to run samtools mpileup to create pileup file. VarScan requires pileup file as input (run mpileup2snp for
# multiple samples).
# @follows(run_platypus)
@transform(["*.bwa.drm.sorted.bam"], suffix(".bwa.drm.sorted.bam"), ".snps.vs2.vcf")
def run_varscan2_snps(infile, outfile):
    os.system("%s mpileup "
              "-B "  # disables probabilistic realignment, reduces false SNPs caused by misalignments
              "-f %s "  # reference fasta
              "%s "
              "| "  # run output through VarScan
              "java -jar %s mpileup2snp "
              "--min-coverage 8 "
              "--min-reads2 2 "
              "--min-avg-qual 10 "
              "--p-value 99e-02 "
              "--output-vcf 1 "
              "--strand-filter 0 > %s"
              % (samtools, genome, infile, varscan, outfile))
# ----------------------------------------------------------------------------------------------------------
