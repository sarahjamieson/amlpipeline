import argparse
from argparse import RawTextHelpFormatter
import os
from ruffus.proxy_logger import *
from ruffus import *
from pybedtools import BedTool
import csv
import datetime
import re
curr_datetime = datetime.datetime.now().isoformat()

bamToBed = '/home/cuser/programs/bedtools2/bin/bamToBed'
pairDiscordants = '/home/cuser/programs/hydra/pairDiscordants.py'
dedupDiscordants = '/home/cuser/programs/hydra/dedupDiscordants.py'
hydra = '/home/cuser/programs/hydra/hydra'


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


def run_varscan2_indels(infile, outfile):
    os.system("%s mpileup "
              "-B "  # disables probabilistic realignment, reduces false SNPs caused by misalignments
              "-f %s "  # reference fasta
              "%s "
              "| "  # run output through VarScan
              "java -jar %s mpileup2indel "
              "--min-coverage 8 "
              "--min-reads2 2 "
              "--min-avg-qual 10 "
              "--p-value 99e-02 "
              "--output-vcf 1 "
              "--strand-filter 0 > %s"
              % (samtools, genome, infile, varscan, outfile))


# ----------------------------------------------------------------------------------------------------------
# BreakDancer to VCF format
def breakdancer_to_vcf(infile, outfile):
    sample = infile[:-22]
    vcf_output = open(outfile, 'w')
    header = ['##fileformat=VCFv4.0\n',
              '##fileDate=%s\n' % curr_datetime,
              '##source=BreakDancer\n',
              '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n',
              '##INFO=<ID=NUM_READS,Number=1,Type=Integer,Description="Number of supporting read pairs">\n',
              '##INFO=<ID=BD_SCORE,Number=1,Type=Integer,Description="BreakDancer Score">\n',
              '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate">\n',
              '##INFO=<ID=END,Number=1,Type=Integer,Description="END coordinate">\n',
              '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">\n',
              '##INFO=<ID=ORIENTATION1,Number=1,Type=String,Description="BD orientation of chr1 in SV">\n',
              '##INFO=<ID=ORIENTATION2,Number=1,Type=String,Description="BD orientation of chr2 in SV">\n',
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
              '##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allele depth">\n',
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % sample]
    with open(infile, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        i = 0
        while i < 5:
            reader.next()
            i += 1
        for line in header:
            vcf_output.write(line)
        for row in reader:
            chrom = row[0]
            pos = row[1]
            id = "."
            ref = "N"
            type = row[6]
            alt = "<%s>" % type
            qual = "."
            filter = "."
            end = row[4]
            svlen = row[7]
            num_reads = row[9]
            score = row[8]
            orientation1 = row[2]
            orientation2 = row[5]
            chr2 = row[3]
            vcf_output.write(
                "%s\t%s\t%s\t%s\t%s\t%s\t%s\tSVLEN=%s;NUM_READS=%s;CHR2=%s;END=%s;SVTYPE=%s;BD_SCORE=%s;ORIENTATION1=%s"
                "ORIENTATION2=%s\tGT:AD\t0/0:%s,0\n" % (
                    chrom, pos, id, ref, alt, qual, filter, svlen, num_reads, chr2, end, type,
                    score, orientation1, orientation2, num_reads))


# ----------------------------------------------------------------------------------------------------------------------
# Run GASV variant caller and convert to VCF format.

bamtogasv = '/home/cuser/programs/gasv/bin/BAMToGASV.jar'
gasv = '/home/cuser/programs/gasv/bin/GASV.jar'


def run_gasv(infile, outfile):
    name = infile[:-19]
    os.system("java -Xms512m -Xmx2048m -jar %s %s -OUTPUT_PREFIX %s -MAPPING_QUALITY 35"
              % (bamtogasv, infile, name))
    os.system("java -jar %s --minClusterSize 10 --batch %s.gasv.in" % (gasv, name))


def gasv_to_vcf(infile, outfile):
    sample = infile[:-19]
    header = ['##fileformat=VCFv4.0\n',
              '##fileDate=%s\n' % curr_datetime,
              '##source=GASVRelease_Oct1_2013\n',
              '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n',
              '##INFO=<ID=NUM_READS,Number=1,Type=Integer,Description="Number of supporting read pairs">\n',
              '##INFO=<ID=LOCAL,Number=1,Type=Float,Description="Breakpoint localization">\n',
              '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate">\n',
              '##INFO=<ID=END,Number=1,Type=Integer,Description="END coordinate">\n',
              '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of SV">\n',
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
              '##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allele depth">\n',
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % sample]
    vcf_output = open(outfile, 'w')
    with open(infile, 'r') as gasv_file:
        reader = csv.reader(gasv_file, delimiter='\t')
        reader.next()
        for line in header:
            vcf_output.write(line)
        for attribute in reader:
            chrom = "chr%s" % attribute[1]
            pos = attribute[2].split(',')[0]
            id = '.'
            ref = 'N'
            type = attribute[7]
            if type == 'D':
                alt = '<DEL>'
            elif re.match("I(.*)", type):
                alt = '<INV>'
            elif type == 'V':
                alt = '<DIV>'
            elif re.match("T(.*)", type):
                alt = '<TRA>'
            else:
                alt = '.'
            qual = "."
            filter = "."
            end = attribute[4].split(',')[0]
            svlen = int(end) - int(pos)
            num_reads = attribute[5]
            local = attribute[6]
            chr2 = "chr%s" % attribute[3]
            vcf_output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\tSVLEN=%s;NUM_READS=%s;LOCAL=%s;SVTYPE=%s;CHR2=%s;END=%s\t"
                             "GT:AD\t0/0:%s,0\n" % (chrom, pos, id, ref, alt, qual, filter, svlen, num_reads, local,
                                                    type, chr2, end, num_reads))
# ----------------------------------------------------------------------------------------------------------------------
