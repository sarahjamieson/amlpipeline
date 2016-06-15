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
from pylatex import Document, Section, Tabular, Package, Command, Figure, SubFigure, Description, Subsection
from pylatex.utils import NoEscape, bold
import re
import csv

genome = '/media/genomicdata/ucsc_hg19/ucsc.hg19.fasta'
delly = '/home/cuser/programs/delly_v0.7.3/delly'
bam = '/home/cuser/PycharmProjects/amldata/02/02-D15-18331-AR-Nextera-Myeloid-Val1-Repeat_S2_L001_.bwa.drm.sorted.bam'
bcftools = '/home/cuser/programs/samtools/bcftools-1.3.1/bcftools'
bamtogasv = '/home/cuser/programs/gasv/bin/BAMToGASV.jar'
gasv = '/home/cuser/programs/gasv/bin/GASV.jar'
bwa = '/home/cuser/programs/bwa/bwa'
samtools = '/home/cuser/programs/samtools/samtools/bin/samtools'
samblaster = '/home/cuser/programs/samblaster/samblaster'

curr_datetime = datetime.datetime.now().isoformat()

'''
@collate("*qfilter.fastq.gz", formatter("([^/]+)R[12]_001.qfilter.fastq.gz$"), "{path[0]}/{1[0]}.bwa.drm.bam")
def align_bwa(infile, outfile):
    fastq1 = infile[0]
    fastq2 = infile[1]
    sample = fastq1[:-20]
    rg_header = '@RG\tID:%s\tSM:%s\tCN:WMRGL\tDS:TruSight_Myeloid_Nextera_v1.1\tDT:%s' % (sample, sample, curr_datetime)

    os.system("%s mem -t 2 "  # mark shorter split hits as secondary, output alignments for SE and unpaired PE
              "-k 18 "  # number of threads
              "-aM "  # min seed length
              "-R \'%s\' "  # read group header
              "%s "  # genome (bwa indexed hg19)
              "%s %s"  # input R1, input R2
              "| sed \'s/@PG\tID:bwa/@PG\tID:%s/\' - "  # "-" requests standard output, useful when combining tools
              "| %s --removeDups 2> "
              "%s.samblaster.log --splitterFile %s.splitreads.sam --discordantFile %s.discordant.sam "
              "| %s view -Sb - > %s"  # -Sb puts output through samtools sort
              % (bwa, rg_header, genome, fastq1, fastq2, sample, samblaster, sample, sample, sample, samtools, outfile))


@follows(align_bwa)
@transform(["02*.bwa.drm.bam"], suffix(".bwa.drm.bam"), ".bwa.drm.sorted.bam")
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
@transform("02*.bwa.drm.sorted.bam", suffix(".bwa.drm.sorted.bam"), ".bwa.drm.sorted.bam.bai")
def index_bam(infile, outfile):
    os.system("%s index %s" % (samtools, infile))


@follows(index_bam)
@transform(["02*.bwa.drm.sorted.bam"], suffix(".bwa.drm.sorted.bam"), ".bam.gasv.in.clusters")
def run_gasv(infile, outfile):
    name = infile[:-19]
    os.system("java -Xms512m -Xmx2048m -jar %s %s" % (bamtogasv, infile))
    os.system("java -jar %s --batch %s.bam.gasv.in" % (gasv, name))
'''
'''
def gasv_to_vcf(gasv_file):
    sample = gasv_file[:-21]
    outfile = '%s.gasv.vcf' % sample
    header = ['##fileformat=VCFv4.0\n',
              '##fileDate=%s\n' % curr_datetime,
              '##source=GASVRelease_Oct1_2013\n',
              '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n',
              '##INFO=<ID=NUM_READS,Number=1,Type=Integer,Description="Number of supporting read pairs">\n',
              '##INFO=<ID=LOCAL,Number=1,Type=Float,Description="Breakpoint localization">\n',
              '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate">\n',
              '##INFO=<ID=END,Number=1,Type=Integer,Description="END coordinate">\n',
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % sample]
    vcf_output = open(outfile, 'w')
    with open(gasv_file, 'r') as gasv:
        reader = csv.reader(gasv, delimiter='\t')
        reader.next()
        for line in header:
            vcf_output.write(line)
        for attribute in reader:
            chrom = attribute[1]
            pos = attribute[2].split(',')[0]
            id = '.'
            ref = '.'
            type = attribute[7]
            if type == 'D':
                alt = '<DEL>'
            elif re.match("I(.*)", type):
                alt = '<INV>'
            elif type == 'V':
                alt = '<DIV>'
            elif type == 'T':
                alt = '<TRA>'
            else:
                alt = '<UNKNOWN>'
            qual = "."
            filter = "."
            end = attribute[4].split(',')[0]
            svlen = int(end) - int(pos)
            num_reads = attribute[5]
            local = attribute[6]
            chr2 = attribute[3]
            vcf_output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\tSVLEN=%s;NUM_READS=%s;LOCAL=%s;CHR2=%s;END=%s\t%s\t%s\n"
                             % (chrom, pos, id, ref, alt, qual, filter, svlen, num_reads, local, chr2, end, " ", " "))

gasv_to_vcf('/home/cuser/programs/02-D15.bwa.drm.sorted.bam.gasv.in.clusters')
'''
listy = ['02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18',
         '19', '20', '21', '22', '23', '24']
for l in listy:
    os.system("mv /home/cuser/PycharmProjects/amlpipeline/%s*vs2.vcf /home/cuser/PycharmProjects/amldata/%s_bd/"
              % (l, l))
