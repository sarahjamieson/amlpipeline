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
import argparse
from argparse import RawTextHelpFormatter

genome = '/media/genomicdata/ucsc_hg19/ucsc.hg19.fasta'
delly = '/home/cuser/programs/delly_v0.7.3/delly'
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
'''

'''
# @follows(index_bam)
@transform(["02*.bwa.drm.sorted.bam"], suffix(".bwa.drm.sorted.bam"), ".gasv.in.clusters")
def run_gasv(infile, outfile):
    name = infile[:-19]
    os.system("java -Xms512m -Xmx2048m -jar %s %s -OUTPUT_PREFIX %s -MAPPING_QUALITY 35"
              % (bamtogasv, infile, name))
    os.system("java -jar %s --minClusterSize 10 --batch %s.gasv.in" % (gasv, name))


@follows(run_gasv)
@transform(["*.gasv.in.clusters"], suffix(".gasv.in.clusters"), ".gasv.vcf")
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
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % sample]
    vcf_output = open(outfile, 'w')
    with open(infile, 'r') as gasv_file:
        reader = csv.reader(gasv_file, delimiter='\t')
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
    os.system("cp %s /media/sf_sarah_share/" % outfile)


pipeline_run(verbose=4, forcedtorun_tasks=[run_gasv, gasv_to_vcf])
'''
'''
listy = ['02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18',
         '19', '20', '21', '22', '23', '24']
for l in listy:
    os.system("mv /home/cuser/PycharmProjects/amlpipeline/%s*vs2.vcf /home/cuser/PycharmProjects/amldata/%s_bd/"
              % (l, l))
'''

'''
def breakdancer_to_vcf(bd_file):
    sample = bd_file[:-22]
    outfile = '%s.bd.vcf' % sample
    bed_out = '%s.bd.bed' % sample
    bed_output = open(bed_out, 'w')
    vcf_output = open(outfile, 'w')
    header = ['##fileformat=VCFv4.0\n',
              '##fileDate=%s\n' % curr_datetime,
              '##source=BreakDancer\n',
              '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n',
              '##INFO=<ID=NUM_READS,Number=1,Type=Integer,Description="Number of supporting read pairs">\n',
              '##INFO=<ID=BD_SCORE,Number=1,Type=Integer,Description="BreakDancer Score">\n',
              '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate">\n',
              '##INFO=<ID=END,Number=1,Type=Integer,Description="END coordinate">\n',
              '##INFO=<ID=ORIENTATION1,Number=1,Type=String,Description="BD orientation of chr1 in SV">\n',
              '##INFO=<ID=ORIENTATION2,Number=1,Type=String,Description="BD orientation of chr2 in SV">\n',
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % sample]
    with open(bd_file, 'r') as f:
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
            ref = "."
            alt = "."
            type = row[6]
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
                "%s\t%s\t%s\t%s\t%s\t%s\t%s\tSVLEN=%s;NUM_READS=%s;CHR2=%s;END=%s;BD_SCORE=%s;ORIENTATION1=%s"
                "ORIENTATION2=%s\t%s\t%s\n" % (chrom, pos, id, ref, alt, qual, filter, svlen, num_reads, chr2, end,
                                               score, orientation1, orientation2, " ", " "))
            bed_output.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, pos, end, ".", ".", "."))
    os.system("mv %s /media/sf_sarah_share/" % outfile)
    os.system("mv %s /media/sf_sarah_share/" % bed_out)


breakdancer_to_vcf(
    '/home/cuser/PycharmProjects/amldata/02_bd/02-D15-18331-AR-Nextera-Myeloid-Val1-Repeat_S2_L001_.'
    'breakdancer_output.sv')


parser = argparse.ArgumentParser(description="Runs pipeline for Illumina MiSeq Nextera AML data.",
                                 formatter_class=RawTextHelpFormatter)

parser.add_argument('-s', '--sheet', action="store", dest='samplesheet', help='An illumina sample sheet for this run.',
                    required=True)
parser.add_argument('-f', '--fastq_dir', action='store', dest='fastqc_dir',
                    help='Directory containing Illumina MiSeq FASTQ files.')
parser.add_argument('-o', '--output_dir', action='store', dest='output_name',
                    help='Name to be given to directory for storing output files e.g. Run name.')
args = parser.parse_args()
'''

'''
def get_delly_output(vcf_file):
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    delly_dict = {}
    loc1 = 0
    loc2 = 0
    for record in vcf_reader:
        for sample in record:
            info_dict = record.INFO
            if info_dict.get("IMPRECISE") is True:
                precision = 'IMPRECISE'
                ad = "%s,%s" % (sample['DR'], sample['DV'])
            else:
                precision = 'PRECISE'
                ad = "%s,%s" % (sample['RR'], sample['RV'])
            filter_str = ";".join(record.FILTER)
            if filter_str == '':
                filter = 'PASS'
            else:
                filter = 'LowQual'
            alt = ",".join(str(a) for a in record.ALT)
            svlen = info_dict.get("END") - record.POS
            variant_key = "%s_%s_%s_%s" % (record.CHROM, record.POS, record.REF, record.ALT[0])

            with open('cytoBand.txt', 'r') as f:
                reader = csv.reader(f, delimiter='\t')
                for row in reader:
                    if (record.CHROM == row[0]) and (int(row[1]) <= record.POS <= int(row[2])):
                        loc1 = row[3]
                        print loc1
                    if (info_dict.get("CHR2") == row[0]) and (int(row[1]) <= info_dict.get("END") <= int(row[2])):
                        loc2 = row[3]
                        print loc2
            if loc1 == loc2:
                loc = '%s%s' % (record.CHROM[3:], loc1)
            else:
                loc = '(%s%s;%s%s)' % (record.CHROM[3:], loc1, info_dict.get("CHR2")[3:], loc2)

            if variant_key in delly_dict:
                pass
            else:
                delly_dict[variant_key] = {"Caller": "Delly",
                                           "Chrom": record.CHROM,
                                           "Position": record.POS,
                                           "Ref": record.REF,
                                           "Alt": alt,
                                           "Qual": sample['GQ'],
                                           "Filter": filter,
                                           "Precision": precision,
                                           "SVTYPE": info_dict.get("SVTYPE"),
                                           "SVLEN": svlen,
                                           "CHR2": info_dict.get("CHR2"),
                                           "END": info_dict.get("END"),
                                           "MAPQ": info_dict.get("MAPQ"),
                                           "PE": info_dict.get("PE"),
                                           "Region": loc,
                                           "GT": sample['GT'],
                                           "AD": ad}
    print delly_dict

get_delly_output('/media/sf_sarah_share/05_tra.vcf')
'''


def get_run_info(csv_file):
    iem = ''
    investigator = ''
    experiment = ''
    run_date = ''
    workflow = ''
    app = ''
    assay = ''
    description = ''
    chemistry = ''
    worksheet = ''
    manifest = ''
    reads = ''
    data_index = 0
    read_index = 0
    manifest_index = 0
    read1 = 0
    read2 = 0
    sample_dict = {}

    with open(csv_file, 'r') as c:
        reader = csv.reader(c, delimiter=',')
        for i, row in enumerate(reader):
            if row[0] == 'IEMFileVersion':
                iem = row[1]
            elif row[0] == "Investigator Name":
                investigator = row[1]
            elif row[0] == 'Experiment Name':
                experiment = row[1]
            elif row[0] == 'Date':
                run_date = row[1]
            elif row[0] == 'Workflow':
                workflow = row[1]
            elif row[0] == 'Application':
                app = row[1]
            elif row[0] == 'Assay':
                assay = row[1]
            elif row[0] == 'Description':
                description = row[1]
            elif row[0] == 'Chemistry':
                chemistry = row[1]
            elif row[0] == 'worksheet':
                worksheet = row[1]
            elif row[0] == '[Manifests]':
                manifest_index = i
            elif row[0] == '[Reads]':
                read_index = i
            elif row[0] == '[Data]':
                data_index = i
            else:
                pass
            if i == (read_index + 1):
                read1 = row[0]
            if i == (read_index + 2):
                read2 = row[0]
            if i == (manifest_index + 1):
                manifest = row[1]
            reads = "(%s,%s)" % (read1, read2)

    run_dict = {
        "IEM": iem,
        "Investigator": investigator,
        "Experiment": experiment,
        "Date": run_date,
        "Workflow": workflow,
        "Application": app,
        "Assay": assay,
        "Description": description,
        "Chemistry": chemistry,
        "worksheet": worksheet,
        "Manifest": manifest,
        "Reads": reads
    }

    df_sample_sheet = pd.read_csv(csv_file, header=None)
    df_data = df_sample_sheet.ix[data_index + 1:]
    for row_index, row in df_data.iterrows():
        lab_id = str(row[1])[3:12]
        sample_id = row[0]
        name = row[1]
        plate = row[2]
        well = row[3]
        index1 = row[5]
        index2 = row[7]
        sample_manifest = row[8]
        project = row[10]

        sample_dict[lab_id] = {
            "Sample_id": sample_id,
            "Name": name,
            "Plate": plate,
            "Well": well,
            "Index1": index1,
            "Index2": index2,
            "Manifest": sample_manifest,
            "Project": project
        }

    print run_dict, sample_dict

get_run_info("ExampleSampleSheet.csv")




































