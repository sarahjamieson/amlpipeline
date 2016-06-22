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
from parse_sample_sheet import ParseSampleSheet
import shelve
from collections import defaultdict

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
'''
'''

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


pipeline_run(verbose=4)
'''
'''
listy = ['02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18',
         '19', '20', '21', '22', '23', '24']
for l in listy:
    os.system(
        "cp /home/cuser/PycharmProjects/amldata/%s_bd/%s*.breakdancer_output.sv "
        "/media/sf_sarah_share/160513_M04103_0019_000000000-ANU1A/outputs/breakdancer/"
        % (l, l))
'''
'''
annovar = '/media/sf_sarah_share/ANNOVAR/annovar/table_annovar.pl'


@transform(["*.breakdance_output.sv"], suffix(".breakdancer_output.sv"), ".breakdancer.vcf")
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
                "%s\t%s\t%s\t%s\t%s\t%s\t%s\tSVLEN=%s;NUM_READS=%s;CHR2=%s;END=%s;SVTYPE=%s;BD_SCORE=%s;ORIENTATION1=%s"
                "ORIENTATION2=%s\tGT:AD\t%s:%s\n" % (
                    chrom, pos, id, ref, alt, qual, filter, svlen, num_reads, chr2, end, type,
                    score, orientation1, orientation2, "0/0", "0,0"))


@follows(breakdancer_to_vcf)
@transform(["*.breakdancer.vcf"], suffix(".breakdancer.vcf"), ".bd.annovar.vcf")
def annotate_vcf(infile, outfile):
    os.system("%s "  # table_annovar.pl
              "%s "  # infile
              "/media/sf_sarah_share/ANNOVAR/annovar/humandb/ "  # directory database files are stored in
              "-buildver hg19 "  # genome build
              "-out %s "  # outfile
              "-remove "  # removes all temporary files
              "-protocol refGene,knownGene,ensgene,esp6500siv2_all,snp138, "  # databases
              "-arg '--exonicsplicing --hgvs --splicing 30',,,, "  # optional arguments, splicing threshold 30bp
              "-operation g,g,g,f,f "  # gene-based or filter-based for protocols
              "-nastring . "  # if variant doesn't have a score from tools (e.g. intronic variant & SIFT), position="."
              "-vcfinput"  # required if input is vcf file
              % (annovar, infile, outfile))

    # rename to match expected output
    os.rename("%s.hg19_multianno.vcf" % outfile, "%s" % outfile)
    # os.system("cp %s /media/sf_sarah_share/160513_M04103_0019_000000000-ANU1A/outputs/annovar/" % outfile)


@follows(annotate_vcf)
@transform(["*.bd.annovar.vcf"], suffix(".bd.annovar.vcf"), ".bd.annovar.xlsx")
def vcf_to_excel(infile, outfile):
    vcf_reader = vcf.Reader(open(infile, 'r'))
    col_list = ['SAMPLE', 'Caller', 'CHROM', 'POS', 'REF', 'ALT', 'CHR2', 'END', 'TYPE', 'SIZE', 'GT', 'DEPTH',
                'Alleles', 'AB', 'GENE', 'FUNC-refGene', 'EXONIC_FUNC-refGene']
    annovar_df = pd.DataFrame(columns=col_list)
    gasv_df = pd.DataFrame(columns=col_list)
    sample_id = infile[:-12]
    for record in vcf_reader:
        for sample in record:
            # PyVCF reader
            chr = record.CHROM
            pos = record.POS
            ref = record.REF
            alt = ",".join(str(a) for a in record.ALT)
            gt = sample['GT']
            ad = sample['AD']
            info_dict = record.INFO
            chr2 = info_dict.get("CHR2")
            end_pos = info_dict.get("END")
            sv_type = info_dict.get("SVTYPE")
            gene = ",".join(str(g) for g in info_dict.get("Gene.refGene"))
            if re.match("NONE(.*)", gene.upper()):
                gene_corr = "."
            else:
                gene_corr = gene
            func = ",".join(str(f) for f in info_dict.get("Func.refGene"))
            exonic_func = ",".join(str(f) for f in info_dict.get("ExonicFunc.refGene"))
            ref_reads = ad[0]
            alt_reads = ad[1]
            total_reads = ref_reads + alt_reads
            if alt_reads != 0 and ref_reads != 0:
                ab = (float(alt_reads) / float(total_reads) * 100)
            else:
                ab = int(0)
            size = info_dict.get("SVLEN")
            ad_str = '%s,%s' % (str(ref_reads), str(alt_reads))
            caller = info_dict.get("Caller")
            if re.match("None", exonic_func):
                exonic_func_mod = "."
            else:
                exonic_func_mod = exonic_func
            if caller == 'Delly':
                caller_prec = "%s_%s" % (caller, info_dict.get("Precision"))
                output_df = pd.DataFrame([[sample_id, caller_prec, chr, pos, ref, alt, chr2, end_pos, sv_type, size, gt,
                                           total_reads, ad_str, ab, gene_corr, func, exonic_func_mod]], columns=col_list)
                annovar_df = annovar_df.append(output_df)
            elif caller == 'GASV':
                output_df = pd.DataFrame([[sample_id, caller, chr, pos, ref, alt, chr2, end_pos, sv_type, size, ".",
                                           total_reads, ".", ".", gene_corr, func, exonic_func_mod]], columns=col_list)
                gasv_df = gasv_df.append(output_df)
            else:
                if re.match("(.*)exonic(.*)", func) or re.match("(.*)splicing(.*)", func):
                    output_df = pd.DataFrame([[sample_id, caller, chr, pos, ref, alt, chr2, end_pos, sv_type, size, gt,
                                               total_reads, ad_str, ab, gene_corr, func, exonic_func_mod]],
                                             columns=col_list)
                    annovar_df = annovar_df.append(output_df)
                else:
                    pass

    writer = ExcelWriter('%s' % outfile)
    annovar_df.to_excel(writer, sheet_name="Variants-Pindel-VS2-Delly", index=False)
    gasv_df.to_excel(writer, sheet_name="GASV", index=False)
    writer.save()

    os.system("cp %s /media/sf_sarah_share/160513_M04103_0019_000000000-ANU1A/outputs/merged_excels/w.gasv/" % outfile)


parser = argparse.ArgumentParser(description="Runs pipeline for Illumina MiSeq Nextera AML data.",
                                 formatter_class=RawTextHelpFormatter)

parser.add_argument('-s', '--sheet', action="store", dest='samplesheet', help='An illumina sample sheet for this run.',
                    required=True)
parser.add_argument('-f', '--fastq_dir', action='store', dest='fastqc_dir',
                    help='Directory containing Illumina MiSeq FASTQ files.')
parser.add_argument('-o', '--output_dir', action='store', dest='output_name',
                    help='Name to be given to directory for storing output files e.g. Run name.')
args = parser.parse_args()



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
'''
freq = "35|22|5|3.1"

hello = freq.split('|')

if int(hello[2]) > 10:
    "print to sheet2"
else:
    pass
