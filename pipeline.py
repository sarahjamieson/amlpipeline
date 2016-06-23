from illuminate import InteropTileMetrics, InteropControlMetrics, InteropErrorMetrics, InteropExtractionMetrics, \
    InteropIndexMetrics, InteropQualityMetrics, InteropCorrectedIntensityMetrics
import os
import datetime
from ruffus import *
import glob
from create_pdf import CreatePDF
import vcf
import pandas as pd
import re
from pandas import ExcelWriter
import parse_vcfs
import argparse
from argparse import RawTextHelpFormatter
import csv
import shelve
from collections import defaultdict

curr_datetime = datetime.datetime.now().isoformat()
'''
parser = argparse.ArgumentParser(description="Runs pipeline for Illumina MiSeq Nextera AML data.",
                                 formatter_class=RawTextHelpFormatter)
parser.add_argument('-s', '--sheet', action="store", dest='sample_sheet', help='An illumina sample sheet for this run.',
                    required=True, default='/media/sf_sarah_share/160513_M04103_0019_000000000-ANU1A/SampleSheet.csv')
parser.add_argument('-d', '--result_dir', action='store', dest='result_dir', required=True,
                    help='Directory containing Illumina MiSeq InterOp and Data folders.',
                    default='/media/sf_sarah_share/160513_M04103_0019_000000000-ANU1A/')
parser.add_argument('-o', '--output_dir', action='store', dest='output_dir', required=True,
                    help='Directory for output to be stored in',
                    default='/media/sf_sarah_share/MiSeq_Nextera_Results/')
args = parser.parse_args()
'''

Trimmomatic = '/home/cuser/programs/Trimmomatic-0.36/trimmomatic-0.36.jar'
trim_adapters = '/home/cuser/programs/Trimmomatic-0.36/adapters/NexteraPE-PE.fa'
# InterOp = '/%s/InterOp/' % args.result_dir
InterOp = '/media/sf_sarah_share/160620_M01622_0286_000000000-AR4UH/InterOp/'
bwa = '/home/cuser/programs/bwa/bwa'
samblaster = '/home/cuser/programs/samblaster/samblaster'
samtools = '/home/cuser/programs/samtools/samtools/bin/samtools'
plot_bamstats = '/home/cuser/programs/samtools/samtools/bin/plot-bamstats'
bam2cfg = '/home/cuser/programs/breakdancer/breakdancer-max1.4.5/bam2cfg.pl'
breakdancer = '/home/cuser/programs/breakdancer/breakdancer-max'
pindel = '/home/cuser/programs/pindel/pindel'
genome = '/media/genomicdata/ucsc_hg19/ucsc.hg19.fasta'
pindel2vcf = '/home/cuser/programs/pindel/pindel2vcf'
annovar = '/media/sf_sarah_share/ANNOVAR/annovar/table_annovar.pl'
varscan = '/home/cuser/programs/VarScan2/VarScan.v2.4.0.jar'
fastqc = '/home/cuser/programs/FastQC/fastqc'
delly = '/home/cuser/programs/delly_v0.7.3/delly'
bcftools = '/home/cuser/programs/samtools/bcftools-1.3.1/bcftools'
bamtogasv = '/home/cuser/programs/gasv/bin/BAMToGASV.jar'
gasv = '/home/cuser/programs/gasv/bin/GASV.jar'
script_dir = os.path.dirname(os.path.abspath(__file__))
# os.system("cp %s/Data/Intensities/BaseCalls/*.fastq.gz %s/" % (args.result_dir, script_dir))
# os.system("cp /media/sf_sarah_share/160620_M01622_0286_000000000-AR4UH/Data/*.fastq.gz %s/" % script_dir)


'''
os.system("mkdir %s%s" % (output_dir, worksheet))
parse_sheet = ParseSampleSheet(args.sample_sheet)
run_dict, sample_dict = parse_sheet.parse_sample_sheet()
worksheet = run_dict.get('worksheet')
'''
worksheet = 'run2'


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
                    corintmetrics, worksheet)
    pdf.create_pdf()  # creates PDF document using LaTeX of quality data.

    os.system('mv %s_InterOp_Results.pdf /media/sf_sarah_share/160620_M01622_0286_000000000-AR4UH/QC/' % "194766")
    # os.system('cp %s_InterOp_Results.pdf %s%s/' % (worksheet args.output_dir, worksheet))

    # Could then have something like:
    # if tilemetrics.mean_cluster_density in range(1100, 1300) or tilemetrics.percent_pf_clusters > 90:
    #       continue with pipeline
    # else:
    #       stop pipeline and return error message

assess_quality()


@collate("*.fastq.gz", formatter("([^/]+)R[12]_001.fastq.gz$"), "{path[0]}/{1[0]}.fastq.gz")
def run_fastqc(infile, outfile):
    fastq1 = infile[0]
    fastq2 = infile[1]
    os.system("%s --extract %s" % (fastqc, fastq1))
    os.system("%s --extract %s" % (fastqc, fastq2))


# Input: all files ending in fastq.gz; formatter collates files with same name ending either R1 or R2; outputs into file
# with common name. Uses Trimmomatic 0.36.
@follows(run_fastqc)
@collate("*.fastq.gz", formatter("([^/]+)R[12]_001.fastq.gz$"),
         "{path[0]}/{1[0]}.qfilter.fastq.gz")  # (input, filter, output)
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

    # glob module: finds path names matching specified pattern (https://docs.python.org/2/library/glob.html).
    filelist = glob.glob("*unpaired*")
    for f in filelist:
        os.remove(f)


@follows(quality_trim)
@collate("*.qfilter.fastq.gz", formatter("([^/]+)R[12]_001.qfilter.fastq.gz$"), "{path[0]}/{1[0]}.fastq.gz")
def run_fastqc_trimmed(infile, outfile):
    fastq1 = infile[0]
    fastq2 = infile[1]
    os.system("%s --extract %s" % (fastqc, fastq1))
    os.system("%s --extract %s" % (fastqc, fastq2))
    joint_sample_name = fastq1[:-23]
    from fastqc import CreateFastQCPDF
    pdf = CreateFastQCPDF(joint_sample_name)
    pdf.create_pdf()


@follows(run_fastqc_trimmed)
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
@transform("*.bwa.drm.sorted.bam", suffix(".bwa.drm.sorted.bam"), ".bwa.drm.sorted.bam.bai")
def index_bam(infile, outfile):
    os.system("%s index %s" % (samtools, infile))


@follows(index_bam)
@transform("*.bwa.drm.sorted.bam", suffix(".bwa.drm.sorted.bam"), ".bwa.drm.sorted.bam.stats")
def run_samtools_stats(infile, outfile):
    os.system("%s stats %s > %s" % (samtools, infile, outfile))
    os.system("%s -p %s %s" % (plot_bamstats, outfile, outfile))
    joint_sample_name = infile[:-19]
    from fastqc import CreateFastQCPDF
    pdf = CreateFastQCPDF(joint_sample_name)
    pdf.create_pdf()


@follows(run_samtools_stats)
@transform("*.bwa.drm.sorted.bam", suffix(".bwa.drm.sorted.bam"), ".breakdancer_config.txt")
def create_breakdancer_config(infile, outfile):
    config_file = open("%s" % outfile, "w")  # write into output file
    config_file.write("map:%s\tmean:160\tstd:50\treadlen:150\tsample:%s\texe:bwa-0.7.12\n" % (infile, infile[:-19]))
    config_file.close()


@follows(create_breakdancer_config)
@transform(create_breakdancer_config, suffix(".breakdancer_config.txt"), ".breakdancer_output.sv")
def run_breakdancer(infile, outfile):
    os.system("%s -q 10 %s > %s" % (breakdancer, infile, outfile))


@follows(run_breakdancer)
@transform(["*.bwa.drm.sorted.bam"], suffix(".bwa.drm.sorted.bam"), ".pindel_config")
def create_pindel_config(infile, outfile):
    config_file = open("%s" % outfile, "w+")  # write into output file
    config_file.write("%s\t240\t%s\n" % (infile, infile[:-19]))  # name of BAM; insert size; sample_label
    config_file.close()


@follows(create_pindel_config)
@transform("*.pindel_config", formatter(), ["{path[0]}/{basename[0]}._BP",
                                            "{path[0]}/{basename[0]}._D",
                                            "{path[0]}/{basename[0]}._INV",
                                            "{path[0]}/{basename[0]}._LI",
                                            "{path[0]}/{basename[0]}._RP",
                                            "{path[0]}/{basename[0]}._SI",
                                            "{path[0]}/{basename[0]}._TD"], "{path[0]}/{basename[0]}.")
def run_pindel(infile, outfile, pindel_prefix):
    start = datetime.datetime.now().isoformat()
    print start
    bd_file = infile[:-14]
    os.system("%s "  # pindel program
              "-f %s "  # reference genome
              "-i %s "
              "-c ALL "  # ??output_prefix
              "-o %s "  # look at all chromosomes
              "--minimum_support_for_event 5 "
              "-T 2 "  # number of threads
              "-w 50 "  # window size
              "-x 4 "  # maximum size of SVs to detect, 4 = 8092    why 4???
              "--sensitivity 0.96 "
              "-v 6 "  # minimum inversion size in bases     why 6???
              "-e 0.01 "  # expected sequencing error rate
              "--report_breakpoints TRUE "
              "--report_long_insertions TRUE "
              "--name_of_logfile %spindel.log "
              "-b %s.breakdancer_output.sv "  # file name with BreakDancer results     or -Q???
              "-Q %s.breakdancer_with_pindel.sv"
              % (pindel, genome, infile, pindel_prefix, pindel_prefix, bd_file, bd_file))
    end = datetime.datetime.now().isoformat()
    print end


# https://github.com/ELIXIR-ITA-training/VarCall2015/blob/master/Tutorials/T5.2_variantcalling_stucturalvariants_tutorial.md
# for pindel2vcf parameters
@follows(run_pindel)
@transform(run_pindel, formatter(), "{path[0]}/{basename[0]}.pindel.merged.vcf", "{path[0]}/{basename[0]}.")
def pindel_to_vcf(infile, outfile, pindel_prefix):
    os.system("%s "  # pindel2vcf program
              "-P %s "  # input file
              "-r %s "  # reference genome on computer
              "-R hg19 "  # reference genome
              "-d 2009 "
              "--min_size 5 "  # minimum size of events (what measurement??)
              "--min_coverage 10 "  # minimum number of reads
              "--min_supporting_reads 5 "
              "--max_size 8092 "
              "-v %s "  # output file
              "--het_cutoff 0.1 "
              "--hom_cutoff 0.85 "
              "-G"
              % (pindel2vcf, pindel_prefix, genome, outfile))


@follows(pindel_to_vcf)
@transform(["*.bwa.drm.sorted.bam"], suffix(".bwa.drm.sorted.bam"), r"\1.bcf")
def call_delly(infile, outfile):
    sample = infile[:-19]
    blacklisted_regions = 'excluded_regions.excl.bed'
    os.system("%s call -t DEL -x %s -o %s.del.bcf -g %s %s" % (delly, blacklisted_regions, sample, genome, infile))
    os.system("%s call -t TRA -x %s -o %s.tra.bcf -g %s %s" % (delly, blacklisted_regions, sample, genome, infile))
    os.system("%s call -t INV -x %s -o %s.inv.bcf -g %s %s" % (delly, blacklisted_regions, sample, genome, infile))
    os.system("%s call -t INS -x %s -o %s.ins.bcf -g %s %s" % (delly, blacklisted_regions, sample, genome, infile))
    os.system("%s call -t DUP -x %s -o %s.dup.bcf -g %s %s" % (delly, blacklisted_regions, sample, genome, infile))


@follows(call_delly)
@transform(["*.bcf"], suffix(".bcf"), r"\1.delly.vcf")
def delly_to_vcf(infile, outfile):
    os.system("%s view %s > %s" % (bcftools, infile, outfile))


@follows(delly_to_vcf)
@collate(["*.vcf"], regex(r"(.+)\.(.+)\.(.+)\.vcf"), r"\1.unified.vcf")
def combine_vcfs(infile, outfile):
    delly_del_data = {}
    delly_dup_data = {}
    delly_inv_data = {}
    delly_ins_data = {}
    delly_tra_data = {}
    pindel_data = {}
    for f in infile:
        if f.endswith(".del.delly.vcf"):
            delly_del_data = parse_vcfs.get_delly_output(f)
        if f.endswith(".dup.delly.vcf"):
            delly_dup_data = parse_vcfs.get_delly_output(f)
        if f.endswith(".inv.delly.vcf"):
            delly_inv_data = parse_vcfs.get_delly_output(f)
        if f.endswith(".ins.delly.vcf"):
            delly_ins_data = parse_vcfs.get_delly_output(f)
        if f.endswith(".tra.delly.vcf"):
            delly_tra_data = parse_vcfs.get_delly_output(f)
        if f.endswith(".pindel.merged.vcf"):
            pindel_data = parse_vcfs.get_pindel_output(f)
    sample = outfile[:-12]
    header = ['##fileformat=VCFv4.0\n',
              '##fileDate=%s\n' % curr_datetime,
              '##source=Delly,Pindel\n',
              '##FILTER=<ID=PASS,Description="All filters passed">\n',
              '##INFO=<ID=Precision,Number=1,Type=String,Description="Precise structural variation">\n',
              '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n',
              '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">\n',
              '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation">\n',
              '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">\n',
              '##INFO=<ID=Caller,Number=1,Type=String,Description="Caller used">\n',
              '##INFO=<ID=MAPQ,Number=1,Type=Integer,Description="Median mapping quality of paired-ends">\n',
              '##INFO=<ID=PE,Number=1,Type=Integer,Description="Paired-end support of the structural variant">\n',
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
              '##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allele depth">\n',
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % sample]
    with open("%s" % outfile, 'w+') as output_file:
        for line in header:
            output_file.write(line)
        for var in delly_del_data.keys():
            if delly_del_data[var]['Filter'] == 'PASS':
                output_file.write(parse_vcfs.print_vcf(delly_del_data, var))
        for var in delly_dup_data.keys():
            if delly_dup_data[var]['Filter'] == 'PASS':
                output_file.write(parse_vcfs.print_vcf(delly_dup_data, var))
        for var in delly_ins_data.keys():
            if delly_ins_data[var]['Filter'] == 'PASS':
                output_file.write(parse_vcfs.print_vcf(delly_ins_data, var))
        for var in delly_inv_data.keys():
            if delly_inv_data[var]['Filter'] == 'PASS':
                output_file.write(parse_vcfs.print_vcf(delly_inv_data, var))
        for var in delly_tra_data.keys():
            if delly_tra_data[var]['Filter'] == 'PASS':
                output_file.write(parse_vcfs.print_vcf(delly_tra_data, var))
        for var in pindel_data.keys():
            output_file.write(parse_vcfs.print_vcf(pindel_data, var))


# the number of arguments for "-protocol", "-arg" and "-operation" should all be equal and match in order.
@follows(combine_vcfs)
@transform(["*.unified.vcf"], suffix(".unified.vcf"), ".annovar.vcf")
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
@transform(["*.annovar.vcf"], suffix(".annovar.vcf"), ".variant_id.db")
def compress_variants(infile, outfile):
    # shelve object containing variant dictionary for each sample, contained within a database file.
    shelf = shelve.open(outfile, writeback=False)
    variant_dict_for_sample = defaultdict(list)
    sample = infile[:-12]

    with open(infile, 'r') as inputfile:
        for line in inputfile:
            line_split = line.strip().split('\t')
            if line_split[0].startswith("#"):
                pass
            else:
                variant_id = "%s_%s_%s_%s" % (line_split[0], line_split[1], line_split[3], line_split[4])
                variant_dict_for_sample[variant_id] = [sample]
        shelf.update(variant_dict_for_sample)
        shelf.close()


@follows(compress_variants)
@merge(compress_variants, "master_variants.txt")
def combine_variants(infile, outfile):
    # produces dictionary where each variant is the key and the value is a list containing the worksheet and sample id.
    variant_dict_for_run = defaultdict(list)
    for inputfile in infile:
        sample_variants = shelve.open(inputfile)
        for k, v in sample_variants.items():
            variant_dict_for_run[k].append((worksheet, ','.join(v)))
        sample_variants.close()  # need to close to save

    master_database = shelve.open("AML_variant_frequency.db")
    # only populate if empty to have something to check against, will add on later if not.
    if len(master_database) == 0:
        master_database.update(variant_dict_for_run)
    else:
        pass

    # if already in master database, add sample to existing sample(s) so it doesn't overwrite; otherwise just add.
    combined_dict = defaultdict(list)
    for k, v in variant_dict_for_run.items():
        if k in master_database:
            combined_dict[k] = master_database[k] + v
        else:
            combined_dict[k] = v

    master_database.update(combined_dict)
    master_database.close()

    # Create new dictionary with unique sample numbers to remove duplicate sample entries which may be caused by
    # running the script multiple times.
    analysis_dict = defaultdict(list)
    for k, v in combined_dict.items():
        value_set = set()
        for item in v:
            value_set.add(item)
        value_list = list(value_set)
        analysis_dict[k] = value_list

    shelf = shelve.open("AML_variant_frequency.db", writeback=False)
    shelf.update(analysis_dict)
    shelf.close()


@follows(combine_variants)
@transform(["*.annovar.vcf"], suffix(".annovar.vcf"), ".annovar.final.vcf")
def add_freq_to_vcf(infile, outfile):
    master_database = shelve.open("AML_variant_frequency.db")
    total_samples_list = []
    for k, v in master_database.items():
        total_samples_list.extend([item[1] for item in v])
        total_samples = len(set(total_samples_list))

    with open("%s" % infile, 'r') as input_vcf:
        with open("%s" % outfile, 'w+') as output_vcf:
            for row in input_vcf:
                if row.startswith("#"):
                    output_vcf.write(row)
                else:
                    row_split = row.strip().split('\t')
                    variant = "%s_%s_%s_%s" % (row_split[0], row_split[1], row_split[3], row_split[4])
                    chrom = row_split[0]
                    pos = row_split[1]
                    id = row_split[2]
                    ref = row_split[3]
                    alt = row_split[4]
                    qual = row_split[5]
                    filter = row_split[6]
                    info = row_split[7]
                    format = row_split[8]
                    format_val = row_split[9]
                    if master_database[variant]:
                        total_runs = len(set(item[0] for item in master_database[variant]))
                        total_vars = len(set(item[1] for item in master_database[variant]))
                        pct = round((float(total_vars)/float(total_samples))*100, 1)
                        freq = "%s|%s|%s|%s" % (total_vars, total_samples, pct, total_runs)
                    else:
                        freq = "."
                    output_vcf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s;FREQ=%s\t%s\t%s\n"
                                     % (chrom, pos, id, ref, alt, qual, filter, info, freq, format, format_val))


@follows(add_freq_to_vcf)
@transform(["*.annovar.final.vcf"], suffix(".annovar.final.vcf"), ".annovar.xlsx")
def vcf_to_excel(infile, outfile):
    vcf_reader = vcf.Reader(open(infile, 'r'))
    col_list = ['SAMPLE', 'Caller', 'CHROM', 'POS', 'REF', 'ALT', 'CHR2', 'END', 'TYPE', 'SIZE', 'GT', 'DEPTH',
                'Alleles', 'AB', 'GENE', 'FUNC-refGene', 'EXONIC_FUNC-refGene', 'FREQ', 'AF']
    annovar_df = pd.DataFrame(columns=col_list)
    gasv_df = pd.DataFrame(columns=col_list)
    bd_df = pd.DataFrame(columns=col_list)
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
            if re.match("None", gene.upper()):
                gene_mod = "."
            else:
                gene_mod = gene
            func = ",".join(str(f) for f in info_dict.get("Func.refGene"))
            exonic_func = ",".join(str(f) for f in info_dict.get("ExonicFunc.refGene"))
            ref_reads = ad[0]
            alt_reads = ad[1]
            total_reads = ref_reads + alt_reads
            freq = ",".join(str(a) for a in info_dict.get("FREQ"))
            freq_split = freq.split('|')
            allele_freq = float(freq_split[2])
            if alt_reads != 0 and ref_reads != 0:
                ab = (float(alt_reads) / float(total_reads) * 100)
            else:
                ab = int(0)
            size = info_dict.get("SVLEN")
            ad_str = '%s,%s' % (str(ref_reads), str(alt_reads))
            caller = info_dict.get("Caller")
            prec = info_dict.get("Precision")
            caller_prec = "%s_%s" % (caller, prec)
            if re.match("None", exonic_func):
                exonic_func_mod = '.'
            else:
                exonic_func_mod = exonic_func
            if re.match("None", func):
                func_mod = '.'
            else:
                func_mod = func
            if caller == 'Delly':
                output_df = pd.DataFrame([[sample_id, caller_prec, chr, pos, ref, alt, chr2, end_pos, sv_type, size, gt,
                                           total_reads, ad_str, ab, gene_mod, func_mod, exonic_func_mod, freq, allele_freq]],
                                         columns=col_list)
                annovar_df = annovar_df.append(output_df)
            else:
                if re.match("(.*)exonic(.*)", func) or re.match("(.*)splicing(.*)", func):
                    output_df = pd.DataFrame([[sample_id, caller, chr, pos, ref, alt, chr2, end_pos, sv_type, size, gt,
                                               total_reads, ad_str, ab, gene_mod, func_mod, exonic_func_mod, freq, allele_freq]],
                                             columns=col_list)
                    annovar_df = annovar_df.append(output_df)
                else:
                    pass

    writer = ExcelWriter('%s' % outfile)
    annovar_df.to_excel(writer, sheet_name="Pindel-VS2-Delly", index=False)
    gasv_df.to_excel(writer, sheet_name="GASV", index=False)
    bd_df.to_excel(writer, sheet_name="BreakDancer", index=False)
    writer.save()

    '''
    os.system("mkdir %s%s/%s/" % (args.output_dir, worksheet, sample_id))
    os.system("mkdir %s%s/%s/Data/" % (args.output_dir, worksheet, sample_id))
    os.system("mkdir %s%s/%s/Results/" % (args.output_dir, worksheet, sample_id))
    os.system("mv %s.annovar.xlsx %s%s/%s/Results/" % (sample_id, args.output_dir, worksheet, sample_id))
    os.system("mv %s_Sample_Quality.pdf %s%s/%s/Results/" % (sample_id, args.output_dir, worksheet, sample_id))
    os.system("mv %s.png %s/quality_figs/" % script_dir)  # to put into app
    os.system("mv %s* %s%s/%s/Data/" % (sample_id, args.output_dir, worksheet, sample_id))
    '''


pipeline_run(verbose=4)
