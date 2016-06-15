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
import fileinput

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
genome = '/media/genomicdata/ucsc_hg19/ucsc.hg19.fasta'
pindel2vcf = '/home/cuser/programs/pindel/pindel2vcf'
annovar = '/media/sf_sarah_share/ANNOVAR/annovar/table_annovar.pl'
varscan = '/home/cuser/programs/VarScan2/VarScan.v2.4.0.jar'
fastqc = '/home/cuser/programs/FastQC/fastqc'
delly = '/home/cuser/programs/delly_v0.7.3/delly'
bcftools = '/home/cuser/programs/samtools/bcftools-1.3.1/bcftools'

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

    os.system('mv /home/cuser/PycharmProjects/amlpipeline/output.pdf /media/sf_sarah_share/MiSeq_quality_outputs/')
    # move file to location.

    # Could then have something like:
    # if tilemetrics.mean_cluster_density in range(1100, 1300) or tilemetrics.percent_pf_clusters > 90:
    #       continue with pipeline
    # else:
    #       stop pipeline and return error message


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
@collate("*qfilter.fastq.gz", formatter("([^/]+)R[12]_001.qfilter.fastq.gz$"), "{path[0]}/{1[0]}.fastq.gz")
def run_fastqc_trimmed(infile, outfile):
    fastq1 = infile[0]
    fastq2 = infile[1]
    os.system("%s --extract %s" % (fastqc, fastq1))
    os.system("%s --extract %s" % (fastqc, fastq2))


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


'''
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
'''


@follows(run_samtools_stats)
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
    # bd_file = infile[:-14]
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
              "--name_of_logfile %spindel.log"  # space here for BD!!!!!
              # "-b %s.breakdancer_output.sv "  # file name with BreakDancer results     or -Q???
              # "-Q %s.breakdancer_with_pindel.sv"
              % (pindel, genome, infile, pindel_prefix, pindel_prefix))


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
@transform(["*.bwa.drm.sorted.bam"], suffix(".bwa.drm.sorted.bam"), ".indels.vs2.vcf")
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


@follows(run_varscan2_indels)
@transform(["*.bwa.drm.sorted.bam"], suffix(".bwa.drm.sorted.bam"), r"\1.bcf")
def call_delly(infile, outfile):
    sample = infile[:-19]
    blacklisted_regions = '/media/sf_sarah_share/excluded_regions_and_blacklisted.sorted.bed'
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
    vs2_data = {}
    for file in infile:
        if file.endswith(".del.delly.vcf"):
            delly_del_data = parse_vcfs.get_delly_output(file)
        if file.endswith(".dup.delly.vcf"):
            delly_dup_data = parse_vcfs.get_delly_output(file)
        if file.endswith(".inv.delly.vcf"):
            delly_inv_data = parse_vcfs.get_delly_output(file)
        if file.endswith(".ins.delly.vcf"):
            delly_ins_data = parse_vcfs.get_delly_output(file)
        if file.endswith(".tra.delly.vcf"):
            delly_tra_data = parse_vcfs.get_delly_output(file)
        if file.endswith(".pindel.merged.vcf"):
            pindel_data = parse_vcfs.get_pindel_output(file)
        if file.endswith(".indels.vs2.vcf"):
            vs2_data = parse_vcfs.get_vs2_output(file)
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
        for var in vs2_data.keys():
            output_file.write(parse_vcfs.print_vcf(vs2_data, var))


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
@transform(["*.annovar.vcf"], suffix(".annovar.vcf"), ".annovar.xlsx")
def vcf_to_excel(infile, outfile):
    vcf_reader = vcf.Reader(open(infile, 'r'))
    col_list = ['SAMPLE', 'Caller', 'CHROM', 'POS', 'REF', 'ALT', 'END', 'TYPE', 'SIZE', 'GT', 'GQ', 'DEPTH',
                'Alleles', 'AB', 'GENE', 'FUNC-refGene', 'EXONIC_FUNC-refGene']
    annovar_df = pd.DataFrame(columns=col_list)
    sample_id = infile[:-12]
    for record in vcf_reader:
        for sample in record:
            # PyVCF reader
            chr = record.CHROM
            pos = record.POS
            ref = record.REF
            alt = ",".join(str(a) for a in record.ALT)
            gt = sample['GT']
            if record.QUAL is None:
                gq = '.'
            else:
                gq = record.QUAL
            ad = sample['AD']
            info_dict = record.INFO
            end_pos = info_dict.get("END")
            sv_type = info_dict.get("SVTYPE")
            gene = ",".join(str(g) for g in info_dict.get("Gene.refGene"))
            func = ",".join(str(f) for f in info_dict.get("Func.refGene"))
            exonic_func = ",".join(str(f) for f in info_dict.get("ExonicFunc.refGene"))
            ref_reads = ad[0]
            alt_reads = ad[1]
            total_reads = ref_reads + alt_reads
            if alt_reads != 0 and ref_reads != 0:
                ab = (float(alt_reads) / float(ref_reads) * 100)
            else:
                ab = int(0)
            size = info_dict.get("SVLEN")
            ad_str = '%s,%s' % (str(ref_reads), str(alt_reads))
            caller = info_dict.get("Caller")
            if exonic_func is None:
                exonic_func_mod = '.'
            else:
                exonic_func_mod = exonic_func
            if caller == 'Delly':
                output_df = pd.DataFrame([[sample_id, caller, chr, pos, ref, alt, end_pos, sv_type, size, gt, gq,
                                           total_reads, ad_str, ab, gene, func, exonic_func_mod]], columns=col_list)
                annovar_df = annovar_df.append(output_df)
            else:
                if re.match("(.*)exonic(.*)", func) or re.match("(.*)splicing(.*)", func):
                    output_df = pd.DataFrame([[sample_id, caller, chr, pos, ref, alt, end_pos, sv_type, size, gt, gq,
                                               total_reads, ad_str, ab, gene, func, exonic_func_mod]], columns=col_list)
                    annovar_df = annovar_df.append(output_df)
                else:
                    pass

    writer = ExcelWriter('%s' % outfile)
    annovar_df.to_excel(writer, index=False)
    writer.save()
    os.system('cp /home/cuser/PycharmProjects/amlpipeline/%s '
              '/media/sf_sarah_share/160513_M04103_0019_000000000-ANU1A/outputs/merged_excels/pindel_only/' % outfile)


pipeline_run(verbose=4)

# Use these databases:
# refGene,knownGene,ensgene,esp6500_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_amr,1000g2014oct_eas,
# 1000g2014oct_eur,1000g2014oct_sas,snp137,clinvar_20140929,cosmic70,exac03
