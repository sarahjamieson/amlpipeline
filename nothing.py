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
'''


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


def get_sample_id_index(df_sample_sheet):
    for column in df_sample_sheet:
        for row_index, row in df_sample_sheet.iterrows():
            if row[column] == 'Sample_ID':
                id_index = row_index
                return id_index


def parse_sample_sheet(csv_file):
    df_sample_sheet = pd.read_csv(csv_file)
    id_index = get_sample_id_index(df_sample_sheet)
    df_final = df_sample_sheet.drop(df_sample_sheet.index[[1, (id_index-1)]])
    print df_final

parse_sample_sheet('SampleSheet.csv')
