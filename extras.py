import argparse
from argparse import RawTextHelpFormatter
import os
import pandas as pd

# use of parser to create arguments for script to run from command line
parser = argparse.ArgumentParser(description="Multiplies given number by 2",
                                 epilog="This script relies on number being integer",
                                 formatter_class=RawTextHelpFormatter)

parser.add_argument('-n', '--number', action="store", dest='number', help=' A number to multiply', required=True,
                    type=int)
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


def parse_sample_sheet(csv_file):
    """Parses a sample sheet from a CSV file to produce...

    :param csv_file: sample sheet file
    :return: dataframe with relevant columns (?dictionary)
    """
    df_final = pd.DataFrame([])
    # 1) Read CSV file
    df_sample_sheet = pd.read_csv(csv_file)
    # 2) Find row index where the cell has "Sample_ID" and remove rows before this.
    for column in df_sample_sheet:
        for row_index, row in df_sample_sheet.iterrows():
            if row[column] == 'Sample_ID':
                id_index = row_index
                df_final = df_sample_sheet.ix[id_index+1:]
            else:
                pass
    # 3) If dataframe is empty "Sample_ID" was not found so raise AssertionError.
    assert not df_final.empty, "No Sample_ID column detected"
    return df_final

parse_sample_sheet('SampleSheet.csv')
