import vcf
import datetime
import csv
import re


def get_delly_output(vcf_file):
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    delly_dict = {}
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
                                           "GT": sample['GT'],
                                           "AD": ad}
    return delly_dict


def get_pindel_output(vcf_file):
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    pindel_dict = {}
    for record in vcf_reader:
        for sample in record:
            info_dict = record.INFO
            alt_temp = ",".join(str(a) for a in record.ALT)
            if alt_temp > 50:
                alt = alt_temp[:50]
            else:
                alt = alt_temp
            if record.QUAL is None:
                qual = "."
            else:
                qual = record.QUAL
            filter_str = ";".join(record.FILTER)
            if filter_str == '':
                filter = "."
            else:
                filter = filter_str
            variant_key = "%s_%s_%s_%s" % (record.CHROM, record.POS, record.REF, record.ALT[0])
            ad = ",".join(str(a) for a in sample['AD'])
            if variant_key in pindel_dict:
                pass
            else:
                pindel_dict[variant_key] = {"Caller": "Pindel",
                                            "Chrom": record.CHROM,
                                            "Position": record.POS,
                                            "Ref": record.REF,
                                            "Alt": alt,
                                            "Qual": qual,
                                            "Filter": filter,
                                            "Precision": ".",
                                            "SVTYPE": info_dict.get("SVTYPE"),
                                            "SVLEN": info_dict.get("SVLEN"),
                                            "CHR2": record.CHROM,
                                            "END": info_dict.get("END"),
                                            "MAPQ": 0,
                                            "PE": 0,
                                            "GT": sample['GT'],
                                            "AD": ad}

    return pindel_dict


def get_vs2_output(vcf_file):
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    vs2_dict = {}
    for record in vcf_reader:
        for sample in record:
            alt = ",".join(str(a) for a in record.ALT)
            ref = ",".join(str(a) for a in record.REF)
            filter_str = ";".join(record.FILTER)
            if filter_str == '':
                filter = "."
            else:
                filter = filter_str
            variant_key = "%s_%s_%s_%s" % (record.CHROM, record.POS, record.REF, record.ALT[0])
            ad = "%s,%s" % (sample['RD'], sample['AD'])
            svlen = len(alt) - len(ref)
            end = record.POS + svlen
            if variant_key in vs2_dict:
                pass
            else:
                vs2_dict[variant_key] = {"Caller": "VarScan2",
                                         "Chrom": record.CHROM,
                                         "Position": record.POS,
                                         "Ref": record.REF,
                                         "Alt": alt,
                                         "Qual": sample['GQ'],
                                         "Filter": filter,
                                         "Precision": ".",
                                         "SVTYPE": ".",
                                         "SVLEN": svlen,
                                         "CHR2": record.CHROM,
                                         "END": end,
                                         "MAPQ": 0,
                                         "PE": 0,
                                         "GT": sample['GT'],
                                         "AD": ad}

    return vs2_dict


def print_vcf(caller_dict, variant):
    variant_in_vcf_format = "%s\t%s\t%s\t%s\t%s\t%s\t%s\tPRECISION=%s;SVTYPE=%s;SVLEN=%s;CHR2=%s;END=%s;MAPQ=%s;" \
                            "Caller=%s;PE=%s\tGT:AD\t%s:%s\n" \
                            % (caller_dict[variant]['Chrom'],
                               caller_dict[variant]['Position'],
                               ".",
                               caller_dict[variant]['Ref'],
                               caller_dict[variant]['Alt'],
                               caller_dict[variant]['Qual'],
                               caller_dict[variant]['Filter'],
                               caller_dict[variant]['Precision'],
                               caller_dict[variant]['SVTYPE'],
                               caller_dict[variant]['SVLEN'],
                               caller_dict[variant]['CHR2'],
                               caller_dict[variant]['END'],
                               caller_dict[variant]['MAPQ'],
                               caller_dict[variant]['Caller'],
                               caller_dict[variant]['PE'],
                               caller_dict[variant]['GT'],
                               caller_dict[variant]['AD'])

    return variant_in_vcf_format


def gasv_to_vcf(gasv_file):
    sample = gasv_file[:-21]
    outfile = '%s.gasv.vcf' % sample
    curr_datetime = datetime.datetime.now().isoformat()
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

