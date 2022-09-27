"""
Author: Michelle Memelink
Date: 01-06-2022
Description: Using this script, it can be checked which SNPs contain the
patients that also appear in the SNP list. this list contains SNP sites found
in the genes of interest. The matching SNP sites are saved in a new file and
the genotype for the sites is noted for each patient.
"""

import gzip
import re
import sys


def main(filename_in_clin, filename_in_snp, filename_out):
    list_snps = read_clin_file(filename_in_clin)
    snp_list_dict = read_snp_file(filename_in_snp, filename_out)
    compare_data(list_snps, snp_list_dict, filename_out)


def read_clin_file(clin_file):
    """
    This function reads the vcf file, which contains the SNPs located on the
    genes encoding the proteins of interest. However, the patient test dataset
    only contains information from chromosome 10, so filtering has been done
    here and only the SNPs from chromosome 10 have been preserved. these SNPs
    are stored in a list (list_snps).
    :param clin_file: it is a vcf file which contains the SNPs of interest.
    These are the SNPs that like to be found in the patients.
    :return list_snps: Is a list containing SNP information of the SNPs located
    on chromosome 10 from the genes of interest. With the following format:
    [['1554841994', 'pathogenic', 'RPS24', 'A/C', 'chr10:78035589'], ....]
    """
    list_snps = []
    try:
        with open(clin_file, "r") as f:
            for line in f:
                split = line.replace("\n", "").split("\t")
                if split[-1].split(":")[0] == "10":
                    chr = "chr" + str(split[-1])
                    lijst = split[0:4]
                    lijst.append(chr)
                    list_snps.append(lijst)
    except FileNotFoundError:
        print("An error has occurred: wrong file of file path")
    except IndexError:
        print("The file is not complete or does not contain the "
              "correct information")
    except Exception as error:
        print("An error has occurred: " + str(error))
    print(list_snps)
    return list_snps


def read_snp_file(snp_file, snp_out):
    """
    This function reads the vcf file containing the SNPs from the patients.
    The information from this is then stored in a dictionary.
    :param snp_file: it is a vcf file containing the SNPs detected in the
    patients. Among other things, the genotype for all SNPs is available per
    patient, provided that genetic data was available for this.
    :param snp_out: This is the output file which is a txt file. This file
    ultimately contains the SNP locations that appear in the list of SNP
    locations of interest as well as among the patients. For all SNPs, the
    genotype for the relevant SNP is noted for each patient.
    :return snp_list_dict: This is a dictionary containing as key the chromosome
    and the exact location where the SNP is located. Lists are given as value.
    These lists contain per patient the genotype for the relevant location of
    the possible SNP. The dictonary has the following format:
    {'chr10:133783692': ['A', 'G', '0|1'], ..}
    """
    snp_list_dict = {}
    index_chrom, index_pos, index_ref, index_alt, index_form = "", "", "", "", ""
    count_snp_list = 0

    try:
        with gzip.open(snp_file, "rt") as file:
            for line in file:
                # print(line)
                if re.search("^#[a-zA-Z]", line):
                    line = line.split("\t")
                    for x in range(len(line)):
                        if line[x].__contains__("CHROM"):
                            index_chrom = x
                        if line[x].__contains__("POS"):
                            index_pos = x
                        if line[x].__contains__("REF"):
                            index_ref = x
                        if line[x].__contains__("ALT"):
                            index_alt = x
                        if line[x].__contains__("FORMAT"):
                            index_form = x

                    list_head = ["id"]
                    for i in range((int(index_form) + 1), len(line)):
                        list_head.append(line[i])
                    file_new = open(snp_out, "w")
                    file_new.write('\t'.join(list_head))

                elif line.startswith("chr"):
                    line = line.split("\t")
                    if index_chrom != "" and index_pos != "" and index_ref != "" \
                            and index_alt != "":
                        key = line[index_chrom] + ":" + line[index_pos]
                        tem_list = [line[index_ref], line[index_alt]]

                        first = int(index_form) + 1
                        last = len(line)
                        for i in range(first, last):
                            split = line[i].split(":")
                            tem_list.append(split[0])

                        snp_list_dict[key] = tem_list
                        count_snp_list += 1
    except FileNotFoundError:
        print("An error has occurred: wrong file of file path")
    except IndexError:
        print("The file is not complete or does not contain the "
              "correct information")
    except Exception as error:
        print("An error has occurred: " + str(error))

    return snp_list_dict


def compare_data(list_snps, snp_list_dict, snp_out):
    """
    This function checks which of the SNP locations found in the patients
    appear in the SNP list of interest. If a match is found here, the SNP is
    saved in a new file and the genotype is determined for the SNP location
    for each patient.
    :param list_snps: Is a list containing SNP information of the SNPs located
    on chromosome 10 from the genes of interest. With the following format:
    [['1554841994', 'pathogenic', 'RPS24', 'A/C', 'chr10:78035589'], ....]
    :param snp_list_dict: This is a dictionary containing as key the chromosome
    and the exact location where the SNP is located. Lists are given as value.
    These lists contain per patient the genotype for the relevant location of
    the possible SNP. The dictonary has the following format:
    {'chr10:133783692': ['A', 'G', '0|1'], ..}
    :param snp_out: This is the output file which is a txt file. This file
    ultimately contains the SNP locations that appear in the list of SNP
    locations of interest as well as among the patients. For all SNPs, the
    genotype for the relevant SNP is noted for each patient.
    :return: -
    """
    temp_list = []

    for key in list_snps:
        if key[-1] in snp_list_dict.keys():
            temp_list = [key[-1]]
            for gt in range(2, len(snp_list_dict[key[-1]])):
                value = snp_list_dict[key[-1]]
                if value[gt] in ["0/0", "0|0"]:
                    # Wildtype
                    temp_list.append("0")
                elif value[gt] in ["0/1", "0|1", "1/2", "1|2"]:
                    # Heterozygous variant
                    temp_list.append("1")
                elif value[gt] in ["1/1", "1|1", "2/2", "2|2"]:
                    # Homozygous variant
                    temp_list.append("2")
            if temp_list:
                file_new = open(snp_out, "a")
                file_new.write('\t'.join(temp_list) + '\n')


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
