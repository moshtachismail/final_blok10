"""
Author: Christel van Haren
The normalised counts file had versions of the gene names,
like ENSG00000260370.1, and removes the versions. So everything that
comes after the '.' will be removed and will look like ENSG00000260370.
Date: 04-08-2022
"""
import sys


def dots(input, output):
    """
    Removes the versions of the gene names.
    :param file: the normalised counts file with the versions of the
    gene names
    :param new_file: the new normalised counts file without the
    versions of the gene names
    :return: the new written file, with new_file in it.
    """
    print("version.py started")
    write_file = open(output, "w")
    with open(input, "r") as c:
        for line in c:
            gene_name = line.split("\t")
            # print(gene_name[0])
            no_dots = gene_name[0].split(".")
            # print(no_dots[0])
            new_column = line.replace(gene_name[0], no_dots[0])
            print(new_column)
            write_file.write(new_column)
    write_file.close()


def main(input, output):
    dots(input, output)
    print("version.py done")


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])