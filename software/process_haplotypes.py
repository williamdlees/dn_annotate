# Read the changeo file and determine which J -genes can act as anchors, i.e. have exactly two alleles called in the file.
# Call rabhit to process the available anchors
# Assumes that the J-calls have been haplotyped

import argparse
import subprocess
import pathlib
import csv

# Read csv file into a list of dicts
def read_csv(file: str, delimiter: str = None):
    """Read a delimited file into a list of dicts (as produced by DictReader)

    :param file: filename of the file
    :type file: str
    :param delimiter: the delimiter (',' by default)
    :type delimiter: str
    :return: the list of dicts
    :rtype: list
    """
    ret = []
    with open(file, 'r') as fi:
        if delimiter:
            reader = csv.DictReader(fi, delimiter=delimiter)
        else:
            reader = csv.DictReader(fi)
        for row in reader:
            ret.append(row)

    return ret



def main():
    parser = argparse.ArgumentParser(description='Call rabhit to process available J anchors called in the sample')
    parser.add_argument('airr_file', help='Annotated sequences in AIRR format (tsv)')
    parser.add_argument('out_prefix', help='Prefix for output files')
    parser.add_argument('chain', help='chain (e.g. IGH, IGK, TRA, TRB..')
    parser.add_argument('-v', '--v_ref', help='V germline file, IMGT gapped (FASTA)')
    parser.add_argument('-d', '--d_ref', help='D germline file(FASTA)')
    args = parser.parse_args()

    python_dir = pathlib.Path(__file__).parent.resolve()
    recs = read_csv(args.airr_file, delimiter='\t')

    j_calls = set()

    for rec in recs:
        if ',' not in rec['j_call']:
            j_calls.add(rec['j_call'])

    j_genes = {}

    for call in j_calls:
        gene, allele = call.split('*')
        if gene not in j_genes:
            j_genes[gene] = []
        j_genes[gene].append(allele)

    for gene, alleles in j_genes.items():
        if len(alleles) == 2:
            alleles.sort()
            print(f"Rscript {python_dir}/to_rabhit.R {args.airr_file} {args.out_prefix} {args.chain} {gene} {alleles[0]} {alleles[1]} {args.v_ref} {args.d_ref if args.d_ref else ''}")

            if args.d_ref:
                subprocess.run([
                    "Rscript",
                    f"{python_dir}/to_rabhit.R",
                    f"{args.airr_file}",
                    f"{args.out_prefix}",
                    f"{args.chain}",
                    f"{gene}",
                    f"{alleles[0]}",
                    f"{alleles[1]}",
                    f"{args.v_ref}",
                    f"{args.d_ref}",
                ], capture_output=False)
            else:
                subprocess.run([
                    "Rscript",
                    f"{python_dir}/to_rabhit.R",
                    f"{args.airr_file}",
                    f"{args.out_prefix}",
                    f"{args.chain}",
                    f"{gene}",
                    f"{alleles[0]}",
                    f"{alleles[1]}",
                    f"{args.v_ref}",
                    f"{args.d_ref}",
                ], capture_output=False)


if __name__ == "__main__":
    main()

