""" HomoloGene update script """
import csv
import gzip
import sys
from typing import Dict, Tuple

from orangecontrib.bioinformatics.ncbi.taxonomy import common_taxids


def gene_history(file_path: str):
    taxonomy_col = 0
    current_id_col = 1
    discontinued_id_col = 2

    with gzip.open(file_path, 'rb') as fp:
        # skip header
        fp.readline()

        # store lines in memory
        out_dict = {}
        for line in fp:
            info = tuple(line.decode().strip().split('\t'))
            out_dict[(info[taxonomy_col], info[discontinued_id_col])] = info[current_id_col]

    return out_dict


def homologene(file_path: str, history: Dict[Tuple[str, str], str]):
    file_name = 'homologene.tab'
    common_tax_ids = common_taxids()

    with open(file_path, 'r') as f:
        content = f.read().strip()

    group_id = 0
    tax_id = 1
    entrez_id = 2

    homologs_out = [['group_id', 'tax_id', 'entrez_id']]
    for line in content.split('\n'):
        columns = line.split('\t')

        # Note: if '-' -> discontinued
        #       if 'some id' -> new entrez id
        #       if None -> valid entrez id
        gene_id_status = history.get((columns[tax_id], columns[entrez_id]), None)

        if not columns[tax_id] in common_tax_ids:
            continue
        elif gene_id_status == '-':
            continue

        if gene_id_status is None:
            gene_id = columns[entrez_id]
        else:
            gene_id = gene_id_status

        homologs_out.append([columns[group_id], columns[tax_id], gene_id])

    with open(f'data/homologene/{file_name}', 'w') as fp:
        csv_writer = csv.writer(fp, delimiter='\t')
        csv_writer.writerows(homologs_out)


if __name__ == "__main__":
    homologs_fp = sys.argv[1]
    gene_history_fp = sys.argv[2]

    homologene(homologs_fp, gene_history(gene_history_fp))
