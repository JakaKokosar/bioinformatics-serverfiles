""" GO update script """
import csv
import gzip
import sys

from shutil import copy2
from collections import defaultdict

from orangecontrib.bioinformatics.ncbi.taxonomy import Taxonomy, common_taxids


def load_gene2go(file_path: str) -> dict:
    genes_by_species = defaultdict(list)

    with gzip.open(file_path, 'rb') as fp:
        # skip header
        fp.readline()

        # store lines in memory
        for line in fp:
            columns = line.decode().strip().split('\t')
            species = to_species.get(columns[0], None)

            if species is not None:
                genes_by_species[species].append(columns)

    return genes_by_species


def gene_annotations(species: str) -> None:
    header = ['tax_id', 'GeneID', 'GO_ID', 'Evidence', 'Qualifier', 'GO_term' 'PubMed', 'Category']
    data = gene2go.get(species, None)

    if data is not None:
        with open(f'data/go/{species}.tab', 'w') as fp:
            csv_writer = csv.writer(fp, delimiter='\t')
            csv_writer.writerow(header)
            csv_writer.writerows(data)


def gene_ontology(file_path: str) -> None:
    copy2(file_path, f'data/go/gene_ontology.obo')


if __name__ == "__main__":
    taxonomy_db = Taxonomy()
    supported_taxonomies = [[tax] + taxonomy_db.get_all_strains(tax) for tax in common_taxids()]
    to_species = {tax: taxonomy_db.get_species(tax) for strains in supported_taxonomies for tax in strains}

    gene2go = load_gene2go(sys.argv[1])
    for tax in common_taxids():
        gene_annotations(tax)

    gene_ontology(sys.argv[2])
