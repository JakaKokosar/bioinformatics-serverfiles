""" Gene Info update script """
import gzip
import json
import sys
from collections import defaultdict

from orangecontrib.bioinformatics.ncbi.taxonomy import Taxonomy, common_taxids, common_taxid_to_name

# columns indexes
# ftp://ftp.ncbi.nlm.nih.gov/gene/README under "gene_info" section
(
    tax_id,
    gene_id,
    symbol,
    locus_tag,
    synonyms,
    db_refs,
    chromosome,
    map_location,
    description,
    type_of_gene,
    symbol_from_nomenclature_authority,
    full_name_from_nomenclature_authority,
    nomenclature_status,
    other_designations,
    modification_date,
) = range(15)


def pipe_delimited_to_list(value: str) -> list:
    return [val for val in value.split('|') if val and val != '-']


def parse_db_refs(value: str) -> dict:
    """ Parse source string from NCBI gene info.

    ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/README:

    pipe-delimited set of identifiers in other databases
    for this gene.  The unit of the set is database:value.
    Note that HGNC and MGI include 'HGNC' and 'MGI', respectively,
    in the value part of their identifier.

    Consequently, dbXrefs for these databases will appear like: HGNC:HGNC:1100
    This would be interpreted as database='HGNC', value='HGNC:1100'


    Args:
        value (str): string of gene sources

    Returns:
        :obj:`dict`: Keys are source names, values are source ids

    """
    external_ids = value.split('|')
    out_dict = {}

    for ref in external_ids:
        if ref != '-':
            source_dest, source_id = ref.split(':', 1)
            out_dict[source_dest] = source_id

    return out_dict


def load_homologs():
    file_name = 'homologene.tab'

    with open(f'data/homologene/{file_name}', 'r') as fp:
        _homologs = {l[2]: tuple(l) for l in [line.strip().split('\t') for line in fp.readlines()]}
        _homologs_by_group = defaultdict(list)

        for _, (hid, tax, gid) in _homologs.items():
            _homologs_by_group[hid].append((hid, tax, gid))

    return _homologs, _homologs_by_group


def gene_info_to_dict(gene_data: tuple):
    homology_group = homologs.get(str(gene_data[gene_id]), [None])[0]
    homolog_genes = {tax: gid for (_, tax, gid) in homologs_by_group.get(homology_group, [])
                     if homology_group and tax != gene_data[tax_id]}
    return {
        'species': common_taxid_to_name(to_species[gene_data[tax_id]]),
        'tax_id': gene_data[tax_id],
        'gene_id': gene_data[gene_id],
        'symbol': gene_data[symbol],
        'synonyms': pipe_delimited_to_list(gene_data[synonyms]),
        'db_refs': parse_db_refs(gene_data[db_refs]),
        'description': gene_data[description] if gene_data[description] != '-' else None,
        'locus_tag': gene_data[locus_tag] if gene_data[locus_tag] != '-' else None,
        'chromosome': gene_data[chromosome] if gene_data[chromosome] != '-' else None,
        'map_location': gene_data[map_location] if gene_data[map_location] != '-' else None,
        'type_of_gene': gene_data[type_of_gene] if gene_data[type_of_gene] != '-' else None,
        'symbol_from_nomenclature_authority':
            gene_data[symbol_from_nomenclature_authority]
            if gene_data[symbol_from_nomenclature_authority] != '-'
            else None,
        'full_name_from_nomenclature_authority':
            gene_data[full_name_from_nomenclature_authority]
            if gene_data[full_name_from_nomenclature_authority] != '-'
            else None,
        'nomenclature_status': gene_data[nomenclature_status] if gene_data[nomenclature_status] != '-' else None,
        'other_designations': pipe_delimited_to_list(gene_data[other_designations]),
        'modification_date': gene_data[modification_date],
        'homology_group_id': homology_group,
        'homologs': homolog_genes
    }


def load_gene_info(file_path: str) -> dict:
    with gzip.open(file_path, 'rb') as info_file:
        # skip header
        info_file.readline()

        # store lines in memory
        genes_by_species = defaultdict(list)
        for line in info_file:
            info = tuple(line.decode().strip().split('\t'))
            species = to_species.get(info[tax_id], None)

            if species is not None:
                genes_by_species[species].append(info)

    return genes_by_species


def to_json(species: str) -> None:
    data = [gene_info_to_dict(gene) for gene in gene_info.get(species, [])]

    with open(f'data/gene/{species}.json', 'w', encoding='utf-8') as fp:
        json.dump(data, fp, ensure_ascii=False, indent=2)  # separators=(',', ':')


if __name__ == "__main__":
    taxonomy_db = Taxonomy()
    supported_taxonomies = [[tax] + taxonomy_db.get_all_strains(tax) for tax in common_taxids()]
    to_species = {tax: taxonomy_db.get_species(tax) for strains in supported_taxonomies for tax in strains}

    gene_info = load_gene_info(sys.argv[1])
    homologs, homologs_by_group = load_homologs()

    for tax in common_taxids():
        to_json(tax)
