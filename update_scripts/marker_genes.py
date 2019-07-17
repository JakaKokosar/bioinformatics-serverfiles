""" Marker genes update script """
import gzip
import sys
from collections import defaultdict

from Orange.data import Table, Domain, StringVariable
from orangecontrib.bioinformatics.ncbi.gene import GeneMatcher, Gene


def panglao_db(file_path: str):
    file_name = 'panglao_gene_markers.tab'
    reference, reference_url = 'PanglaoDB', 'https://panglaodb.se/'

    with gzip.open(file_path, 'rb') as f:
        content = f.read().decode('utf-8').strip()

    species = 0
    gene_symbol = 1
    cell_type = 2
    genes_by_organism = defaultdict(list)
    organism_mapper = {'Mm': 'Mouse', 'Hs': 'Human'}

    def _gene_function_table(desc_col: StringVariable, gm_results: GeneMatcher):
        _domain = Domain([], metas=[desc_col])
        _data = [[str(gene.description) if gene.description else ''] for gene in gm_results.genes]
        return Table(_domain, _data)

    for line in content.split('\n'):
        columns = line.split('\t')

        for org in columns[species].split(' '):
            if org in organism_mapper.keys():
                gene_entry = [organism_mapper[org], columns[gene_symbol], columns[cell_type], reference, reference_url]
                genes_by_organism[organism_mapper[org]].append(gene_entry)

    domain = Domain(
        [],
        metas=[
            StringVariable('Organism'),
            StringVariable('Name'),
            StringVariable('Cell Type'),
            StringVariable('Reference'),
            StringVariable('URL'),
        ],
    )

    entrez_id_column = StringVariable('Entrez ID')
    description_column = StringVariable('Function')

    # construct data table for mouse
    gm_mouse = GeneMatcher('10090', case_insensitive=True)
    mouse_table = Table(domain, genes_by_organism['Mouse'])
    mouse_table = gm_mouse.match_table_column(mouse_table, 'Name', entrez_id_column)
    mouse_table = Table.concatenate([mouse_table, _gene_function_table(description_column, gm_mouse)])

    # construct data table for human
    gm_human = GeneMatcher('9606', case_insensitive=True)
    human_table = Table(domain, genes_by_organism['Human'])
    human_table = gm_human.match_table_column(human_table, 'Name', entrez_id_column)
    human_table = Table.concatenate([human_table, _gene_function_table(description_column, gm_human)])

    # return combined tables
    Table.concatenate([mouse_table, human_table], axis=0).save(f'data/marker_genes/{file_name}')


def cell_marker_db(file_path: str):
    file_name = 'cellMarker_gene_markers.tab'
    reference, reference_url = '{}', 'https://www.ncbi.nlm.nih.gov/pubmed/?term={}'

    with open(file_path, 'r') as f:
        content = f.read().strip()

    species_col = 0
    cell_name_col = 5
    gene_symbol_col = 8
    entrez_id_col = 9
    pubmed_id_col = 13
    organisms = ['Human', 'Mouse']
    data = list()

    domain = Domain(
        [],
        metas=[
            StringVariable('Organism'),
            StringVariable('Name'),
            StringVariable('Cell Type'),
            StringVariable('Reference'),
            StringVariable('URL'),
            StringVariable('Entrez ID'),
        ],
    )

    unique_rows = set()

    for line in content.split('\n'):
        columns = line.split('\t')
        organism = columns[species_col]

        if organism in organisms:
            symbols = columns[gene_symbol_col].replace('[', '').replace(']', '').split(', ')
            entrez_ids = columns[entrez_id_col].replace('[', '').replace(']', '').split(', ')

            for symbol, entrez_id in zip(symbols, entrez_ids):
                try:
                    int(entrez_id)
                except ValueError:
                    continue

                if (columns[cell_name_col], entrez_id) in unique_rows:
                    continue

                ref = reference.format(columns[pubmed_id_col])

                try:
                    int(ref)
                    ref_link = reference_url.format(columns[pubmed_id_col])
                except ValueError:
                    # its not pubmed_id
                    ref_link = '?'

                unique_rows.add((columns[cell_name_col], entrez_id))
                gene_entry = [organism, symbol, columns[cell_name_col], ref, ref_link, entrez_id]

                data.append(gene_entry)

    table = Table(domain, data)
    genes = [Gene(name) for name in table.get_column_view('Entrez ID')[0]]
    for gene in genes:
        gene.ncbi_id = gene.input_name
        gene.load_ncbi_info()

    description_column = StringVariable('Function')
    domain = Domain([], metas=table.domain.metas + (description_column,))
    table = table.transform(domain)
    table[:, description_column] = [[gene.description] for gene in genes]
    table.save(f'data/marker_genes/{file_name}')


if __name__ == "__main__":
    panglao_db(sys.argv[1])
    cell_marker_db(sys.argv[2])
