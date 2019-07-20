import re
import requests

from plumbum import local
from plumbum.cmd import wget, FG


def fetch_html(url: str) -> str:
    # todo: handle exceptions
    return requests.get(url).content.decode('utf-8')


def panglao_db_filename() -> str:
    panglao_url = 'https://panglaodb.se/markers.html'
    re_pattern = r'<a href=\"markers/(.*?)\"'
    m = re.search(re_pattern, fetch_html(panglao_url))
    return m.group(1)


panglao_fn = panglao_db_filename()

# Download links
gene_info = 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz'
gene_history = 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_history.gz'
gene_to_go = 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz'
homolog_genes = 'ftp://ftp.ncbi.nlm.nih.gov/pub/HomoloGene/current/homologene.data'
gene_ontology = 'http://purl.obolibrary.org/obo/go/go-basic.obo'
cellmarker = 'http://biocc.hrbmu.edu.cn/CellMarker/download/all_cell_markers.txt'
panglaodb = f"https://panglaodb.se/markers/{panglao_fn}"

# dictybase
base_url = 'http://dictybase.org/db/cgi-bin/dictyBase/download/download.pl?area=mutant_phenotypes&ID={}'
all_mutants = 'all-mutants.txt'
null_mutants = 'null-mutants.txt'
overexpression_mutants = 'overexpression-mutants.txt'
developmental_mutants = 'developmental-mutants.txt'
multiple_mutants = 'multiple-mutants.txt'
other_mutants = 'other-mutants.txt'


python = local['python']
download_process = wget['-P', 'temp/']

print(f'downloading gene info from {gene_info}')
download_process(gene_info)

print(f'downloading gene history from {gene_history}')
download_process(gene_history)

print(f'downloading gene to go from {gene_to_go}')
download_process(gene_to_go)

print(f'downloading homolog genes from {homolog_genes}')
download_process(homolog_genes)

print(f'downloading gene ontology from {gene_ontology}')
download_process(gene_ontology)

print(f'downloading cellmarker gene markers from {cellmarker}')
download_process(cellmarker)

print(f'downloading panglao gene markers from {panglaodb}')
download_process(panglaodb)


print('Download mutants from dictybase ...')
download_process['-O', f'temp/{all_mutants}'](base_url.format(all_mutants))
download_process['-O', f'temp/{null_mutants}'](base_url.format(null_mutants))
download_process['-O', f'temp/{overexpression_mutants}'](base_url.format(overexpression_mutants))
download_process['-O', f'temp/{developmental_mutants}'](base_url.format(developmental_mutants))
download_process['-O', f'temp/{multiple_mutants}'](base_url.format(multiple_mutants))
download_process['-O', f'temp/{other_mutants}'](base_url.format(other_mutants))

python('update_scripts/gene.py', 'temp/gene_info.gz')
python('update_scripts/marker_genes.py', f'temp/{panglao_fn}', 'temp/all_cell_markers.txt')
python('update_scripts/homologene.py', 'temp/homologene.data', 'temp/gene_history.gz')
python('update_scripts/go.py', 'temp/gene2go.gz', 'temp/go-basic.obo')
python(
    'update_scripts/dictybase.py',
    f'temp/{all_mutants}',
    f'temp/{null_mutants}',
    f'temp/{overexpression_mutants}',
    f'temp/{multiple_mutants}',
    f'temp/{developmental_mutants}',
    f'temp/{other_mutants}',
)
python('update_scripts/gene_sets.py')
