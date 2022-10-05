import os
import gzip
import json
import shutil
import tarfile
import serverfiles

from http.server import HTTPServer, SimpleHTTPRequestHandler
from multiprocessing import Process
from datetime import datetime
from collections import OrderedDict
from typing import Generator, List
from pathlib import Path, PosixPath

from plumbum import local

from orangecontrib.bioinformatics.ncbi.taxonomy import Taxonomy, common_taxid_to_name, shortname
from orangecontrib.bioinformatics.ncbi.taxonomy.utils import TaxonomyDB
from orangecontrib.bioinformatics.geneset import filename_parse


taxonomy_db = Taxonomy()

sqlite_utils_insert = local['sqlite_utils']['insert']
sqlite_utils_fts = local['sqlite_utils']['enable-fts']
sqlite_utils_optimize = local['sqlite_utils']['optimize']


def list_files(path: str) -> Generator[PosixPath, None, None]:
    path = Path(path)
    for file in (entry for entry in path.iterdir() if entry.is_file() and not entry.name.startswith('.')):
        yield file


class ServerFile:

    info_file_schema = {
        'domain': None,
        'filename': None,
        'source': None,
        'title': None,
        'tags': [],
        'size': None,
        'datetime': None,
        # used only if files are compressed
        'uncompressed': None,
        'compression': None,
    }

    def __init__(self, domain: str, file_path: str, title: str, description: str, tags: List[str], compression='gz'):
        self.domain = domain
        self.file_path = file_path
        self.file_name = os.path.basename(file_path)
        self.title = title
        self.description = description
        self.tags = tags
        self.compression = compression

    def create_info_file(self, **kwargs):
        info_dict = OrderedDict(self.info_file_schema)

        info_dict.update(**kwargs)
        info_dict['datetime'] = '{0:%Y-%m-%d %H:%M:%S.%f}'.format(datetime.today())
        info_dict['size'] = os.stat(self.file_path).st_size
        info_dict['source'] = 'server_file'

        with open(self.file_path + '.info', 'wt') as f:
            json.dump(info_dict, f)

    def to_server_format(self):
        uncompressed_size = os.stat(self.file_path).st_size

        with open(self.file_path, 'rb') as f_in:
            with gzip.open(self.file_path + '.temp', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

        os.remove(self.file_path)
        os.rename(self.file_path + '.temp', self.file_path)

        self.create_info_file(domain=self.domain,
                              filename=self.file_name,
                              title=self.title,
                              tags=self.tags,
                              uncompressed=uncompressed_size,
                              compression=self.compression)


def prepare_gene_serverfiles():
    Path('serverfiles/gene/').mkdir(parents=True, exist_ok=True)

    for file in list_files('data/gene/'):
        tax_id = file.name.replace('.json', '')
        species_name = taxonomy_db.get_entry(tax_id).name

        # dump json to sqlite
        sqlite_dump_path = f"serverfiles/gene/{tax_id}.sqlite"
        sqlite_utils_insert(sqlite_dump_path, 'gene_info', file.absolute(), '--pk=gene_id')

        # enable full text search (FTS-5)
        sqlite_utils_fts(sqlite_dump_path, 'gene_info', ['gene_id', 'symbol', 'synonyms', 'db_refs', 'description',
                                                         'locus_tag', 'symbol_from_nomenclature_authority'])
        # optimize sqlite database
        sqlite_utils_optimize(sqlite_dump_path)

        # preparse serverfiles
        server_file = ServerFile(domain='gene',
                                 file_path=sqlite_dump_path,
                                 title=f'Gene Info: {species_name}',
                                 description=f'Gene Info for {species_name} and all its strains.',
                                 tags=['NCBI', 'genes', 'info', 'gene info', species_name, tax_id])

        server_file.to_server_format()


def prepare_marker_genes_serverfiles():
    Path('serverfiles/marker_genes/').mkdir(parents=True, exist_ok=True)

    for file in list_files('data/marker_genes/'):
        source_path = str(file.absolute())
        dest_path = f'serverfiles/marker_genes/{file.name}'
        shutil.copy2(source_path, dest_path)

        source = 'Panglao' if 'panglao' in file.name else 'CellMarker'
        server_file = ServerFile(domain='marker_genes',
                                 file_path=dest_path,
                                 title=f'Marker Genes: {source}',
                                 description=f'Marker genes from {source} database.',
                                 tags=[source, 'genes', 'markers', 'marker genes'])

        server_file.to_server_format()


def prepare_go_serverfiles():
    Path('serverfiles/go/').mkdir(parents=True, exist_ok=True)

    for file in list_files('data/go/'):
        if file.name == 'gene_ontology.obo':
            title = 'Gene Ontology (GO)'
            tags = ['gene', 'ontology', 'GO']
            description = 'Basic version of the GO, filtered such that the graph is guaranteed to be ' \
                          'acyclic and annotations can be propagated up the graph.'
        else:
            tax_id = file.name.replace('.tab', '')
            species_name = common_taxid_to_name(tax_id)  # taxonomy_db.get_entry(tax_id).name
            title = f'GO Annotations for {species_name}'
            tags = ['gene', 'annotation', 'ontology', 'GO', tax_id, species_name]
            description = 'The gene association file ingested from GO Consortium members.'

        source_path = str(file.absolute())
        dest_path = f'serverfiles/go/{file.name}'
        shutil.copy2(source_path, dest_path)

        # preparse serverfiles
        server_file = ServerFile(domain='go',
                                 file_path=dest_path,
                                 title=title,
                                 description=description,
                                 tags=tags)

        server_file.to_server_format()


def prepare_homologene_serverfiles():
    Path('serverfiles/homologene/').mkdir(parents=True, exist_ok=True)

    for file in list_files('data/homologene/'):
        source_path = str(file.absolute())
        dest_path = f'serverfiles/homologene/{file.name}'
        shutil.copy2(source_path, dest_path)

        server_file = ServerFile(domain='homologene',
                                 file_path=dest_path,
                                 title='Homolog Genes',
                                 description='Homoloh genes from NCBI HomoloGene database.',
                                 tags=['NCBI', 'genes', 'homologs', 'homolog genes', 'homologene'])

        server_file.to_server_format()


def prepare_dictybase_serverfiles():
    Path('serverfiles/dictybase/').mkdir(parents=True, exist_ok=True)

    for file in list_files('data/dictybase/'):
        source_path = str(file.absolute())
        dest_path = f'serverfiles/dictybase/{file.name}'
        shutil.copy2(source_path, dest_path)

        server_file = ServerFile(domain='dictybase',
                                 file_path=dest_path,
                                 title='dictyBase mutant phenotypes',
                                 description='All curated mutants with phenotypes',
                                 tags=['Dictyostelium discoideum', 'mutant', 'dictyBase', 'phenotype'])

        server_file.to_server_format()


def prepare_gene_sets_serverfiles():
    Path('serverfiles/gene_sets/').mkdir(parents=True, exist_ok=True)

    for file in list_files('data/gene_sets/'):
        source_path = str(file.absolute())
        dest_path = f'serverfiles/gene_sets/{file.name}'
        shutil.copy2(source_path, dest_path)

        hierarchy, tax_id = filename_parse(file.name)
        title = f"Gene sets: {', '.join(hierarchy)} ({common_taxid_to_name(tax_id)})"
        tags = list(hierarchy) + ['gene sets', common_taxid_to_name(tax_id)] + shortname(tax_id)

        server_file = ServerFile(domain='gene_sets',
                                 file_path=dest_path,
                                 title=title,
                                 description='',
                                 tags=tags)

        server_file.to_server_format()


def prepare_taxonomy_serverfiles():
    Path('serverfiles/taxonomy/').mkdir(parents=True, exist_ok=True)

    TaxonomyDB.download('serverfiles/taxonomy/')
    TaxonomyDB.init_db('serverfiles/taxonomy/taxonomy.sqlite',
                       tarfile.open('serverfiles/taxonomy/taxdump.tar.gz'))

    server_file = ServerFile(domain='taxonomy',
                             file_path='serverfiles/taxonomy/taxonomy.sqlite',
                             title='NCBI Taxonomy database',
                             description='',
                             tags=['NCBI', 'taxonomy', 'organism', 'taxid'])

    server_file.to_server_format()
    Path('serverfiles/taxonomy/taxdump.tar.gz').unlink()


def start_server():
    class Handler(SimpleHTTPRequestHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory='serverfiles/', **kwargs)

    httpd = HTTPServer(('localhost', 8000), Handler)
    httpd.serve_forever()


if __name__ == '__main__':
    p = Process(target=start_server)
    p.start()

    prepare_gene_serverfiles()
    prepare_marker_genes_serverfiles()
    prepare_homologene_serverfiles()
    prepare_dictybase_serverfiles()
    prepare_go_serverfiles()
    prepare_gene_sets_serverfiles()
    prepare_taxonomy_serverfiles()

    info = list(serverfiles.ServerFiles(server='http://localhost:8000/').allinfo().items())
    with open('serverfiles/__INFO__', 'wt') as f:
        json.dump(info, f)

    shutil.make_archive('serverfiles_dump', 'gztar', root_dir='serverfiles/')
    p.terminate()
