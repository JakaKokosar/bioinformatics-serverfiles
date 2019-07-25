""" Gene sets update script """
import os
import io
import tempfile
from zipfile import ZipFile
from typing import Optional
from urllib.request import urlopen

from orangecontrib.bioinformatics import go, kegg
from orangecontrib.bioinformatics.omim import OMIM, diseases, disease_genes
from orangecontrib.bioinformatics.dicty import phenotypes
from orangecontrib.bioinformatics.kegg import caching
from orangecontrib.bioinformatics.ncbi import taxonomy
from orangecontrib.bioinformatics.ncbi.gene import GeneMatcher
from orangecontrib.bioinformatics.geneset import (
    GeneSet,
    GeneSets,
    filename
)


data_path = os.path.abspath(os.path.join(__file__, '../../data'))


def go_gene_sets(tax_id: str) -> None:
    domain = 'go'
    ontology = go.Ontology(filename=f'{data_path}/{domain}/gene_ontology.obo')
    annotations = go.Annotations(tax_id, filename=f'{data_path}/{domain}/{tax_id}.tab', ontology=ontology)

    def to_gene_set(term: go.Term) -> Optional[GeneSet]:
        genes = annotations.get_genes_by_go_term(term.id)

        if len(genes) > 0:
            return GeneSet(
                gs_id=term.id,
                name=term.name,
                genes=genes,
                hierarchy=('GO', term.namespace),
                organism=tax_id,
                link=f'http://amigo.geneontology.org/amigo/term/{term.id}'
            )

    gene_sets = GeneSets([gs for gs in [to_gene_set(term) for term in ontology.terms.values()] if gs is not None])

    for gs_group in gene_sets.split_by_hierarchy():
        hierarchy = gs_group.common_hierarchy()
        gs_group.to_gmt_file_format(f'{data_path}/gene_sets/{filename(hierarchy, tax_id)}')


def dicty_mutant_gene_sets(org):
    """ Return dicty mutant phenotype gene sets from Dictybase
    """
    if org == '352472':
        gene_sets = []
        gene_matcher = GeneMatcher('352472')

        for phenotype, mutants in phenotypes.phenotype_mutants().items():
            phenotype = phenotype.replace(",", " ")
            gene_symbols = [phenotypes.mutant_genes(mutant)[0] for mutant in mutants]
            gene_matcher.genes = gene_symbols
            genes = []

            for gene in gene_matcher.genes:
                if gene.gene_id is not None:
                    genes.append(int(gene.gene_id))

            if len(gene_symbols) != len(genes):
                print(len(gene_symbols), len(genes))

            gs = GeneSet(gs_id=phenotype,
                         name=phenotype,
                         genes=genes,
                         hierarchy=('Dictybase', 'Phenotypes'),
                         organism='352472',
                         link='')

            gene_sets.append(gs)

        return GeneSets(gene_sets)


def kegg_gene_sets(tax_id: str) -> None:
    """ Returns gene sets from KEGG pathways.
    """
    caching.clear_cache()
    kegg_org = kegg.KEGGOrganism(taxonomy.name(tax_id))
    ncbi_id_mapper = kegg_org.kegg_to_ncbi_mapper()
    genesets = []

    for id in kegg_org.pathways():
        pway = kegg.KEGGPathway(id)
        hier = ('KEGG', 'Pathways')

        if pway.pathway_attributes():
            kegg_names = kegg_org.get_genes_by_pathway(id)
            mapped_genes = []
            for gene in kegg_names:
                try:
                    mapped_genes.append(ncbi_id_mapper[gene.upper()])
                except KeyError:
                    # some kegg names can not be matched to ncbi ids
                    # they are included in geneset anyway
                    # remove prefix, that specifies kegg organism
                    # mapped_genes.append(gene.split(':')[-1])
                    pass

            gs = GeneSet(gs_id=id,
                         name=pway.title,
                         genes=mapped_genes,
                         hierarchy=hier,
                         organism=tax_id,
                         link=pway.link)
            genesets.append(gs)

    for gs_group in GeneSets(genesets).split_by_hierarchy():
        hierarchy = gs_group.common_hierarchy()
        gs_group.to_gmt_file_format(f'{data_path}/gene_sets/{filename(hierarchy, tax_id)}')


def cytoband_gene_sets(tax_id: str) -> None:
    """ Create cytoband gene sets from Stanford Microarray Database
    """
    if tax_id == '9606':
        download_link = 'http://statweb.stanford.edu/~tibs/GSA/cytobands-stanford.gmt'
        gene_matcher = GeneMatcher('9606')

        with urlopen(download_link) as stream:
            data = stream.read().splitlines()
            genesets = []

            for band in data:
                b = band.decode().split('\t')
                gene_symbols = b[2:]
                gene_matcher.genes = gene_symbols

                genes = []
                for gene in gene_matcher.genes:
                    if gene.gene_id is not None:
                        genes.append(gene.gene_id)

                genesets.append(GeneSet(gs_id=b[0], name=b[1], genes=genes if b[2:] else [],
                                        hierarchy=('Cytobands',), organism='9606', link=''))

        for gs_group in GeneSets(genesets).split_by_hierarchy():
            hierarchy = gs_group.common_hierarchy()
            gs_group.to_gmt_file_format(f'{data_path}/gene_sets/{filename(hierarchy, tax_id)}')


def reactome_gene_sets(tax_id: str) -> None:
    """ Prepare human pathways gene sets from reactome pathways
    """
    if tax_id == '9606':
        download_link = 'http://www.reactome.org/download/current/ReactomePathways.gmt.zip'
        file_name = 'ReactomePathways.gmt'
        detail_link = 'https://reactome.org/content/detail/{}'

        gene_matcher = GeneMatcher('9606')

        with urlopen(download_link) as url:
            memfile = io.BytesIO(url.read())

            with ZipFile(memfile, 'r') as myzip:
                f = myzip.open(file_name)
                content = f.read().decode().splitlines()
                genesets = []

                for path in content:
                    gene_symbols = path.split('\t')[2:] if path.split('\t')[2:] else []
                    gene_matcher.genes = gene_symbols
                    genes = []

                    for gene in gene_matcher.genes:
                        if gene.gene_id is not None:
                            genes.append(str(gene.gene_id))

                    pathway = path.split('\t')[0].replace(',', ' ')
                    pathway_id = path.split('\t')[1].replace(',', ' ')

                    gs = GeneSet(gs_id=pathway_id,
                                 name=pathway,
                                 genes=genes,
                                 hierarchy=('Reactome', 'pathways'),
                                 organism='9606',
                                 link=detail_link.format(pathway_id))

                    genesets.append(gs)

        for gs_group in GeneSets(genesets).split_by_hierarchy():
            hierarchy = gs_group.common_hierarchy()
            gs_group.to_gmt_file_format(f'{data_path}/gene_sets/{filename(hierarchy, tax_id)}')


if __name__ == "__main__":

    for common_tax_id in taxonomy.common_taxids():
        reactome_gene_sets(common_tax_id)
        cytoband_gene_sets(common_tax_id)
        # dicty_mutant_gene_sets(common_tax_id) TODO: re-enable this after orange3-bioinformatics release

        try:
            kegg_gene_sets(common_tax_id)
        except taxonomy.utils.UnknownSpeciesIdentifier as e:
            # KEGG organism code not found
            pass
        try:
            go_gene_sets(common_tax_id)
        except FileNotFoundError as e:
            # Organism is not supported in Gene Ontology module
            pass
