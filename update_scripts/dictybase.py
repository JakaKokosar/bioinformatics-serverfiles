import sys
import json

from typing import List, Dict, Any


def mutant_to_dict(mutant: List[str], mutant_types: Dict[str, bool]) -> Dict[str, Any]:
    return {
        'systematic_name': mutant[0].strip(),
        'strain_descriptor': mutant[1].strip(),
        'associated_genes': mutant[2].split(' | '),
        'phenotypes': mutant[3].split(' | '),
        'mutant_types': mutant_types,
    }


def phenotypes_for_dictyostelium_mutants(
    all_mutants: str,
    null_mutants: str,
    overexpression_mutants: str,
    multiple_mutants: str,
    developmental_mutants: str,
    other_mutants: str,
):

    file_name = 'mutant_phenotypes.json'

    with open(null_mutants, 'r') as fp:
        null_mutants = set([line.strip().split('\t')[0] for line in fp])

    with open(overexpression_mutants, 'r') as fp:
        overexpression_mutants = set([line.strip().split('\t')[0] for line in fp])

    with open(multiple_mutants, 'r') as fp:
        multiple_mutants = set([line.strip().split('\t')[0] for line in fp])

    with open(developmental_mutants, 'r') as fp:
        developmental_mutants = set([line.strip().split('\t')[0] for line in fp])

    with open(other_mutants, 'r') as fp:
        other_mutants = set([line.strip().split('\t')[0] for line in fp])

    with open(all_mutants, 'r') as fp:
        # skip header
        fp.readline()

        all_mutants = []
        for line in fp.readlines():
            mutant = line.strip().split('\t')

            if len(mutant) == 4:
                mutant_id = mutant[0]

                types = {
                    'null': mutant_id in null_mutants,
                    'overexpression': mutant_id in overexpression_mutants,
                    'multiple ': mutant_id in multiple_mutants,
                    'other': mutant_id in other_mutants,
                    'developmental': mutant_id in developmental_mutants,
                }

                all_mutants.append(mutant_to_dict(mutant, types))

    with open(f'data/dictybase/{file_name}', 'w', encoding='utf-8') as fp:
        json.dump(all_mutants, fp, ensure_ascii=False, indent=4)


if __name__ == "__main__":
    phenotypes_for_dictyostelium_mutants(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
