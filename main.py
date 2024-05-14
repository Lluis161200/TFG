import json

import search_in_blast


def main():
    data = json.load(open("inicializacion.json"))
    blast_parameters = data["blast_parameters"]
    e_value = blast_parameters["e_value"]
    min_coverage = blast_parameters["coverage"]
    protein_list = data["protein_accession_list"]
    max_hits = 10
    taxonomic_id = data["taxonomic_id"]
    bd = blast_parameters["bd"]
    hits = []

    for protein in protein_list.values():
        for protein_accesion in protein:
            hits = search_in_blast.search_blast(protein_accesion, bd, taxonomic_id, e_value, max_hits, 3, min_coverage)


main()
