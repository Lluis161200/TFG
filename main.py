import json

from Bio import SeqIO

import search_in_blast
from filtros import sim_filter
from get_nucleotide import get_nuc_rec_from_prot, get_nuc_seq


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
    min_len = 10
    max_percent_similarity = 0.75
    file_name = "prueba"

    for protein in protein_list.values():
        for protein_accesion in protein:
            hits = search_in_blast.search_blast(protein_accesion, bd, taxonomic_id, e_value, max_hits, 3, min_coverage)
            nucleotide_records = []
            for hit in hits:
                nucleotide_records.append((hit[0], hit[1], get_nuc_rec_from_prot(hit[1])))

            nuc_seq_rec = []

            # For each record, fetch  the sequence record that is requested
            for record in nucleotide_records:
                # Get the sequence record requested
                seq = get_nuc_seq(record, 250, 50,)
                if len(seq) > min_len:
                    nuc_seq_rec.append(seq)
                else:
                    print("|~> Record did not meet min_len requirement")

            filtered_seqs = sim_filter(nuc_seq_rec, max_percent_similarity)
            o_file = file_name + ".fasta"
            print("|~> Writing to output file")
            SeqIO.write(filtered_seqs, o_file, "fasta")


main()
