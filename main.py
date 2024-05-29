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
    max_hits = 50
    taxonomic_id = data["taxonomic_id"]
    bd = blast_parameters["bd"]
    min_len = 10
    max_percent_similarity = 0.75
    file_name = "prueba"
    file_path = "blast_results.xml"

    # Initialize filtered_seqs to store sequences from all iterations
    filtered_seqs = []

    for protein in protein_list.values():
        for protein_accession in protein:
            hits = search_in_blast.search_blast(protein_accession, bd, taxonomic_id, e_value, max_hits, 3, min_coverage)
            #hits = search_in_blast.parse_blast_xml(protein_accession, file_path, min_coverage)
            nucleotide_records = []
            for hit in hits:
                nucleotide_records.append((hit[0], hit[1], get_nuc_rec_from_prot(hit[1])))

            nuc_seq_rec = []

            # For each record, fetch the sequence record that is requested
            for record in nucleotide_records:
                # Get the sequence record requested
                seq_record = get_nuc_seq(record)
                print(seq_record)
                if seq_record is not None:
                    seq_length = len(seq_record.seq)
                    if seq_length > min_len:
                        nuc_seq_rec.append(seq_record)
                    else:
                        print("|~> Record did not meet min_len requirement")

            # Filter sequences and append to the filtered_seqs list
            filtered_seqs.extend(sim_filter(nuc_seq_rec, max_percent_similarity))

    # Write all filtered sequences to the output file
    o_file = file_name + ".fasta"
    print("|~> Writing to output file")
    SeqIO.write(filtered_seqs, o_file, "fasta")

main()
