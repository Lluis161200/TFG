

from Bio import Align
from Bio.Seq import Seq

def get_percent_matches(seq1_in, seq2_in):
    aligner = Align.PairwiseAligner()
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -1
    aligner.target_end_gap_score = 0.0
    aligner.query_end_gap_score = 0.0

    # Perform global alignment
    aligner.mode = 'global'
    alignments = aligner.align(seq1_in, seq2_in)

    percent_matches = []

    for align in alignments:
        seq1_aligned = align.seqA
        seq2_aligned = align.seqB
        matches = 0

        align_size = len(seq1_aligned)
        size_adj = 0

        for x in range(align_size):
            # gapped positions
            if seq1_aligned[x] == "-" or seq2_aligned[x] == "-":
                size_adj += 1
                continue
            # matches
            if seq1_aligned[x] == seq2_aligned[x]:
                matches += 1

        # Exclude gaps
        current_percent_matches = (matches / (align_size - size_adj))
        percent_matches.append(current_percent_matches)

    # Calculate average of all alignments
    if percent_matches:
        average_percent_similar = sum(percent_matches) / len(percent_matches)
    else:
        average_percent_similar = 0

    print("The average percent match is:", average_percent_similar)
    return average_percent_similar


def sim_filter(input_seqs, percent_ident):
    print("|~> Filtering the sequences")
    print("Sequencia Inicial", input_seqs)

    # Create a dictionary to store the status of each sequence using the sequence string as the key
    seq_status = {str(seq.seq): False for seq in input_seqs}
    filtered_list = []
    seq_list = list(seq_status.keys())

    for i, seq_str in enumerate(seq_list):
        if not seq_status[seq_str]:
            for other_seq_str in seq_list[i+1:]:
                if not seq_status[other_seq_str]:
                    if get_percent_matches(seq_str, other_seq_str) >= percent_ident:
                        print("Sequencia eliminada:", other_seq_str)
                        seq_status[other_seq_str] = True

    for seq_str in seq_status:
        if not seq_status[seq_str]:
            # Find the original SeqRecord corresponding to the sequence string
            for seq in input_seqs:
                if str(seq.seq) == seq_str:
                    filtered_list.append(seq)
                    break

    return filtered_list