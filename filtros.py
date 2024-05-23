from Bio import Seq, pairwise2
from Bio.SubsMat.MatrixInfo import blosum62


def get_percent_matches(seq1_in, seq2_in):
    '''
    Performs a pairwise alignment between 2 sequences and returns the percent
    similarity 

    % Similarity is calculated as: #matches / (#matches + #mismataches); gaps are not included

    Parameters
    ----------
    seq1_str : string
        Sequence 1.
    seq2_str : string
        Sequence 2.
    seq_type : string, optional
        "dna" or "prot" for defining the sequence type. The default is "dna".

    Returns
    -------
    average_percent_similar : double
        The average % gapless matches between the set of optimal alignments
        generated for the two sequences.

    '''

    seq1 = Seq.Seq(seq1_in)
    seq2 = Seq.Seq(seq2_in)

    # Alignments will hold all global optimal alignments


        # Perform DNA alignment with the following scores:
        # Match = +2
        # Mismatch = -1
        # Gap-opening = -1
        # Gap-extending = -.5

    alignments = pairwise2.align.globalms(seq1, seq2, 2, -1, -1, -.5)


    # If there is more than one global optimal alignment, the % identity is
    # calculated for each alignemnt and averaged. The average is compared to
    # the threshold. 
    # Holds the percent matches for each alignment
    percent_matches = []

    for align in alignments:

        # The format of an alignment:
        # [0] : first sequence aligned
        # [1] : second sequenced aligned
        # [2] : alignment score
        # [3] : starting position of the alignment
        # [4] : ending position of the alignment

        seq1_aligned = align[0]
        seq2_aligned = align[1]
        matches = 0

        # For a global alignment, the start position is always 0
        align_size = align[4]
        size_adj = 0

        for x in range(align_size):

            # size_adj is subtracted from the length of the alignment to remove
            # gapped positions. 
            if seq1_aligned == "-" or seq2_aligned == "-":
                size_adj += 1
                continue

            if seq1_aligned[x] == seq2_aligned[x]:
                matches += 1

        # The size of the alignment is adjusted for the gapped positions
        current_percent_matches = (matches / (align_size - size_adj))

        percent_matches.append(current_percent_matches)
        # print("The percent match: ", current_percent_matches)

    # The average percent match is calculated for the ensemble of alignments
    # that was returned. 
    total_match_percentage = 0

    for match_percentage in percent_matches:
        total_match_percentage += match_percentage
    if len(percent_matches) != 0:
        average_percent_similar = total_match_percentage / len(percent_matches)
    else:
        average_percent_similar = 0

    # print("The average percent match is: ", average_percent_similar)

    return average_percent_similar
def sim_filter(input_seqs, percent_ident):

    print("|~> Filtering the sequences")

    seq_status = {seq: False for seq in input_seqs}
    filtered_list = []
    i = 0
    for sequencia in seq_status:
        i = i+1
        if not seq_status[sequencia]:
            for seq in seq_status[i+1:]:
                if not seq_status[seq]:
                    if get_percent_matches(str(sequencia), str(seq)) >= percent_ident:
                        seq_status[seq] = True

    for sequ in seq_status:
        if not seq_status[sequ]:
            filtered_list.append(sequ)

    return filtered_list
