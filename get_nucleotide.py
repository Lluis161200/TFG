import time

from Bio import Entrez, SeqIO, Seq
from Bio.SeqRecord import SeqRecord

REQUEST_LIMIT = 5
SLEEP_TIME = .5

def get_nuc_rec_from_prot(prot_id):
    '''
    Gets the the most complete nucelotide record for a given protein accession.

    Parameters
    ----------
    prot_id : string
        Accession number for the protein record of interest.

    Returns
    -------
    max_p_record : string
        The accession for the nucleotide record with the most complete genome.

    '''

    print("|~> Getting nucelotide records for " + str(prot_id))

    # Fetch the protein record
    print("\t|~> Fetching the protein records")
    max_p_record=0

    for i in range(REQUEST_LIMIT):
        try:
            records = Entrez.read(Entrez.efetch(db="protein", id=prot_id, rettype="ipg",
                                                retmode="xml"))
            time.sleep(SLEEP_TIME)
            break

        except:

            print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")

            if i == (REQUEST_LIMIT - 1):
                print("\t\tCould not download record after " + str(REQUEST_LIMIT) + " attempts")
                print("\t|~> No gene records found for " + prot_id)
                return None

    # The priority scores for the types of gene records available
    p_scores = {"NC_": 7, "AC_": 7,
                "AE": 6, "CP": 6, "CY": 6,
                "NZ_": 5, "NT_": 5, "NW_": 5,
                "AAAA-AZZZ": 4,
                "U": 3, "AF": 3, "AY": 3, "DQ": 3}

    # Gene records are appended to this list
    genome_list = []

    print("\t|~> Recording the gene records with priorities")
    if 'ProteinList' in records['IPGReport'].keys():
        for idprotein in records['IPGReport']['ProteinList']:
            if 'CDSList' in idprotein.keys():
                for cds in idprotein['CDSList']:
                    cds_acc = cds.attributes['accver']
                    cds_start = cds.attributes['start']
                    cds_stop = cds.attributes['stop']
                    cds_strand = cds.attributes['strand']
                    cds_scr = 0
                    # assign priority
                    for key in p_scores:
                        if cds_acc.startswith(key):
                            cds_scr = p_scores[key]
                            if max_p_record < cds_scr:
                                max_p_record = cds_scr
                                # create and append record
                                cds_rec = {'acc': cds_acc, 'start': cds_start,
                                           'stop': cds_stop, 'strand': cds_strand,
                                           'p_score': cds_scr}

            else:
                continue
    else:
        print("\t|~> No gene records found for " + prot_id)
        return None

    # Finds the genome with the max p-score
    if bool(cds_rec):
        return cds_rec
    else:
        print("\t|~> No gene records found for " + prot_id)
        return None


def get_nuc_seq(nuc_rec_in, start_adj=250, stop_adj=3, isolate_promoters=False):
    '''
    Returns the nucleotide sequence of nucelotide record.

    Parameters
    ----------
    nuc_rec_in : list[(string, string, string)]
        First element in the tuple is the input protein record the hit came from.
        The second element is the protein accession from the BLAST search of the input protein.
        The third element is the nucleotide record for the hit.
    start_adj : int, optional
        The number of nucleotides to go upstream of the start site. The default is 250.
    stop_adj : int, optional
        The number of nucelotides to go downstream of the start site. The default is 3.
    isolate_promoters : bool, optional
        Set whether to return promoter sequences or the raw sequence. The default is False.

    Returns
    -------
    SeqRecord
        Holds the promoter or nucleotide sequence.

    '''

    nuc_rec = nuc_rec_in[2]

    if nuc_rec is None:
        print("\t|~>Record is not valid")
        return None

    print("|~> Get nucleotide sequence for " + str(nuc_rec['acc']))

    print("\t|~> Adjusting start and stop positions")

    if nuc_rec['strand'] == '+':
        s_start = int(nuc_rec['start']) - start_adj
        s_stop = int(nuc_rec['start']) + stop_adj
        s_strand = 1
    else:
        s_stop = int(nuc_rec['stop']) + start_adj
        s_start = int(nuc_rec['stop']) - stop_adj
        s_strand = 2

    if isolate_promoters:

        print("\t|~> Getting genbank record")

        # Fetch and read the annotated GenBank record

        for i in range(REQUEST_LIMIT):

            try:

                handle = Entrez.efetch(db="nuccore", id=nuc_rec['acc'], strand=s_strand,
                                       seq_start=s_start, seq_stop=s_stop,
                                       rettype='gbwithparts', retmode="XML")

                genome_record = Entrez.read(handle, "xml")

                time.sleep(SLEEP_TIME)
                break

            except:

                print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")

                if i == (REQUEST_LIMIT - 1):
                    print("\t\tCould not download record after " + str(REQUEST_LIMIT) + " attempts")

        print("\t|~> Parsing intervals for coding regions")
        # Find all coding regions in the returned GenBank sequence.
        coding_intervals = []

        sequence = genome_record[0]['GBSeq_sequence']

        for feature in genome_record[0]['GBSeq_feature-table']:
            if feature['GBFeature_key'] == 'gene':
                if "GBInterval_from" in feature['GBFeature_intervals'][0]:
                    coding_start = feature['GBFeature_intervals'][0]['GBInterval_from']
                    coding_end = feature['GBFeature_intervals'][0]['GBInterval_to']
                    coding_intervals.append((coding_start, coding_end))

        # The FASTA ID for the promoter sequence is in the following format:
        # p_NucleotideRecord
        print("\t|~> Returning promoter sequence")
        return_id = "p_" + str(nuc_rec['acc'])

        # Setting up the description for the FASTA record
        return_description = "Original_Query_Protein " + str(nuc_rec_in[0]) + " BLAST_Hit_Accession " + str(
            nuc_rec_in[1])

        # If there is only one coding region in the selected sequence, then
        # the sequence is returned unmodified.
        if len(coding_intervals) == 1:
            # Appends information to record description
            print("\t\t|~>No-additional-coding-regions-found-Returning-full-sequence")
            return SeqRecord(Seq.Seq(sequence), id=return_id, description=return_description)

        # If no coding intervals are indentified, None is returned.
        elif len(coding_intervals) == 0:
            print("\t\t|~> No coding intervals found for record: " + str(nuc_rec['acc']) + ".")
            return None

        # The start of the promoter is set to the start/end of the upstream gene
        # based on the directionality. ( --> --> or <-- -->)
        # If there was no downstream adjustment, then the last record in the list is upstream from the feature of interest.
        # If there was a downstream adjustment, then the second to last record in the list is upstream from the feature of interest.

        if stop_adj > 0:
            promoter_start = max(int(coding_intervals[-2][0]), int(coding_intervals[-2][1]))
        else:
            promoter_start = max(int(coding_intervals[-1][0]), int(coding_intervals[-1][1]))

        # Everything upstream of the promoter start is clipped off the
        # sequence and the substring is returned.
        return_seq = str(sequence[promoter_start:])
        # Appends information to record description
        print("\t\t|~>Successfully-clipped-off-upstream-coding-regions")
        return SeqRecord(Seq.Seq(return_seq), id=return_id, description=return_description)

    # If promoters aren't being isolated
    else:

        print("\t|~> Getting FASTA record")
        # Fetch the requested nucleotide sequence and return without any
        # modification.

        for i in range(REQUEST_LIMIT):

            try:

                handle = Entrez.efetch(db="nuccore", id=nuc_rec['acc'], strand=s_strand,
                                       seq_start=s_start, seq_stop=s_stop,
                                       rettype='fasta', retmode="txt")

                time.sleep(SLEEP_TIME)
                break

            except:

                print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")

                if i == (REQUEST_LIMIT - 1):
                    print("\t\tCould not download record after " + str(REQUEST_LIMIT) + " attempts")

        print("\t|~> Returnig sequence")
        return SeqIO.read(handle, "fasta")