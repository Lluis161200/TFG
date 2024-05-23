import time

from Bio import Entrez, SeqIO, Seq
from Bio.SeqRecord import SeqRecord

REQUEST_LIMIT = 5
SLEEP_TIME = .5
Entrez.email = "lluistfg@gmail.com"
def get_nuc_rec_from_prot(prot_id):


    print("|~> Getting nucelotide records for " + str(prot_id))

    # Fetch the protein record
    print("\t|~> Fetching the protein records")

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
                "U": 3, "AF": 3, "AY": 3, "DQ": 3, "KT": 3, "OM": 3,}

    return process_gene_records(p_scores, records, prot_id)



def process_gene_records(p_scores, records, prot_id):
    print("\t|~> Processing gene records")
    max_p_record = 0

    if 'ProteinList' in records['IPGReport'].keys():

        for idprotein in records['IPGReport']['ProteinList']:

            if 'CDSList' in idprotein.keys():
                print("\t|~> Fetching the CDS")

                for cds in idprotein['CDSList']:
                    cds_acc = cds.attributes['accver']
                    for key in p_scores:
                        if cds_acc.startswith(key):
                            if max_p_record < p_scores[key]:
                                max_p_record = p_scores[key]
                                # create the record
                                cds_rec = {'acc': cds_acc, 'start': cds.attributes['start'],
                                           'stop': cds.attributes['stop'], 'strand': cds.attributes['strand'],
                                           'p_score': p_scores[key]}

    if max_p_record != 0:
        return cds_rec
    else:
        print("\t|~> No gene records found for " + prot_id)
        return None

def get_nuc_seq(nuc_rec_in, start_adj=250, stop_adj=3):
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



