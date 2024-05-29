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
                                print()
                                cds_rec = {'acc': cds_acc, 'start': cds.attributes['start'],
                                           'stop': cds.attributes['stop'], 'strand': cds.attributes['strand'],
                                           'p_score': p_scores[key]}
                                print(cds_rec)

    if max_p_record != 0:
        return cds_rec
    else:
        print("\t|~> No gene records found for " + prot_id)
        return None


def get_nuc_seq(nuc_rec_in):
    nuc_rec = nuc_rec_in[2]

    if nuc_rec is None:
        print("\t|~>Record is not valid")
        return None

    print("|~> Get nucleotide sequence for " + str(nuc_rec['acc']))

    print("\t|~> Adjusting start and stop positions")

    # Initialize variables
    acc = nuc_rec['acc']
    strand = nuc_rec['strand']
    start = int(nuc_rec['start'])
    stop = int(nuc_rec['stop'])

    # Get the previous gene information
    try:
        handle = Entrez.elink(dbfrom="nuccore", id=acc, linkname="nuccore_nuccore")
        record = Entrez.read(handle)
        handle.close()

        if not record[0]["LinkSetDb"]:
            print("\t|~> No linked records found")
            return None

        prev_gene_id = record[0]["LinkSetDb"][0]["Link"][0]["Id"]

        handle_summary = Entrez.esummary(db="nuccore", id=prev_gene_id)
        summary = Entrez.read(handle_summary)
        handle_summary.close()

        prev_start = int(summary[0]['Start'])
        prev_stop = int(summary[0]['Stop'])

    except Exception as e:
        print(f"\t|~> Error fetching previous gene: {e}")
        return None

    # Adjust the start and stop positions based on the strand
    if strand == '+':
        s_start = prev_stop + 1
        s_stop = start - 1
        s_strand = 1
    else:
        s_start = stop + 1
        s_stop = prev_start - 1
        s_strand = 2

    print("\t|~> Getting FASTA record")

    for i in range(REQUEST_LIMIT):
        try:
            handle = Entrez.efetch(db="nuccore", id=acc, strand=s_strand,
                                   seq_start=s_start, seq_stop=s_stop,
                                   rettype='fasta', retmode="text")
            time.sleep(SLEEP_TIME)
            break
        except Exception as e:
            print(f"\t\tNCBI exception raised on attempt {i}\n\t\treattempting now... {e}")

            if i == (REQUEST_LIMIT - 1):
                print("\t\tCould not download record after " + str(REQUEST_LIMIT) + " attempts")
                return None

    print("\t|~> Returning sequence")
    return SeqIO.read(handle, "fasta")




