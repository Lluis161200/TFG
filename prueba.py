# import two Bio.Blast modules required
from Bio.Blast import NCBIWWW, NCBIXML

# import sequence and SeqRecord functionality
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# import time library to "pause" between consecutive BLAST attempts
from time import sleep


def blastp_search(query, cutoff, nhits, db, tax_id, request_limit=5):

    # form the txidXXX[orgn] syntax from the numeric tax_id
    taxon = 'txid' + str(tax_id) + '[orgn]'

    # try block: attempts to grab the sequence from query_sequence
    # if the code fails, it is because query_sequence is not a SeqRecord object
    try:
        query_sequence = query.seq
    # if not a SeqRecord object, we simply inform the user and return None
    except:
        print('Incorrect SeqRecord object provided. Exiting function.')
        return (None)

    print('Running BLASTP search for ', query.id, ' on database ', db)
    print('This may take a while...')

    # performs "request_limit" attempts at BLASTing our sequence
    for i in range(request_limit):
        try:
            # perform BLASTp search and parse results
            handleresults = NCBIWWW.qblast(program='blastp', database=db,
                                           sequence=query_sequence, entrez_query=taxon,
                                           expect=cutoff, hitlist_size=nhits)
            blast_record = NCBIXML.read(handleresults)
            print('\tNCBI blastp search completed')
            # wait for one second before returning (in case we call the function right away)
            sleep(1)
            # if everything went well, we "break" the for loop (we don't do more iterations)
            break
        # handles the exemption (just inform the user, and reattempt in this case)
        except:
            print('\tNCBI exception raised on attempt ' + str(i) + '\n\treattempting now...')
            # wait three seconds before re-attempting
            sleep(3)
            if i == (request_limit - 1):
                print('\tCould not download record after ' + str(request_limit) + ' attempts')

    return (blast_record)


#we start by creating a sequence record, which is what our BLAST function expects
myprot = SeqRecord(Seq('MLNKIAMLGTEKTAEAVGVDKSQISRWKRDWIPKFSMLLAVLEWGVVDDDMARLARQVAAILTNKKRPAATERSEQIQMEF'),\
                         id='ALA45761.1',\
                         description='repressor protein C1 [Escherichia phage Lambda]')

#we next define the parameters for the BLAST search that we are going to use on our BLAST function
threshold = 1E-5
nhits = 10
database = 'nr'
tax_id = 2
request_limit=3

#and finally we call the BLAST function and store the result
blast_result = blastp_search(myprot, threshold, nhits, database, tax_id, request_limit)

#we want to retrieve the ID of the query that NCBI BLAST assigned to our search
print('The BLAST ID for your query was: ', blast_result.query_id)

#we are curious about how many descriptions have been returned
print('You have ', len(blast_result.descriptions), ' descriptions')

#we want to figure out what the function of our query sequence might be
#to do so, we look at the first hit
print('\nOur first hit is annotated as: ', blast_result.descriptions[0].title)
print('Based on the hit e-value', blast_result.descriptions[0].e, \
      ' we are pretty confident that this is the function of our query sequence')

print('\nOur mysterious protein is most likely a ', blast_result.descriptions[0].title.split('|')[-1:][0])

#let's first check that the alignment for our first hit does indeed correspnds to the description
print('\nOur first alignment is against ', blast_result.alignments[0].title)

#we'll assess now how many hsp's the alignment has
print('The alignment has ', len(blast_result.alignments[9].hsps), 'HSPs')

#and now we can visualize the alignment's single HSP
print('\nHere is your alignment! (well, technically your HSP)')
print(blast_result.alignments[9].hsps[0])

#we can make further inquiries into the hsp
#for instance, we can ask where in the query did the alignment start and where did it end
print('\nThe alignment started at query position: ',blast_result.alignments[9].hsps[0].query_start)
print('The alignment ended at query position: ',blast_result.alignments[9].hsps[0].query_end)
