def search_blast(input_records, max_hits=50, e_cutoff=10E-10, tax_limit=None, min_cover=None, db="nr"):

    print("|~> BLAST search:")



    # The final list of BLAST hits
    hits = []

    # Gets the accession numbers for all the hits in the BLAST search and
    # appends hits[] with every unique record. 
    for input_record in input_records:

        # Fetches the protein record based off the accession number
        print("\t|~> Getting protein record")

        for i in range(REQUEST_LIMIT):

            try:

                handle = Entrez.efetch("protein", id=input_record, rettype="fasta",
                                       retmode="text")

                time.sleep(SLEEP_TIME)
                break

            except:

                print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")

                if i == (REQUEST_LIMIT - 1):
                    print("\t\tCould not download record after " + str(REQUEST_LIMIT) + " attempts")

        # Fetches the protein sequence to be used in the BLAST search
        print("\t|~> Getting protein sequence")
        input_seq = (SeqIO.read(handle, "fasta")).seq

        print("\t|~> Performing BLAST search: " + str(input_record))
        # Performs the appropriate BLAST search based
        if tax_limit != None:

            # Holds the parameter for the taxonomic limitation set
            taxon = ""

            if len(tax_limit) == 1:

                taxon = "txid" + str(tax_limit[0]) + "[orgn]"

            else:

                # Goes through each of the taxa limits appends to the overall entrez_querey parameter
                for i in range(len(tax_limit) - 1):
                    taxon = taxon + "txid" + str(tax_limit[0]) + "[orgn]" + " AND "

                taxon = taxon + "txid" + str(tax_limit[-1]) + "[orgn]"

            for i in range(REQUEST_LIMIT):

                try:

                    result_handle = NCBIWWW.qblast("blastp", "nr", input_seq,
                                                   entrez_query=taxon, expect=e_cutoff,
                                                   hitlist_size=max_hits)

                    # Parses the resulting hits as a list
                    print("\t|~> Getting records")
                    blast_records = list(NCBIXML.parse(result_handle))

                    time.sleep(SLEEP_TIME)
                    break

                except:

                    print("\t\tNCBI exception raised on attempt " + str(i) + "\n\t\treattempting now...")

                    if i == (REQUEST_LIMIT - 1):
                        print("\t\tCould not download record after " + str(REQUEST_LIMIT) + " attempts")



        '''
        #Parses the resulting hits as a list
        print("\t|~> Getting records")
        blast_records = list(NCBIXML.parse(result_handle))'''

        # Adds each unique accession number to hits[]
        for record in blast_records[0].alignments:

            curr_hit_rec = record.hit_id.split('|')[-2]
            print("\t\t|~> Analyzing hit " + str(curr_hit_rec))

            for hit in record.hsps:

                # Checks if hit meets the minimum coverage if provided


                    cov = (hit.query_end - hit.query_start + 1) / (len(input_seq))

                    if (cov >= min_cover):

                        # Check if the hit is already in the return list
                        if len(hits) == 0:
                            print("\t\t|~> Adding first hit (Coverage = " + str(cov) + "): " + str(curr_hit_rec))
                            hits.append((input_record, curr_hit_rec))
                        elif (not (curr_hit_rec in list(zip(*hits))[1])):
                            print("\t\t|~> Adding hit (Coverage = " + str(cov) + "): " + str(curr_hit_rec))
                            hits.append((input_record, curr_hit_rec))

                    # Prints error if the minimum coverage is not met
                    else:
                        print("\t\t|~> Hit did not meet coverage (Coverage = " + str(cov) + ") requirement: " + str(
                            curr_hit_rec))


    print("\t|~> Returning " + str(len(hits)) + " unique hits")
    return hits