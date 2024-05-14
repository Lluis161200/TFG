from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

Entrez.email = "lluistfg@gmail.com"
def obtener_secuencia_proteina_por_id(id_proteina):
    handle = Entrez.efetch(db="protein", id=id_proteina, rettype="fasta", retmode="text")
    fasta_data = handle.read()
    handle.close()

    # Parsear el formato FASTA para obtener el ID y la secuencia
    lines = fasta_data.split('\n')
    seq_id = lines[0][1:]
    seq_data = ''.join(lines[1:])

    # Crear un objeto SeqRecord
    seq_record = SeqRecord(Seq(seq_data), id=seq_id, description="")

    return seq_record
def search_blast(id_proteina, db, tax_id, cutoff, nhits, request_limit, min_coverage):

    hits=[]

    # Obtener la secuencia de la proteína utilizando su ID
    print("Obteniendo sequencia de la proteina:", id_proteina)

    query_sequence = obtener_secuencia_proteina_por_id(id_proteina)

    print("Iniciando busqueda en blast de la proteina:", id_proteina)

    #Creamos el id taxonomico con el formato correcto
    taxon = str(tax_id) + '[orgn]'


    # Realizar la búsqueda BLAST en NCBI
    for i in range(request_limit):
        try:

            result_handle = NCBIWWW.qblast(program='blastp', database=db,
                                           sequence=query_sequence, entrez_query=taxon,
                                           expect=cutoff, hitlist_size=nhits)

            print("Busqueda finalizada")
            blast_records = list(NCBIXML.parse(result_handle))
            break

        except:

            print("Erro en la busqueda de la proteina en blast")

    # Comprobación si cada hit cumple con el min_coverage
    for record in blast_records[0].alignments:

        curr_hit_rec = record.hit_id.split('|')[-2]
        print("Comprobación del coverage de : ", str(curr_hit_rec))

        for hit in record.hsps:
            cov = (hit.query_end - hit.query_start + 1) / (len(query_sequence))

            if cov >= min_coverage:

                #Guardamos solo los hits unicos (comprobamos si el id del hit ya existe en la lista final de hits si es asi no se guarda):
                second_components = [hit[1] for hit in hits]

                if len(hits) == 0:
                    print("Primer hit:", str(curr_hit_rec))
                    hits.append((id_proteina, curr_hit_rec))
                elif curr_hit_rec in second_components:
                    print("Hit ya existente")
                else:
                    print("Nuevo hit:" + str(curr_hit_rec))
                    hits.append((id_proteina, curr_hit_rec))
            else:
                print("El hit no cumple el min_coverage: ", min_coverage, ">", cov)
    return hits






