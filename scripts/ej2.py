from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML
import argparse

# Crear el analizador de argumentos
parser = argparse.ArgumentParser(description='Realizar una búsqueda BLAST utilizando una secuencia de consulta.')
parser.add_argument('input_file', help='Archivo de entrada en formato FASTA')
parser.add_argument('-x', '--xml_output', default='blast_results.xml', help='Archivo de salida para los resultados BLAST en formato XML')
parser.add_argument('-f', '--fasta_output', default='top_sequence.fasta', help='Archivo de salida para la secuencia de consulta en formato FASTA')

# Obtener los argumentos del comando
args = parser.parse_args()

# Ruta del archivo de entrada
input_file = args.input_file

# Cargar la secuencia de consulta desde el archivo FASTA
record_dict = SeqIO.index(input_file, "fasta")
lista_ids = list(record_dict.keys())
secuencia = record_dict[lista_ids[0]].format('fasta')

# Realizar la búsqueda BLAST remota con la base de datos Swissprot
result_handle = NCBIWWW.qblast("blastp", "swissprot", secuencia)

# Guardar los resultados en el archivo XML
with open(args.xml_output, "w") as out_handle:
    out_handle.write(result_handle.read())

result_handle.close()

# Procesar los resultados BLAST del archivo XML
result_handle = open(args.xml_output)
blast_record = NCBIXML.read(result_handle)

# Imprimir los resultados
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        print(f"Sequence: {alignment.title}")
        print(f"Length: {alignment.length}")
        print(f"Score: {hsp.score}")
        print(f"E-value: {hsp.expect}")
        print(f"Gaps: {hsp.gaps}")
        print(f"Query: {hsp.query}")
        print(f"Match: {hsp.match}")
        print(f"Subject: {hsp.sbjct}")
        print("\n")

# Guardar la secuencia de consulta en formato FASTA
with open(args.fasta_output, "w") as out_handle:
    SeqIO.write(record_dict[lista_ids[0]], out_handle, "fasta")
