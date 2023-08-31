from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import argparse

# Crear el analizador de argumentos
parser = argparse.ArgumentParser(description='Procesar los resultados BLAST y generar un archivo FASTA.')
parser.add_argument('-x', '--xml_file', help='Archivo de entrada BLAST en formato XML')
parser.add_argument('-f', '--fasta_file', help='Archivo de salida FASTA')
parser.add_argument('-n', '--num_results', type=int, default=10, help='Número de resultados BLAST a guardar en el archivo FASTA')

# Obtener los argumentos del comando
args = parser.parse_args()

# Ruta del archivo de entrada BLAST
input_file = args.xml_file

# Ruta del archivo de salida FASTA
output_file = args.fasta_file

# Número de resultados BLAST a guardar en el archivo FASTA
num_results = args.num_results

# Parsear los resultados BLAST del archivo XML
blast_records = NCBIXML.parse(open(input_file))

# Crear una lista para almacenar las secuencias y descripciones
sequences = []

# Leer el archivo FASTA existente y agregar la primera secuencia a la lista
existing_sequences = list(SeqIO.parse(output_file, "fasta"))
if existing_sequences:
    first_sequence = existing_sequences[0]
    sequences.append(first_sequence)

# Extraer información de los registros BLAST
for blast_record in blast_records:
    for alignment in blast_record.alignments[:num_results]:
        sequence_id = alignment.title
        sequence = alignment.hsps[0].sbjct
        description = alignment.hit_def
        
        print("Sequence ID:", sequence_id)
        print("Sequence:", sequence)
        print("Description:", description)
        
        seq_record = SeqRecord(Seq(sequence), id=sequence_id, description=description)
        sequences.append(seq_record)

# Guardar todas las secuencias en el archivo FASTA
SeqIO.write(sequences, output_file, "fasta")

