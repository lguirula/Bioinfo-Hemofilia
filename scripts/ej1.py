import argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def complement_dna(sequence):
    complement = ""
    for base in sequence:
        if base == "A":
            complement += "T"
        elif base == "T":
            complement += "A"
        elif base == "G":
            complement += "C"
        elif base == "C":
            complement += "G"
        else:
            complement += base
    return complement

def get_seq(DNA_seq, strand):
    origin = Seq(DNA_seq)
    reading_frames = [0, 1, 2]
    frame_array = []
    array = []
    for frame in reading_frames:
        sequence = origin[frame:]
        if len(sequence) % 3 != 0:
            sequence += "N" * (3 - (len(sequence) % 3))
        protein = sequence.translate()
        sequences_array = protein.split("*")
        longest = max(sequences_array, key=len)
        array.append((longest, frame, strand))
        if len(sequence) % 3 != 0:
            sequence += "N" * (3 - (len(sequence) % 3))
    frame_array = max(array, key=lambda x: len(x[0]))
    return frame_array

def file_convert(list_seq, output_file):
    sequences = []
    for i, seq_info in enumerate(list_seq):
        protein_seq = seq_info[0]
        if protein_seq:
            seq_record = SeqRecord(Seq(protein_seq), id=f'Unknown protein sequence {i}', description=f'Longest protein sequence from frame {seq_info[1]} strand {seq_info[2]}')
            sequences.append(seq_record)
    if sequences:
        SeqIO.write(sequences, output_file, 'fasta')

# Configurar el parser de argumentos
parser = argparse.ArgumentParser(description='Descripción del script')
parser.add_argument('input_file', help='Archivo de entrada')
parser.add_argument('-o', '--output_file', default='protein_sequences.fasta', help='Archivo de salida')

# Parsear los argumentos de línea de comandos
args = parser.parse_args()

# Obtener los nombres de archivo de los argumentos
archivo_entrada = args.input_file
archivo_salida = args.output_file

# Obtener la secuencia del archivo GenBank
record = SeqIO.read(archivo_entrada, "genbank")
origin = record.seq

# Traducir la secuencia de origen
s = []
s.append(get_seq(origin, '+'))

# Obtener la cadena complementaria
complement = complement_dna(origin)
s.append(get_seq(complement, '-'))

print(f'Longest proteins: {s}')

ss = []
ss.append(max(s, key=lambda x: len(x[0])))

# Convertir y guardar en el archivo de salida
file_convert(ss, archivo_salida)

