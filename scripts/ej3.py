import argparse
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO


# Crear el analizador de argumentos
parser = argparse.ArgumentParser(description='Realizar alineamiento de secuencias y generar un archivo HTML con colores.')
parser.add_argument('input_file', help='Archivo de entrada en formato FASTA')
parser.add_argument('-o', '--output_file', default='aligned_sequences_with_colors.html', help='Archivo de salida en formato HTML')
parser.add_argument('-f', '--alignment_file', default='top_blast_results_aligned.fasta', help='Archivo de salida para el alineamiento en formato FASTA')
parser.add_argument('--muscle_path', help='Ruta completa del ejecutable Muscle')

# Obtener los argumentos del comando
args = parser.parse_args()

# Nombre del archivo de entrada
input_file = args.input_file

# Nombre del archivo de salida HTML con colores
output_file = args.output_file

# Nombre del archivo de salida FASTA para el alineamiento
alignment_file = args.alignment_file

# Ruta completa del ejecutable Muscle
muscle_path = args.muscle_path

# Alinear secuencias con MUSCLE
muscle_cline = MuscleCommandline(cmd=muscle_path , input=input_file, out=alignment_file, diags=True, maxiters=1)
muscle_cline()

# Leer el alineamiento desde el archivo FASTA
alignment = AlignIO.read(alignment_file, "fasta")

# Definir un diccionario para asignar colores a los aminoácidos según sus propiedades
amino_acid_colors = {
    'A': 'green',   # Alanina
    'R': 'blue',    # Arginina
    'N': 'cyan',    # Asparagina
    'D': 'red',     # Ácido aspártico
    'C': 'yellow',  # Cisteína
    'Q': 'cyan',    # Glutamina
    'E': 'red',     # Ácido glutámico
    'G': 'green',   # Glicina
    'H': 'magenta', # Histidina
    'I': 'green',   # Isoleucina
    'L': 'green',   # Leucina
    'K': 'blue',    # Lisina
    'M': 'green',   # Metionina
    'F': 'green',   # Fenilalanina
    'P': 'green',   # Prolina
    'S': 'green',   # Serina
    'T': 'green',   # Treonina
    'W': 'green',   # Triptófano
    'Y': 'green',   # Tirosina
    'V': 'green'    # Valina
}

# Abrir el archivo de salida HTML para escribir
with open(output_file, "w") as f:
    # Escribir el alineamiento con aminoácidos coloreados en el archivo
    f.write("<html>\n")
    f.write("<body>\n")
    for record in alignment:
        f.write('<pre style="margin: 0;">\n')  # Usar <pre> para preservar el formato y eliminar el margen
        f.write("<b>" + record.id + "</b><br>")
        for char in record.seq:
            # Asignar un color según el aminoácido
            if char in amino_acid_colors:
                color = amino_acid_colors[char]
                f.write('<span style="color: ' + color + ';">' + char + '</span>')
            else:
                f.write(char)
        f.write("\n</pre>\n")
    f.write("</body>\n")
    f.write("</html>\n")

print("Alineación con colores guardada en:", output_file)

