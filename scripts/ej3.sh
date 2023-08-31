#!/bin/bash
registro="$1"
error="$2"
# Función para escribir en el archivo de registro
function escribir_registro {
    echo "$(date +"%Y-%m-%d %H:%M:%S") $1" >> "$registro"
}
echo

input_file="files/top_sequences.fasta"
html_output="files/aligned_sequences_with_colors.html"
fasta_output="files/aligned_sequences.fasta"
script_file="scripts/ej3.py"
script_file2="scripts/ej3_1.py"
muscle_executable="scripts/muscle"

escribir_registro "Verificar la existencia y no vaciedad del archivo ($input_file)"
if [ -s "$input_file" ]; then
    escribir_registro "El archivo XML ($input_file) existe y no está vacío."
else
   escribir_registro "Error: El archivo XML ($input_file) no existe o está vacío."
   echo "Error 4300" >> "$error"
   echo "Error 4300: El archivo FASTA ($input_file) no existe o está vacío."
    exit 1
fi

escribir_registro "Verificar la extensión del archivo de entrada"
if [[ ! "$input_file" =~ \.fasta$ ]]; then
    echo "Error 3300: El archivo de entrada ($input_file) no tiene la extensión .fasta."
    escribir_registro "El archivo de entrada ($input_file) no tiene la extensión .fasta."
    echo "Error 3300" >> "$error"
    exit 2
else
	escribir_registro "La extensión del archivo ($input_file) es correcta."
fi

escribir_registro "Verificar la existencia y permisos del archivo de salida ($html_output)"
if [ ! -f "$html_output" ]; then
    touch "$html_output"
    escribir_registro "El archivo de salida ($html_output) se creó con éxito."
fi
chmod +w "$html_output"

escribir_registro "Verificar la existencia y permisos del archivo de salida ($fasta_output)"
if [ ! -f "$fasta_output" ]; then
    touch "$fasta_output"
    escribir_registro "El archivo de salida ($fasta_output) se creó con éxito."
fi
escribir_registro "Verificar la existencia y permisos del ejecutable ($muscle_executable)"
if [ -f "$muscle_executable" ]; then
    # Obtener la ruta completa del ejecutable "muscle"
    muscle_path=$(realpath "$muscle_executable")
    chmod +x "$muscle_executable"
else
    echo "Error 2303: El ejecutable ($muscle_executable) no existe en el directorio 'scripts'."
    escribir_registro "Error: El ejecutable ($muscle_executable) no existe en el directorio 'scripts'."
    echo "Error 2303" >> "$error"
    exit 1
fi


chmod +w "$fasta_output"
touch files/MSA.log
chmod 777 files/MSA.log

escribir_registro "Verificar la existencia del archivo de ($script_file)"

if [ -e "$script_file" ]; then
    escribir_registro "El archivo ($script_file) existe"
    escribir_registro "Inicia de ejecución en python de ($script_file)"
    
    python3 "$script_file" "$input_file" -o "$html_output" -f "$fasta_output" --muscle_path "$muscle_path" >> files/MSA.log
    exit_code=$?
    if [ $exit_code -eq 0 ]; then
        escribir_registro "La ejecución del archivo ($script_file) fue exitosa."
    else
        echo "Error 1301: La ejecución del archivo ($script_file) falló."
    	echo "Error 1301" >> "$error"
    	escribir_registro "La ejecución del archivo ($script_file) falló."
        exit 1
    fi
else
     echo "Error 2302: El archivo ($script_file) no existe en el directorio scripts."
     echo "Error 2302" >> "$error"
     escribir_registro "El archivo ($script_file) no existe en el directorio scripts."
     exit 2
fi

escribir_registro "Verificar la existencia y no vaciedad del archivo ($fasta_output)"
if [ -s "$fasta_output" ]; then
    escribir_registro "El archivo  ($fasta_output) existe y no está vacío."
else
   escribir_registro "Error: El archivo  ($fasta_output) no existe o está vacío."
   echo "Error 4300" >> "$error"
   echo "Error 4300: El archivo ($fasta_output) no existe o está vacío."
    exit 1
fi

escribir_registro "Verificar la existencia del archivo de ($script_file2)"

if [ -e "$script_file2" ]; then
    escribir_registro "El archivo ($script_file2) existe"
    escribir_registro "Inicia de ejecución en python de ($script_file2)"
    
    python3 "$script_file2" "$fasta_output" >> files/MSA.log
    exit_code=$?
    if [ $exit_code -eq 0 ]; then
        escribir_registro "La ejecución del archivo ($script_file2) fue exitosa."
    else
        echo "Error 1301: La ejecución del archivo ($script_file2) falló."
    	echo "Error 1301" >> "$error"
    	escribir_registro "La ejecución del archivo ($script_file2) falló."
        exit 1
    fi
else
     echo "Error 2302: El archivo ($script_file2) no existe en el directorio scripts."
     echo "Error 2302" >> "$error"
     escribir_registro "El archivo ($script_file2) no existe en el directorio scripts."
     exit 2
fi
