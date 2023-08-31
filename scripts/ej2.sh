#!/bin/bash
registro="$1"
error="$2"
# Función para escribir en el archivo de registro
function escribir_registro {
    echo "$(date +"%Y-%m-%d %H:%M:%S") $1" >> "$registro"
}
echo


input_file="files/protein_sequences.fasta"
xml_output="files/blast_results.xml"
fasta_output="files/top_sequences.fasta"
script_file="scripts/ej2.py"
script_file2="scripts/ej2_2.py"

escribir_registro "Verificación de la existencia del archivo de entrada en el directorio files"
if [ ! -f "$input_file" ]; then
    echo "Error 2201: El archivo de entrada ($input_file) no existe en el directorio files."
    escribir_registro "El archivo de entrada ($input_file) no existe en el directorio files."
    echo "Error 2201" >> "$error"
    exit 2
fi

# Verificar la extensión del archivo de entrada
if [[ ! "$input_file" =~ \.fasta$ ]]; then
    echo "Error 3200: El archivo de entrada ($input_file) no tiene la extensión .fasta."
    escribir_registro "El archivo de entrada ($input_file) no tiene la extensión .fasta."
    echo "Error 3200" >> "$error"
    exit 2
else
	escribir_registro "La extensión del archivo ($input_file) es correcta."
fi


escribir_registro "Verificar la existencia y permisos del archivo de salida ($xml_output)"
if [ ! -f "$xml_output" ]; then
    touch "$xml_output"
    escribir_registro "El archivo de salida ($xml_output) se creó con éxito."
fi
chmod +w "$xml_output"

escribir_registro "Verificar la existencia y permisos del archivo de salida ($fasta_output)"
if [ ! -f "$fasta_output" ]; then
    touch "$fasta_output"
    escribir_registro "El archivo de salida ($fasta_output) se creó con éxito."
fi
chmod +w "$fasta_output"

touch files/blast.log
chmod 777 files/blast.log

escribir_registro "Verificar la existencia del archivo de ($script_file)"
if [ -e "$script_file" ]; then
    escribir_registro "El archivo ($script_file) existe"
    chmod +x "$script_file"
    escribir_registro "Inicia de ejecución en python de ($script_file)"
    python3 "$script_file" "$input_file" -x "$xml_output" -f "$fasta_output" >> files/blast.log
    exit_code=$?
    if [ $exit_code -eq 0 ]; then
        escribir_registro "La ejecución del archivo ($script_file) fue exitosa."
    else
        echo "Error 1201: La ejecución del archivo ($script_file) falló."
    	echo "Error 1201" >> "$error"
    	escribir_registro "La ejecución del archivo ($script_file) falló."
        exit 1
    fi
else
     echo "Error 2202: El archivo ($script_file) no existe en el directorio scripts."
     echo "Error 2202" >> "$error"
     escribir_registro "El archivo ($script_file) no existe en el directorio scripts."
     exit 2
fi



escribir_registro "Verificar la existencia y no vaciedad del archivo FASTA"
if [ -s "$fasta_output" ]; then
    escribir_registro "El archivo FASTA ($fasta_output) existe y no está vacío."
else
    echo "Error 4200: El archivo FASTA ($fasta_output) no existe o está vacío."
    echo "Error 4200" >> "$error"
    escribir_registro "Error: El archivo FASTA ($fasta_output) no existe o está vacío."
    exit 1
fi


escribir_registro "Verificar la existencia del archivo de ($script_file2)"
if [ -e "$script_file" ]; then
    escribir_registro "El archivo ($script_file2) existe"
    chmod +x "$script_file"
    escribir_registro "Inicia de ejecución en python de ($script_file2)"
    python3 "$script_file2" -x "$xml_output" -f "$fasta_output" >> files/blast.log
    exit_code=$?
    if [ $exit_code -eq 0 ]; then
        escribir_registro "La ejecución del archivo ($script_file2) fue exitosa."
    else
        echo "Error 1201: La ejecución del archivo ($script_file2) falló."
    	echo "Error 1201" >> "$error"
    	escribir_registro "La ejecución del archivo ($script_file2) falló."
        exit 1
    fi
else
     echo "Error 2202: El archivo ($script_file2) no existe en el directorio scripts."
     echo "Error 2202" >> "$error"
     escribir_registro "El archivo ($script_file2) no existe en el directorio scripts."
     exit 2
fi

