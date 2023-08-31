registro="$1"
error="$2"
# Función para escribir en el archivo de registro
function escribir_registro {
    echo "$(date +"%Y-%m-%d %H:%M:%S") $1" >> "$registro"
}
echo

escribir_registro "Verificación de existencia de GenBank"

archivo="files/hemophiliaB.gb"
extension="${archivo##*.}"
archivo_out="files/protein_sequences.fasta"
script="scripts/ej1.py"
if [ -f "$archivo" ]; then
    escribir_registro "El archivo ($archivo) existe en el directorio 'files'."
    
else
    escribir_registro "Error: El archivo ($archivo) no existe en el directorio 'files'."
    echo "Error 2101: El archivo ($archivo) no existe en el directorio 'files'."
    echo "Error 2101">> "$error"
    exit 2
fi
if [ "$extension" = "gb" ]; then
    escribir_registro "La extensión del archivo ($archivo) es correcta."
else
    escribir_registro "El archivo ($archivo) no tiene la extensión esperada. La extención debe ser '.gb'"
    
    echo "Error 3100: La extención del archivo debe ser '.gb'"
    echo "Error 3100">> "$error"
    exit 2
fi


# Verificar la existencia del archivo "protein_sequences.fasta"
if [ ! -f "$archivo_out" ]; then
    touch "$archivo_out"
    escribir_registro "El archivo de salida ($archivo_out) se creó con éxito."
fi

# Otorgar permisos de escritura al archivo
chmod +w "$archivo_out"

# Verificar la existencia del archivo 'ej1.py' en el directorio 'scripts'
escribir_registro "Verificación de existencia de ($script)"

touch files/fasta.log
chmod 777 files/fasta.log

if [ -e "$script" ]; then
    escribir_registro "El archivo ($script) existe"
    escribir_registro "Inicia de ejecución en python de ($script)"
    chmod +x "$script"
    python3 "$script"  "$archivo" -o "$archivo_out" >> files/fasta.log 
    exit_code=$?
    if [ $exit_code -eq 0 ]; then
        escribir_registro "La ejecución del archivo ($script) fue exitosa."
    else
        echo "Error 1301: La ejecución del archivo ($script) falló."
    	echo "Error 1301" >> "$error"
    	escribir_registro "La ejecución del archivo ($script) falló."
        exit 1
    fi
else
     echo "Error 2302: El archivo ($script) no existe en el directorio scripts."
     echo "Error 2302" >> "$error"
     escribir_registro "El archivo ($script) no existe en el directorio scripts."
     exit 2
fi

