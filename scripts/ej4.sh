#/bin/bash
registro="$1"
error="$2"
# Función para escribir en el archivo de registro
function escribir_registro {
    echo "$(date +"%Y-%m-%d %H:%M:%S") $1" >> "$registro"
}
echo
# Directorio de salida
archivo_out1="files/orfs_emboss.fasta"
archivo_out2="files/orfs_ordenados_emboss.fasta"
archivo_out3="files/longest_orf.fasta"
archivo_out4="files/report_secstructure.txt"
archivo_out5="files/primstructure.txt"
archivo_out6="files/pepinfo.pdf"

# Verificar la existencia de los archivos de salida
if [ ! -f "$archivo_out1" ]; then
    touch "$archivo_out1"
    escribir_registro "El archivo de salida ($archivo_out1) se creó con éxito."
fi

if [ ! -f "$archivo_out2" ]; then
    touch "$archivo_out2"
    escribir_registro "El archivo de salida ($archivo_out2) se creó con éxito."
fi

if [ ! -f "$archivo_out3" ]; then
    touch "$archivo_out3"
    escribir_registro "El archivo de salida ($archivo_out3) se creó con éxito."
fi

if [ ! -f "$archivo_out4" ]; then
    touch "$archivo_out4"
    escribir_registro "El archivo de salida ($archivo_out4) se creó con éxito."
fi

if [ ! -f "$archivo_out5" ]; then
    touch "$archivo_out5"
    escribir_registro "El archivo de salida ($archivo_out5) se creó con éxito."
fi

touch files/emboss.log
chmod 777 files/emboss.log

escribir_registro "Comienza obtencion de ORFs"

getorf -sequence files/hemophiliaB.gb -minsize 150 -outseq "$archivo_out1" 2>> files/emboss.log

escribir_registro "Comienza orden de ORFs"

sizeseq -sequences "$archivo_out1" -outseq "$archivo_out2" -descending N 2>> files/emboss.log

escribir_registro "Seleccion de ORF mas largo"

longest_orf=$(cat "$archivo_out2" | grep ">" -n | awk -F ":" {'print $1'} | tail -n1) 2>> files/emboss.log

tail -n +$longest_orf "$archivo_out2" > "$archivo_out3"

escribir_registro "Analizamos estructuras primarias y secundarias"

garnier -sequence "$archivo_out3" -outfile "$archivo_out4" 2>> files/emboss.log
pepinfo -sequence "$archivo_out3" -graph pdf "$archivo_out5" 1>> files/emboss.log 2>> files/emboss.log



mv pepinfo.pdf files/pepinfo.pdf


