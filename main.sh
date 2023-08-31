#!/bin/bash
touch Errores.err
chmod 777 Errores.err
error="Errores.err"

touch registro.log
chmod 777 registro.log
registro="registro.log"

# Función para escribir en el archivo de registro
function escribir_registro {
    echo "$(date +"%Y-%m-%d %H:%M:%S") $1" >> "$registro"
}
echo

#Archivos
script_file="scripts/ej1.sh"
script_file2="scripts/ej2.sh"
script_file3="scripts/ej3.sh"
script_file4="scripts/ej4.sh"


escribir_registro "Verificación de exitencia de directorios"
# Verificar la existencia del directorio "files"
if [ -d "files" ]; then
    escribir_registro "El directorio 'files' existe."
else
    escribir_registro "El directorio 'files' no existe. Creando directorio..."
    mkdir "files"
    escribir_registro "Directorio 'files' creado."
fi

# Verificar la existencia del directorio "scripts"
if [ -d "scripts" ]; then
    escribir_registro "El directorio 'scripts' existe."
else
    escribir_registro "El directorio 'scripts' no existe. Creando directorio..."
    mkdir "files"
    escribir_registro "Directorio 'scripts' creado."
fi

escribir_registro "Finalizó verificación de exitencia de directorios"


# Verificar la existencia del archivo ej1.sh
escribir_registro "Verificación de la existencia del archivo $script_file en el directorio scripts."
if [ -e "$script_file" ]; then
    escribir_registro "El archivo ($script_file) existe"
    chmod +x "$script_file"
    escribir_registro "Inicia de ejecución en python de ($script_file)"
    bash "$script_file" "$registro" "$error"
    exit_code=$?
    if [ $exit_code -eq 0 ]; then
        escribir_registro "La ejecución del archivo ($script_file) fue exitosa."
    else
        echo "Error 1100: La ejecución del archivo ($script_file) falló."
    	echo "Error 1100" >> "$error"
    	escribir_registro "La ejecución del archivo ($script_file) falló."
        exit 1
    fi
else
     echo "Error 2100: El archivo ($script_file) no existe en el directorio scripts."
     echo "Error 2100" >> "$error"
     escribir_registro "El archivo ($script_file) no existe en el directorio scripts."
     exit 2
fi


# Verificar la existencia del archivo ej2.sh
escribir_registro "Verificación de la existencia del archivo $script_file2 en el directorio scripts."
if [ -e "$script_file2" ]; then
    escribir_registro "El archivo ($script_file2) existe"
    chmod +x "$script_file2"
    escribir_registro "Inicia de ejecución en python de ($script_file2)"
    bash "$script_file2" "$registro" "$error"
    exit_code=$?
    if [ $exit_code -eq 0 ]; then
        escribir_registro "La ejecución del archivo ($script_file2) fue exitosa."
    else
        echo "Error 1200: La ejecución del archivo ($script_file2) falló."
    	echo "Error 1200" >> "$error"
    	escribir_registro "La ejecución del archivo ($script_file2) falló."
        exit 1
    fi
else
     echo "Error 2200: El archivo ($script_file2) no existe en el directorio scripts."
     echo "Error 2200" >> "$error"
     escribir_registro "El archivo ($script_file2) no existe en el directorio scripts."
     exit 2
fi


# Verificar la existencia del archivo ej3.sh
escribir_registro "Verificación de la existencia del archivo $script_file3 en el directorio scripts."
if [ -e "$script_file3" ]; then
    escribir_registro "El archivo ($script_file3) existe"
    chmod +x "$script_file3"
    escribir_registro "Inicia de ejecución en python de ($script_file3)"
    bash "$script_file3" "$registro" "$error"
    exit_code=$?
    if [ $exit_code -eq 0 ]; then
        escribir_registro "La ejecución del archivo ($script_file3) fue exitosa."
    else
        echo "Error 1300: La ejecución del archivo ($script_file3) falló."
    	echo "Error 1300" >> "$error"
    	escribir_registro "La ejecución del archivo ($script_file3) falló."
        exit 1
    fi
else
     echo "Error 2300: El archivo ($script_file3) no existe en el directorio scripts."
     echo "Error 2300" >> "$error"
     escribir_registro "El archivo ($script_file3) no existe en el directorio scripts."
     exit 2
fi

# Verificar la existencia del archivo ej4.sh
escribir_registro "Verificación de la existencia del archivo $script_file4 en el directorio scripts."
if [ -e "$script_file4" ]; then
    escribir_registro "El archivo ($script_file4) existe"
    chmod +x "$script_file4"
    escribir_registro "Inicia de ejecución en python de ($script_file4)"
    bash "$script_file4" "$registro" "$error"
    exit_code=$?
    if [ $exit_code -eq 0 ]; then
        escribir_registro "La ejecución del archivo ($script_file4) fue exitosa."
    else
        echo "Error 1400: La ejecución del archivo ($script_file4) falló."
    	echo "Error 1400" >> "$error"
    	escribir_registro "La ejecución del archivo ($script_file4) falló."
        exit 1
    fi
else
     echo "Error 2400: El archivo ($script_file4) no existe en el directorio scripts."
     echo "Error 2400" >> "$error"
     escribir_registro "El archivo ($script_file4) no existe en el directorio scripts."
     exit 2
fi
