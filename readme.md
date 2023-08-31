# Análisis de secuencia de ADN con scripts en Bash y Python

Este repositorio contiene una serie de scripts en Bash y Python que permiten realizar el análisis de una secuencia de ADN a partir de un archivo GenBank llamado `hemofilia.gb`.


## Requisitos
1. Los scripts en Bash ej1.sh, ej2.sh, ej3.sh y ej4.sh deben encontrarse en el directorio *scripts*.
2. Los scripts en Python ej1.py, ej2.py, ej2_2.py, ej3.py y ej3_1.py también deben encontrarse en el directorio *scripts*.
3. El archivo GenBank a analizar debe estar ubicado en el directorio *files*.
4. En la sección de alineamientos se utiliza el ejecutable MUSCLE, este viene incluido con el repositorio pero con la versión para Linux-64. Si se quiere correr en windows se debe reemplazar por el ejecutable de MUSCLE correspondiente.
5. Debido a que la sección de BLAST se ejecuta de manera remota, se debe tener una conexión WiFi o la red.

## Librerías
Si al intentar ejecutar el proyecto te encuentras con errores relacionados a la falta de librerías o dependencias, puedes instalarlas fácilmente utilizando el archivo requirements.txt. 
En la consola de comando posicionarse con el directorio en la carpeta descargada y utilice la siguiente línea de código:
 ```bash
    pip install -r requirements.txt
 ```
## Instrucciones de ejecución
1. Otorgar permisos de ejecución al script principal:
   ```bash
   chmod 777 main.sh
2. Ejecutar script principal:
   ```bash
   ./main.sh

## Resumen de scripts
### main.sh
Este script es el punto de entrada principal que ejecuta los demás scripts y registra los resultados en archivos de registro.

1. Crea el archivo "Errores.err" y lel archivo “registro.log”
3. Verifica la existencia de los directorios "files" y "scripts".
4. Ejecuta el archivo "ej1.sh" ubicado en el directorio "scripts" y registra la salida en el archivo "registro.log".
5. Ejecuta el archivo "ej2.sh" ubicado en el directorio "scripts" y registra la salida en el archivo "registro.log".
6. Ejecuta el archivo "ej3.sh" ubicado en el directorio "scripts" y registra la salida en el archivo "registro.log".
7. Ejecuta el archivo "ej4.sh" ubicado en el directorio "scripts".
8. Muestra errores y mensajes informativos en caso de que los archivos o directorios no existan o si la ejecución de los scripts falla.

### ej1.sh
Este script verifica la existencia del archivo GenBank, crea un archivo de salida protein_sequences.fasta si no existe y ejecuta el script ej1.py para traducir la secuencia de ADN a proteínas.

1. Verifica la existencia del archivo GenBank y su extensión correcta.
2. Crea el archivo de salida "protein_sequences.fasta" si no existe.
3. Verifica la existencia del archivo "ej1.py" en el directorio "scripts".
4. Ejecuta el archivo "ej1.py" y registro de los resultados en el archivo "fasta.log".
5. Muestra errores y mensajes informativos en caso de que los archivos o directorios no existan o si la ejecución de los scripts falla.


### ej1.py
Este script en Python lee una secuencia de ADN desde un archivo GenBank, traduce la secuencia a proteínas y guarda la secuencia de proteínas más larga en un archivo de salida en formato FASTA.

1. Importa las bibliotecas necesarias y define las funciones `complement_dna`, `get_seq` y `file_convert` para realizar operaciones en secuencias de ADN y proteínas.
2. Configura un parser de argumentos para recibir el archivo de entrada y el archivo de salida como argumentos de línea de comandos.
3. Parsea los argumentos de línea de comandos para obtener los nombres de archivo de entrada y salida.
4. Lee la secuencia de ADN desde el archivo GenBank de entrada.
5. Traduce la secuencia de ADN original y su cadena complementaria utilizando la función `get_seq`.
6. Imprime la secuencia de proteína más larga encontrada.
7. Convierte y guarda la secuencia de proteína más larga en el archivo de salida especificado en formato FASTA.


### ej2.sh
Este script verifica la existencia del archivo de entrada en formato FASTA, crea archivos de salida en formato XML y FASTA, y ejecuta los scripts ej2.py y ej2_2.py para realizar una búsqueda BLAST y guardar los resultados.

1. Verifica la existencia de un archivo de entrada FASTA en el directorio "files" y su correcta extensión. Si el archivo no existe, registra un error.
2. Verifica la existencia y los permisos de los archivos de salida XML y FASTA. Si el archivo no existe, se crea.
3. Crea un archivo de registro llamado "blast.log" para cargar las salidas de los scripts de python.
4. Verifica la existencia de un script de Python (ej2.py). Si el script existe, otorga permisos de ejecución y lo ejecuta utilizando Python con algunos argumentos de línea de comandos. La salida se redirige al archivo "blast.log". Si la ejecución del script falla, registra un error y finaliza. Si el archivo del script ej2.py no existe, registra un error y finaliza.
5. Verifica la existencia y el estado no vacío del archivo de salida FASTA. Si el archivo no existe o está vacío, registra un error y finaliza.
6. Verifica la existencia de otro script de Python (ej2_2.py). Si el script existe, otorga permisos de ejecución y lo ejecuta utilizando Python con algunos argumentos de línea de comandos. La salida se redirige al archivo "blast.log". Si la ejecución del script falla, registra un error y finaliza. Si el archivo del script no existe, registra un error y finaliza.
7. Muestra errores y mensajes informativos en caso de que los archivos o directorios no existan o si la ejecución de los scripts falla.

### ej2.py
Este script en Python realiza una búsqueda BLASTp utilizando la base de datos SwissProt y la matriz de sustitución BLOSSUM 62 (default); y guarda los resultados en un archivo XML y la secuencia query en un archivo FASTA.

1. Importa los módulos necesarios de la biblioteca Biopython: NCBIWWW, SeqIO y NCBIXML. Además, se importa “argparse” para analizar los argumentos de línea de comandos.
2. Crea un analizador de argumentos que acepta el archivo de entrada en formato FASTA, así como opciones para el archivo de salida en formato XML y FASTA. Los argumentos se obtienen del comando utilizando argparse.
3. Carga la secuencia de consulta desde el archivo FASTA especificado. Realiza una búsqueda BLAST remota utilizando NCBIWWW.qblast(), utilizando la base de datos Swissprot y el tipo de búsqueda "blastp".
4. Guarda los resultados en el archivo XML especificado. Procesa los resultados del archivo XML y se imprimen en pantalla los detalles de las alineaciones y los HSPs (High-Scoring Segment Pairs).
5. Guarda la secuencia de consulta en el archivo FASTA especificado.

### ej2_2.py
Este script en Python toma los resultados del archivo XML generado por ej2.py, selecciona las 10 mejores secuencias y las guarda en un archivo FASTA junto con la secuencia query utilizada en el BLAST.

1. Importa los módulos necesarios de la biblioteca Biopython: NCBIXML, Seq, SeqRecord y SeqIO. Además, se importa “argparse” para analizar los argumentos de línea de comandos.
2. Crea un analizador de argumentos que acepta el archivo XML de entrada, el archivo FASTA de salida y el número de resultados BLAST a guardar. Los argumentos se obtienen del comando utilizando argparse.
3. Lee el archivo FASTA existente y se agrega la primera secuencia a la lista.
4. Extrae información de los registros BLAST y se agrega las primeras 10 secuencias a la lista de secuencias.
5. Imprime la identificación de la secuencia, la secuencia en sí y su descripción.
6. Guarda todas las secuencias en el archivo FASTA especificado utilizando SeqIO.write(), reemplazando las secuencias existentes si el archivo ya contiene secuencias.

### ej3.sh
Este script verifica la existencia de archivos de entrada y salida, y ejecuta los scripts ej3.py y ej3_1.py para realizar el alineamiento de secuencias.

1. Define las rutas de los archivos de entrada y salida, así como las rutas de los scripts y el ejecutable necesarios.
2. Verifica la existencia y no vaciedad del archivo de entrada fasta y los archivos de salida fasta y html.
3. Verifica la existencia y los permisos del ejecutable "muscle".
4. Ejecuta los scripts ‘ej3.py’ y ‘ej3_1.py’  Python, pasando los archivos de entrada y salida, así como la ruta del ejecutable "muscle" según corresponda. Los resultados se redirigen a un archivo de registro “MSA.log”.
5. Se registran los resultados en el archivo de registro y se muestra un mensaje de éxito o error según corresponda.

### ej3.py
El código realiza el alineamiento de secuencias utilizando Muscle, asigna colores a los aminoácidos y genera dos archivos: un HTML con el alineamiento coloreado y un FASTA con los alineamientos sin colorear. Al final, se muestra un mensaje con la ubicación del archivo de salida HTML.

1. Importa las librerías necesarias, incluyendo argparse para el manejo de argumentos de línea de comandos, y módulos de BioPython para realizar el alineamiento de secuencias.
2. Define un analizador de argumentos para recibir el archivo de entrada, archivos de salida y la ruta del ejecutable Muscle.
3. Obtienen los argumentos ingresados por el usuario.
4. Utiliza el ejecutable Muscle para alinear las secuencias del archivo de entrada y se guarda el resultado en un archivo de salida en formato FASTA.
5. Lee el archivo de alineamiento y se asigna un color a cada aminoácido según sus propiedades. Luego, se escribe el alineamiento con los aminoácidos coloreados en un archivo HTML de salida.

### ej3_1.py
El código calcula los puntajes de conservación de las posiciones en un alineamiento, determina la identidad de secuencia entre pares de secuencias y construye un árbol filogenético basado en las distancias entre las secuencias. El resultado final es un archivo PNG que representa el árbol filogenético.

1. Importa las librerías necesarias, incluyendo `argparse` para el manejo de argumentos de línea de comandos, módulos de BioPython para trabajar con alineamientos y construcción de árboles filogenéticos, y `matplotlib` para la visualización.
2. Define un analizador de argumentos para recibir el archivo de entrada que contiene las secuencias alineadas en formato FASTA.
3. Carga las secuencias alineadas desde el archivo especificado.
4. Calcula los puntajes de conservación para cada posición en el alineamiento y se imprimen.
5. Calcula la identidad de secuencia entre cada par de secuencias y se imprimen los resultados. 
6. Construye un árbol filogenético basado en la matriz de distancias y se guarda como un archivo PNG.

### ej4.sh
Este script utiliza herramientas de Emboss para realiza el análisis de secuencias biológicas, incluyendo la búsqueda y selección del ORF más largo, el análisis de estructuras primarias y secundarias, y la generación de un gráfico en formato PDF con propiedades fisicoquímicas de proteínas. 

1. getorf: Es un programa que encuentra y extrae marcos abiertos de lectura (ORFs) en una secuencia de nucleótidos. Se utiliza con el archivo "hemophiliaB.gb" como entrada y se especifica un tamaño mínimo de ORF de 150 nucleótidos. El resultado se guarda en el archivo especificado por la variable "orfs_emboss.fasta".
2. sizeseq: Es un programa que analiza y modifica secuencias biológicas. En este caso, se utiliza para ordenar las secuencias según su longitud, de forma descendente (parámetro "descending N"). La entrada es el archivo generado anteriormente ("orfs_emboss.fasta") y el resultado se guarda en el archivo especificado por la variable "orfs_ordenados_emboss.fasta".
3. longest_orf: Esta línea utiliza una combinación de comandos como cat, grep, awk y tail para encontrar la línea que contiene el ORF más largo en el archivo "orfs_ordenados_emboss.fasta". Luego, utiliza tail para extraer las líneas desde la línea del ORF más largo hasta el final del archivo, y guarda el resultado en el archivo "longest_orf.fasta".
4. garnier: Es un programa que realiza un análisis de predicción de estructura secundaria de proteínas utilizando el método de Garnier-Osguthorpe-Robson (GOR). Se utiliza con el archivo "longest_orf.fasta" como entrada y se guarda el resultado en el archivo "$report_secstructure.txt".
5. pepinfo: Es un programa que realiza un análisis de propiedades fisicoquímicas de proteínas. Se utiliza con el archivo "longest_orf.fasta" como entrada y se genera un gráfico en formato PDF con los resultados, que se guardan también aparte en el archivo "primstructure.txt".
5. mv: Es un comando para mover o renombrar archivos. En este caso, se utiliza para mover el archivo "pepinfo.pdf" al directorio "files/".


## Errores

A continuación se enumeran los diferentes tipos de errores que puedes encontrar durante la ejecución:

1. Errores de ejecución de códigos (1000  - 1100 -1200 - 1300 -1400):
   - **1100**: Error de ejecución de código Bash de ejercicio 1.
   - **1101**: Error de ejecución de código Python de ejercicio 1.
   - **1200**: Error de ejecución de código Bash de ejercicio 2.
   - **1201**: Error de ejecución de código Python de ejercicio 2.
   - **1300**: Error de ejecución de código Bash de ejercicio 3.
   - **1301**: Error de ejecución de código Python de ejercicio 3.
   - **1400**: Error de ejecución de código Bash de ejercicio 4.
   - **1401**: Error de ejecución de código Python de ejercicio 4.
2. Errores de directorio de archivos (2000  - 2100 -2200 - 2300 -2400):
   - **2100**: Error de directorio de archivo Bash de ejercicio 1.
   - **2101**: Error de directorio de archivo de entrada de ejercicio 1.
   - **2102**: Error de directorio de archivo de Python de ejercicio 1.
   - **2200**: Error de directorio de archivo Bash de ejercicio 2.
   - **2201**: Error de directorio de archivo de entrada de ejercicio 2.
   - **2202**: Error de directorio de archivo de Python de ejercicio 2.
   - **2300**: Error de directorio de archivo Bash de ejercicio 3.
   - **2301**: Error de directorio de archivo de entrada de ejercicio 3.
   - **2302**: Error de directorio de archivo de Python de ejercicio 3.
   - **2400**: Error de directorio de archivo Bash de ejercicio 4.
   - **2401**: Error de directorio de archivo de entrada de ejercicio 4.
   - **2402**: Error de directorio de archivo de Python de ejercicio 4.
3. Errores de extensión (3100 - 3200 - 3300 - 3400):
   - **3100**: Error de extensión de archivo de ejercicio 1.
   - **3200**: Error de extensión de archivo de ejercicio 2.
   - **3300**: Error de extensión de archivo de ejercicio 3.
   - **3400**: Error de extensión de archivo de ejercicio 4.
4. Errores de archivos vacíos (4100 - 4200 - 4300 - 4400):
   - **4100**: Error de archivo vacío de ejercicio 1.
   - **4200**: Error de archivo vacío de ejercicio 2.
   - **4300**: Error de archivo vacío de ejercicio 3.
   - **4400**: Error de archivo vacío de ejercicio 4.

### Recomendaciones para manejo de errores



1. Errores de ejecución de códigos (1000  - 1100 -1200 - 1300 -1400):

Verifica que el código Bash que da error esté correctamente escrito y no contenga errores sintácticos, que se le esten dando los permisos de ejecución. Asegúrate de que se esté utilizando la versión correcta del intérprete de Bash. Asegúrate que no haya errores de otro estilo dentro del código de Bash.

Verifica que el código Python esté correctamente escrito y no contenga errores sintácticos. Asegúrate de que se esté utilizando la versión correcta del intérprete de Python. Asegúrate que las librerías se encuentren correctamente instaladas y los archivos de entrada y/o salida sean de la extensión correcta ynse estén pasando correctamente.

2. Errores de directorio de archivos (2000  - 2100 -2200 - 2300 -2400):

Asegúrate de que el archivo Bash se encuentre en el directorio especificado. Verifica la ruta del archivo y corrígela si es necesario.

Asegúrate de que el archivo de entrada se encuentre en el directorio especificado. Verifica la ruta del archivo y corrígela si es necesario.

Asegúrate de que el archivo Python del ejercicio 1 se encuentre en el directorio especificado. Verifica la ruta del archivo y corrígela si es necesario.

3. Errores de extensión (3100 - 3200 - 3300 - 3400):

Verifica que el archivo tenga la extensión correcta. Asegúrate de que el archivo tenga la extensión adecuada según el tipo de archivo esperado.

4. Errores de archivos vacíos (4100 - 4200 - 4300 - 4400):

Asegúrate de que el archivo contenga contenido. Verifica que el archivo no esté vacío y que contenga la información necesaria para la ejecución del código. Si el archivo es generado por un código de python anterior, asegúrate que esté ejecutándose correctamente.

