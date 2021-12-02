Esta es la práctica 2 de ALF (2º grado de UM)

## FASTA:

- \>{identificador_de_enzima}{5_espacios_blanco}{nº_nucleotido}{espacio}{"nt"}({espacio}{fragment})?

- El número de nucleótidos es siempre % 3

## PROGRAMA:

1. Solicitamos el nombre del fichero con el ADN de genes en FASTA descargado de la página [REBASE](http://rebase.neb.com/rebase/rebase.seqs.html) (comprobano que finaliza con la cadena `DNA.txt` (o precedido de un directorio)).

2. Generar un fichero en el mismo directorio con las cadenas de aminoacidos traduciendolo a proteinas (terminado en `Protein.txt`)

3. (Opcional) -> Comprobar la corrección de la traducción

4. Detección de errores:
    1. El nombre del fichero de entrada no tiene el sufijo DNA.txt. -> Solicitar nombre del fichero de entrada
    2. El fichero de entrada no existe.  -> Solicitar nombre del fichero de entrada
    3. El fichero de entrada no respeta el formato FASTA. 
    4. El número de nucleótidos de un gen no es múltiplo de 3. 
    5. El número de nucleótidos de un gen no coincide con el indicado en la primera línea. 
    6. La secuencia del ARN obtenido en el paso 1 no comienza con el codón de inicio AUG. 
    7. La secuencia del ARN obtenido en el paso 1 no finaliza con alguno de los codones de fin. 
    8. La secuencia del ARN contiene un codón de parada en alguna posición distinta del final. 

    En caso de que el error no sea ninguno de estos ^ (Pog) hay que indicar la línea y clase de error.