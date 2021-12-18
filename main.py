import regex as re
import math
import os

## Comprueba si la proteina introducida (en regex) es válida siguiendo las normas
## de la parte 1

def comprobar_proteina(proteina, archivo_string):

    proteina = validar_entera_regex(proteina, r">(?P<nombre>.+?)     (?P<num>\d+) nt( fragment)?\n(((\w{10} ){4}(\w{10})\n)*(\w{10} )*\w+)")

    if not(proteina):
        print("ERROR 3: La proteína no cumple con el formato FASTA.")
        print("  Línea de la proteína: " + numero_de_linea(proteina.splitlines()[0], archivo_string))
        return False

    nombre = proteina.group("nombre")

    # Cogemos los nucleótidos del match
    nucleotidos = proteina.group(4)

    # Quitamos los espacios y saltos de linea para mejor manejo de la cadena
    nucleotidos_stripped = reemplazar_regex(nucleotidos, "", r" |\n")

    # Traducimos las "T"s en "U"s
    nucleotidos_traducidos = reemplazar_regex(nucleotidos_stripped, "U", r"T")

    # Comprobamos que el número de nucleótidos es módulo de 3
    if not(validar_regex(nucleotidos_traducidos, r"^(.{3})+$")):
        print("ERROR 4: El número de nucleótidos de un gen no es múltiplo de 3.")
        print("  Línea de la proteína: " + numero_de_linea(nombre, archivo_string))
        return False    

    # Comprobamos que el número de nucleótidos coincide con el número indicado
    # Se podría hacer con regex el tema de len(string) = numero
    elif not(comprobar_nucleotidos(nucleotidos_traducidos, int(proteina.group("num")))):
        print("ERROR 5: El número de nucleótidos de un gen no coincide con el indicado en la primera línea.")
        print("  Línea de la proteína: " + numero_de_linea(nombre, archivo_string))
        return False
    
    # Comprobar que 'AUG' está presente al principio
    elif not(validar_regex(nucleotidos_traducidos, r"^AUG")):
        print("ERROR 6: La secuencia de ARN obtenido en el paso 1 no comienza con el codón de inicio 'AUG'.")
        print("  Línea de la proteína: " + numero_de_linea(nombre, archivo_string))
        return False

    # Comprobamos que termina con algún codón de final
    elif not(validar_regex(nucleotidos_traducidos, r"^.+?((UAA)|(UAG)|(UGA))$")):
        print("ERROR 7: La secuencia de ARN obtenido en el paso 1 no finaliza con alguno de los codones de fin.")
        print("  Línea de la proteína: " + numero_de_linea(nombre, archivo_string))
        return False

    # Comprobamos que el codón de final solo está presente al final
    elif validar_regex(nucleotidos_traducidos, r"^(.{3})*?((UAA)|(UAG)|(UGA)).+"):
        print("ERROR 8: La secuencia de ARN contiene un codón de parada en alguna posición distinta del final.")
        print("  Línea de la proteína: " + numero_de_linea(nombre, archivo_string))
        return False

    return True

# Busca una busqueda en una string y retorna la linea en la que se encuentra
def numero_de_linea(busqueda, string):
    for num, linea in enumerate(string.splitlines(), 1):
        if linea.__contains__(busqueda):
            return str(num)
        

def comprobar_nucleotidos(nucl, num):
    if len(nucl) == num:
        return True
    return False

# Funciones de Regex
def validar_regex(entrada, regex):
    regex_compilado = re.compile(regex)
    return regex_compilado.match(entrada)

def reemplazar_regex(entrada, sustitucion, regex):
    regex_compilado = re.compile(regex)
    return regex_compilado.sub(sustitucion, entrada)

def validar_entera_regex(entrada, regex):
    regex_compilado = re.compile(regex)
    return regex_compilado.fullmatch(entrada)

def validar_iter_regex(entrada, regex):
    regex_compilado = re.compile(regex)
    return regex_compilado.finditer(entrada)

def traducir_nucleotidos(nucleotidos):
    if validar_entera_regex(nucleotidos, r"UU[UC]"):
        return "F"
    elif validar_entera_regex(nucleotidos, r"UU[AG]|CU."):
        return "L"
    elif validar_entera_regex(nucleotidos, r"UC.|AG[CU]"):
        return "S"
    elif validar_entera_regex(nucleotidos, r"UA[UC]"):
        return "Y"
    elif validar_entera_regex(nucleotidos, r"UA[AG]"):
        return ""
    elif validar_entera_regex(nucleotidos, r"UG[UC]"):
        return "C"
    elif validar_entera_regex(nucleotidos, r"UGA"):
        return ""
    elif validar_entera_regex(nucleotidos, r"UGG"):
        return "W"
    elif validar_entera_regex(nucleotidos, r"CC."):
        return "P"
    elif validar_entera_regex(nucleotidos, r"CA[UC]"):
        return "H"
    elif validar_entera_regex(nucleotidos, r"CA[AG]"):
        return "Q"
    elif validar_entera_regex(nucleotidos, r"CG.|AG[AG]"):
        return "R"
    elif validar_entera_regex(nucleotidos, r"AU[ACU]"):
        return "I"
    elif validar_entera_regex(nucleotidos, r"AUG"):
        return "M"
    elif validar_entera_regex(nucleotidos, r"AC."):
        return "T"
    elif validar_entera_regex(nucleotidos, r"AA[CU]"):
        return "N"
    elif validar_entera_regex(nucleotidos, r"AA[GA]"):
        return "K"
    elif validar_entera_regex(nucleotidos, r"GU."):
        return "V"
    elif validar_entera_regex(nucleotidos, r"GC."):
        return "A"
    elif validar_entera_regex(nucleotidos, r"GA[CU]"):
        return "D"
    elif validar_entera_regex(nucleotidos, r"GA[GA]"):
        return "E"
    elif validar_entera_regex(nucleotidos, r"GG."):
        return "G"
    return "-"
    
def traducir_a_proteina(proteina):

    if proteina.group(3):
        fragment = proteina.group(3)
    else:
        fragment = ""

    cabezera = ">" + proteina.group("nombre") + "     " + str(math.floor(float(proteina.group("num")) / 3.0 - 1.0)) + " aa" + fragment

    nucleotidos = proteina.group(4)
    nucleotidos_stripped = reemplazar_regex(nucleotidos, "", r" |\n")
    nucleotidos_traducidos = reemplazar_regex(nucleotidos_stripped, "U", r"T")

    nucleotidos_nuevos = ""

    for grupo in validar_iter_regex(nucleotidos_traducidos, r".{3}"):
        nuevo = traducir_nucleotidos(grupo[0])
        if nuevo == '-':            # Hay un caso en el que hay dos nucleótidos que no tienen traducción, en concreto "NNN" y "NN_".
                                    #     Para este caso, hemos decidido no contar la proteína para no sacarla como resultado.
            return ""
        nucleotidos_nuevos += nuevo

    nucleotidos_nuevos_traducidos_1 = ""

    for linea in validar_iter_regex(nucleotidos_nuevos, r".{1,50}"):
        nucleotidos_nuevos_traducidos_1 += linea[0] + "\n"

    nucleotidos_nuevos_traducidos_2 = ""

    for linea in validar_iter_regex(nucleotidos_nuevos_traducidos_1, r".{1,10}\n?"):
        if validar_regex(linea[0], r".{1,10}\n"):
            nucleotidos_nuevos_traducidos_2 += linea[0]
        else:
            nucleotidos_nuevos_traducidos_2 += linea[0] + " "

    # Por cada grupo de 3, traducimos a una y hacemos append a una salida

    return cabezera + "\n" + nucleotidos_nuevos_traducidos_2

if __name__ == "__main__":

    # Parte 1.2

    debug = True

    if not(debug):
        camino_fichero = input("Hola buenos dias. Te solicito el nombre del fichero: ")
        formato_archivo = validar_entera_regex(camino_fichero, r"(\.?(?:.+\/)*)(.+)(DNA\.txt)$")
        existe_archivo = os.path.exists(camino_fichero)
        while not(
            formato_archivo and                   # Que tenga el formato '/dir1/.../nombre.DNA.txt'
            os.path.exists(camino_fichero)        # Que ese 'archivo' exista
            ):

            if not(formato_archivo):
                print("ERROR 1: El fichero no finaliza con 'DNA.txt'.")

            elif not(existe_archivo):
                print("ERROR 2: El fichero no existe.")

            camino_fichero = input("Hola buenos dias de nuevo. A ver como te lo explico... que este no vale socio.")
            formato_archivo = validar_entera_regex(camino_fichero, r"(\.?(?:.+\/)*)(.+)(DNA\.txt)$")
    else:
        camino_fichero = "/home/lunaesomnia/Documents/GitHub/ALFPrac2/All_C_genes_DNA.txt"
        formato_archivo = validar_entera_regex(camino_fichero, r"(\.?(?:.+\/)*)(.+)(DNA\.txt)$")
        existe_archivo = os.path.exists(camino_fichero)

    archivo = open(camino_fichero, "r").read()


    directorios = formato_archivo.group(1)     # directorios
    nombre_archivo = formato_archivo.group(2)  # nombre del fichero (sin DNA.txt)
    extension = formato_archivo.group(3)       # extension (DNA.txt)

    cadenas = []

    # Por cada proteína que cumple el formato FASTA, hacemos los checks
    for prot in validar_iter_regex(archivo, r">(?P<nombre>.+?)     (?P<num>\d+) nt( fragment)?\n(((\w{10} ){4}(\w{10})\n)*(\w{10} )*\w+)"):
        if comprobar_proteina(prot[0], archivo):
            cadenas.append(prot)

    proteinas = []

    # Por cada cadena de adn, generar su proteina
    for cadena in cadenas:
        nueva_proteina = traducir_a_proteina(cadena)
        if nueva_proteina != "":
            proteinas.append(nueva_proteina)

    # Generar un archivo en el directorio llamado ".... Protein.txt" y guardar los resultados
    with open(directorios + nombre_archivo + "Protein.txt", "w") as f:
        for proteina in proteinas:
            f.write(proteina + "\n\n")

# Parte 2


    
    
    