import os
import shutil

def Obtain(fdf_file):
    Species = {}
    with open(fdf_file, 'r') as f:
        in_block = False

        for line in f:
            line = line.strip()

            if line == "%block ChemicalSpeciesLabel":
                in_block = True
                continue

            if line == "%endblock ChemicalSpeciesLabel":
                break

            if in_block:
                number, Z, symbol = line.split()
                Species[symbol] = {"Number": int(number), "Z": int(Z)}
    return Species

def psml_find(fdf_file,dir):
    
    file_path = os.path.join(os.getcwd(), fdf_file)
    elements = []

    in_block = False

    with open(file_path,"r") as f:
        line = f.read().strip().split("\n")
        for info in line:
            if "%block ChemicalSpeciesLabel" in info:
                in_block = True
                continue
            
            if "%endblock ChemicalSpeciesLabel" in info:
                in_block = False
                break

            if in_block:
                elements.append(info.strip().split()[-1])

    dir_psml = os.path.expanduser("~/Siesta_Standart_psml/Scalar_Relativistic/")

    for element in elements:

        file_psml = os.path.join(dir_psml,f"{element}.psml")
        
        shutil.copy(file_psml, dir)
    return
