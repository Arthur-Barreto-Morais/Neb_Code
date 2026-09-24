import os
import shutil

fdf_file = "Grafino.fdf"

dir = os.getcwd()

file_path = os.path.join(dir,fdf_file)

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
