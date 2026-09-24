#############################
##-Code Idea-##
#############################
#1. Ajustar as configurações iniciais:
#   1.1. Arquivos de entrada: initial.XSF e final.XSF
#   1.2. Label name do arquivo fdf
#   1.3. Número de imagens intermediárias N_Images + 2.
#   1.4. Nome do diretório onde ocorre as contas.
#2. Criar o primeiro caminho inicial (first_path.traj).
#3. Construir os arquivos fdf para cada imagem.
#4. Rodar o SIESTA para cada imagem.
#5. Ler os arquivos de saída do SIESTA (Total Energy e Forces).
#6. Rodar o NEB com os dados obtidos.
#7. 
#---------------------------------------------------------------------------#

##################
##-Imports Used-##
##################

import os 
import shutil
from ase.io import read, write as ase_read, ase_write
from ase.calculators.emt import EMT
from ase.mep import NEB

#---------------------------------------------------------------------------#

###############
##-Main Code-##
###############

dir_run = os.getcwd()

file_name = 'File'

fdf_file = 'File' + ".fdf"

num_images = 'Value'

img_initial = 'initial.XSF'
xsf_initial = ase_read(img_initial, format='xsf')
xsf_initial.pbc = [True, True, False]
img_final = 'final.XSF'
xsf_final = ase_read(img_final, format='xsf')
xsf_final.pbc = [True, True, False]

Images = []

Images.append(xsf_initial)

for i in range(num_images):
    image = xsf_initial.copy()
    Images.append(image)

Images.append(xsf_final)

neb_file = NEB(Images, k = 0.10, climb = False, method = 'improvedtangent', remove_rotation_and_translation= True)

neb_file.interpolate(Images, mic = True, interpolate_cell = False, method = 'idpp')

#Obtaining the Chemical Species

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

dir_neb = 'neb'

if os.path.exists(dir_neb):
    shutil.rmtree(dir_neb)
os.makedirs(dir_neb, exist_ok = True)

for i, image in enumerate(Images):
    dir_neb_image = os.path.join(dir_neb,f'image_{i}')
    os.makedirs(dir_neb_image, exist_ok= True)
    
