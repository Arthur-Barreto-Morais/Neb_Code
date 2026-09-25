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

import chemical
import my_calc
import os 
import shutil
import argparse as arg


from ase.io import read as ase_read, write as ase_write
from ase.calculators.emt import EMT
from ase.mep import NEB

#---------------------------------------------------------------------------#

###############
##-Main Code-##
###############

parser = arg.ArgumentParser()
parser.add_argument("-f", dest = "label", type = str, required = True)
parser.add_argument("-n", dest = "num_images", type = int, required = True)
parser.add_argument("-ci", dest = "CI_NEB", required= False, action = "store_true")
args = parser.parse_args()

label_file = args.label
num_images = args.num_images
want_CI = args.CI_NEB

fdf_file = label_file + ".fdf"


dir_run = os.getcwd()


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
chemical.chemical_Obtain(fdf_file)

dir_neb = 'neb'

if os.path.exists(dir_neb):
    shutil.rmtree(dir_neb)
os.makedirs(dir_neb, exist_ok = True)

for i, image in enumerate(Images):
    dir_neb_image = os.path.join(dir_neb,f'image_{i}')
    
    os.makedirs(dir_neb_image, exist_ok= True)
    
    shutil.copy(fdf_file,dir_neb_image)
    
    chemical.psml_find(fdf_file)

for i, image in enumerate(Images):
    dir_neb_image = os.path.join(dir_neb,f'image_{i}')
    image.calc = my_calc.Calc()
