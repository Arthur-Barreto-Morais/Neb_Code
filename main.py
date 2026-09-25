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
from ase.mep import NEB, DyNEB
from ase.optimize import BFGS

#---------------------------------------------------------------------------#

###############
##-Main Code-##
###############

parser = arg.ArgumentParser()
parser.add_argument("-f", dest = "label", type = str, required = True)
parser.add_argument("-n", dest = "num_images", type = int, required = True)
parser.add_argument("-np", dest = "n_proc", type = int, required = True)
parser.add_argument("-ci", dest = "CI_NEB", required= False, action = "store_true")
args = parser.parse_args()

label_file = args.label
num_images = args.num_images -2
n_proc = args.n_proc
want_CI = args.CI_NEB

fdf_file = label_file + ".fdf"

os.getcwd()

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

neb_file.interpolate(mic = True, interpolate_cell = False, method = 'idpp')

#Obtaining the Chemical Species
Species = chemical.chemical_Obtain(fdf_file)

dir_neb = 'neb'

if os.path.exists(dir_neb):
    shutil.rmtree(dir_neb)
os.makedirs(dir_neb, exist_ok = True)

for i, image in enumerate(Images):
    dir_image = os.path.join(dir_neb,f'image_{i}')
    
    os.makedirs(dir_image, exist_ok= True)
    
    shutil.copy(fdf_file,dir_image)
    
    chemical.psml_find(fdf_file,dir_image)

command_run =["/usr/bin/mpirun", "-np", str(n_proc), "/home/usr/bin/siesta-5.0.0/MPICH2/bin/siesta"]

#Neb calculation
for i, image in enumerate(Images):
    dir_image = os.path.join(dir_neb,f'image_{i}')
    
    image.calc = my_calc.Siesta(image_dir = dir_image, label_name = label_file, species = Species, command = command_run)

for file in ['neb_final.traj', 'neb.log', 'neb_relax.traj']:
    if os.path.exists(file):
        os.remove(file)


opt = BFGS(neb_file,logfile = 'neb.log',trajectory = "neb_relax.traj")

opt.run(fmax = 1.0, steps = 300)

ase_write('neb_final.traj',Images)

#Dy-Neb calculation
Dy_Images = [image.copy() for image in Images]

for i, image in enumerate(Dy_Images):
    dir_image = os.path.join(dir_neb,f'image_{i}')

    image.calc = my_calc.Siesta(image_dir = dir_image, label_name = label_file, species = Species, command = command_run)

Dy_neb_file = DyNEB(Dy_Images, k = 0.10, climb = False, dynamic_relaxation = True, method = 'improvedtangent', remove_rotation_and_translation= True)

for file in ['Dy_neb_final.traj', 'Dy_neb.log', 'Dy_neb_relax.traj']:
     if os.path.exists(file):
         os.remove(file)

opt = BFGS(Dy_neb_file,logfile = 'Dy_neb.log',trajectory = "Dy_neb_relax.traj")

opt.run(fmax = 0.1, steps = 300)

ase_write('Dy_neb_final.traj',Dy_Images)

if want_CI:
    #CI-Neb calculation
    CI_Images = [image.copy() for image in Dy_Images]

    for i, image in enumerate(CI_Images):
        dir_image = os.path.join(dir_neb,f'image_{i}')
        
        image.calc = my_calc.Siesta(image_dir = dir_image, label_name = label_file, species = Species, command = command_run)

    CI_neb_file = NEB(CI_Images, k = 0.10, climb = True, method = 'improvedtangent', remove_rotation_and_translation= True)

    for file in ['CI_neb_final.traj', 'CI_neb.log', 'CI_neb_relax.traj']:
        if os.path.exists(file):
            os.remove(file)

    CI_opt = BFGS(CI_neb_file,logfile = 'CI_neb.log',trajectory = "CI_neb_relax.traj")

    CI_opt.run(fmax = 0.05, steps = 300)

    ase_write('CI_neb_final.traj', CI_Images)
