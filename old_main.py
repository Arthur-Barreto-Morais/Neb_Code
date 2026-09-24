#Base do Corpo do código
#
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

import os
import shutil
import subprocess
import numpy as np
from ase.io import read, write
from ase.mep import NEB
from ase.io import write
from ase.visualize import view
from ase.calculators.singlepoint import SinglePointCalculator
from ase.calculators.calculator import Calculator, all_changes
from ase.optimize import BFGS
from ase.mep import NEB, NEBTools
import matplotlib.pyplot as plt

class SiestaCalculator(Calculator):

    implemented_properties = ['energy', 'forces']

    def __init__(self, image_dir, label_name, species, command, **kwargs):

        super().__init__(**kwargs)

        self.image_dir = image_dir
        self.label_name = label_name
        self.species = species
        self.command = command

    def calculate(
        self,
        atoms=None,
        properties=['energy'],
        system_changes=all_changes
    ):

        super().calculate(atoms, properties, system_changes)

        # ============================================================
        # 1. Caminho do FDF
        # ============================================================

        Dir_Aux = os.path.join(
            self.image_dir,
            self.label_name + '.fdf'
        )

        # ============================================================
        # 2. Criar as novas posições
        # ============================================================

        New_Positions = ''

        for atom in atoms:

            New_Positions += (
                f"{float(atom.position[0]):>12.9f} "
                f"{float(atom.position[1]):>12.9f} "
                f"{float(atom.position[2]):>12.9f} "
                f"{self.species[atom.symbol]['Number']:3d} "
                f"{int(atom.index + 1):3d}   "
                f"{atom.symbol:<2}\n"
            )

        # ============================================================
        # 3. Substituir AtomicCoordinates no FDF
        # ============================================================

        with open(Dir_Aux, 'r') as f:
            lines = f.read()

        begin = lines.find(
            '%block AtomicCoordinatesAndAtomicSpecies'
        )

        end = lines.find(
            '%endblock AtomicCoordinatesAndAtomicSpecies'
        )

        lines = (
            lines[:begin]
            + '%block AtomicCoordinatesAndAtomicSpecies\n'
            + New_Positions
            + '%endblock AtomicCoordinatesAndAtomicSpecies\n'
            + lines[
                end
                + len('%endblock AtomicCoordinatesAndAtomicSpecies'):
            ]
        )

        with open(Dir_Aux, 'w') as f:
            f.write(lines)

        # ============================================================
        # 4. Rodar SIESTA
        # ============================================================

        os.environ["OMP_NUM_THREADS"] = "1"

        Dir_Aux_log = os.path.join(
            self.image_dir,
            self.label_name + '.out'
        )

        with open(Dir_Aux_log, "w") as file_out, \
             open(Dir_Aux, "r") as file_in:
            subprocess.run(
                self.command,
                cwd=self.image_dir,
                stdin=file_in,
                stdout=file_out,
                stderr=subprocess.STDOUT,
                check=True
            )
        # ============================================================
        # 5. Ler energia
        # ============================================================

        Energy = None

        with open(Dir_Aux_log, 'r') as f:

            for line in f:

                if 'Total =' in line:

                    Energy = float(line.split()[3])

        if Energy is None:

            raise RuntimeError(
                f"Não foi possível encontrar a energia em {Dir_Aux_log}"
            )

        # ============================================================
        # 6. Ler forças
        # ============================================================

        FA_File = os.path.join(
            self.image_dir,
            self.label_name + '.FA'
        )

        FA_data = np.loadtxt(
            FA_File,
            skiprows=1
        )

        Forces = FA_data[:, 1:4]

        # ============================================================
        # 7. Entregar os resultados para o ASE
        # ============================================================

        self.results['energy'] = Energy
        self.results['forces'] = Forces

########################
#INITIAL CONFIGURATIONS#
########################

Image_Initial = 'initial.XSF'
Image_Final = 'final.XSF'

Directory = os.getcwd()

#Verify the fdf file exists in the directory:
#while True:
#    Label_Name = input("Enter the Label Name: ")
#
#    File_Name = Label_Name + '.fdf'
#    if os.path.isfile(os.path.join(Directory, File_Name)):
#        break
#    print(f"File {File_Name} not found in {Directory} . Please try again.")

Label_Name = 'Grafino'

Fdf_File = Label_Name + '.fdf'

#Fdf_File print test:
#print(Fdf_File)

N_Images = 5


initial = read(Image_Initial, format='xsf')
final = read(Image_Final, format='xsf')

initial.pbc = [True, True, False]
final.pbc = [True, True, False]

Images = [initial]

for i in range(0, N_Images):
    Images.append(initial.copy())
Images.append(final)

print(len(Images))	

Initial_Original = Images[0].positions.copy()
Final_Original = Images[-1].positions.copy()

neb=NEB(Images, method='improvedtangent', climb=False)

neb.interpolate()

#Visualization of initial path
#write('first_path.traj', Images)
#view(Images)

#Pega as espécies quimicas
with open(Fdf_File, 'r') as f:
    lines = f.readlines()

Species = {}
In_block = False
for line in lines:
    if "%block ChemicalSpeciesLabel" in line:
        In_block = True
        continue
    if "%endblock ChemicalSpeciesLabel" in line:   
        In_block = False
        break
    if In_block:
        Species_Number, Species_Z, Species_Symbol = line.split()
        Species[Species_Symbol] = {
                                   'Number': int(Species_Number),
                                   'Z': int(Species_Z)
                                  }
#Verify the Chemical Species Dictionary:
#print(Species)

#Create the directory for SIESTA calculations
Dir_Path = 'neb' 

if os.path.exists(Dir_Path):
    shutil.rmtree(Dir_Path)
os.makedirs(Dir_Path, exist_ok=True)

#Move the psml and fdf files to the directory and replace the positions for each image:
for i, image in enumerate(Images):

    Image_Dir = os.path.join(
        Dir_Path,
        f'image_{i}'
    )

    os.makedirs(Image_Dir, exist_ok=True)

    for Species_Symbol in Species.keys():

        shutil.copy(
            Species_Symbol + '.psml',
            Image_Dir
        )

    Dir_Aux = os.path.join(
        Image_Dir,
        Label_Name + '.fdf'
    )

    shutil.copy(
        Fdf_File,
        Dir_Aux
    )       
os.environ["OMP_NUM_THREADS"] = "1"

command = [
    "/usr/bin/mpirun",
    "-np",
    "10",
    "/home/usr/bin/siesta-5.0.0/MPICH2/bin/siesta"
]
for i, image in enumerate(Images):

    Image_Dir = os.path.join(
        Dir_Path,
        f'image_{i}'
    )

    image.calc = SiestaCalculator(
        image_dir=Image_Dir,
        label_name=Label_Name,
        species=Species,
        command=command
    )


neb_forces = neb.get_forces()

print(neb_forces)

for file in ['neb_final.traj', 'neb.log']:
    if os.path.exists(file):
        os.remove(file)
        print(f"Arquivo antigo removido: {file}")
        
opt = BFGS(neb,logfile='neb.log')

opt.run(
    fmax=0.1,
    steps=300
)

write('neb_final.traj', Images)

Images_NEB = [image.copy() for image in Images]

CI_Images = [image.copy() for image in Images]

for i, image in enumerate(CI_Images):

    Image_Dir = os.path.join(
        Dir_Path,
        f'image_{i}'
    )

    image.calc = SiestaCalculator(
        image_dir=Image_Dir,
        label_name=Label_Name,
        species=Species,
        command=command
    )


CI_neb=NEB(CI_Images, method='improvedtangent', climb=True)

for file in ['CI_neb_final.traj', 'CI_neb.log']:
    if os.path.exists(file):
        os.remove(file)
        print(f"Arquivo antigo removido: {file}")

CI_opt = BFGS(CI_neb,logfile='CI_neb.log')

CI_opt.run(
    fmax=0.05,
    steps=300
)

write('CI_neb_final.traj', CI_Images)
