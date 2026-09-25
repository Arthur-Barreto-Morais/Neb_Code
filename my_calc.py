import os
import subprocess
import numpy as np
from ase.calculators.calculator import Calculator, all_changes

class Siesta(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, image_dir, label_name, species, command, **kwargs):
        Calculator.__init__(self,**kwargs)

        self.image_dir = image_dir
        self.label_name = label_name
        self.species = species
        self.command = command

        return
    
    def calculate(self, atoms = None, properties = ['energy'], system_changes = all_changes):
        Calculator.calculate(self,atoms, properties, system_changes)

        self.create(atoms)

        self.account_run()

        Energy, Forces = self.get_variables()

        self.results['energy'] = Energy
        self.results['forces'] = Forces

        return
    
    #Auxiliary Functions        
    def create(self, atoms):
        fdf_aux = os.path.join(self.image_dir,self.label_name + ".fdf")

        new_pos = ''
        
        for atom in atoms:
            new_pos += ((
                f"{float(atom.position[0]):>12.9f} "
                f"{float(atom.position[1]):>12.9f} "
                f"{float(atom.position[2]):>12.9f} "
                f"{self.species[atom.symbol]['Number']:3d} "
                f"{int(atom.index + 1):3d}   "
                f"{atom.symbol:<2}\n"
                       ))
        
        with open(fdf_aux, 'r') as f:
            lines = f.read()

        begin = lines.find('%block AtomicCoordinatesAndAtomicSpecies')
        end = lines.find('%endblock AtomicCoordinatesAndAtomicSpecies')

        lines = (
            lines[:begin] +
            '%block AtomicCoordinatesAndAtomicSpecies\n' +
            new_pos +
            '%endblock AtomicCoordinatesAndAtomicSpecies\n' +
            lines[end + len('%endblock AtomicCoordinatesAndAtomicSpecies'):]
                )
        with open(fdf_aux, 'w') as f:
            f.write(lines)
        
        return
    
    def account_run(self):
        os.environ["OMP_NUM_THREADS"] = "1"

        fdf_aux = os.path.join(self.image_dir,self.label_name + ".fdf")
        out_aux = os.path.join(self.image_dir,self.label_name + ".out")
        
        with open(out_aux, 'w') as f_out, open(fdf_aux, 'r') as f_in:
            subprocess.run(
                self.command,
                cwd =self.image_dir,
                stdin =f_in,
                stdout =f_out,
                stderr =subprocess.STDOUT,
                check =True
                          )
        
        return
    
    def get_variables(self):
        out_aux = os.path.join(self.image_dir,self.label_name + ".out")

        Energy = None

        with open(out_aux, 'r') as f:
            for line in f:
                if 'Total =' in line:
                    Energy = float(line.split()[3])
        FA_aux = os.path.join(self.image_dir,self.label_name + ".FA")

        Forces = np.loadtxt(FA_aux, skiprows = 1)[:, 1:4]

        return Energy, Forces
