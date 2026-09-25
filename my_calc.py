import os
import subprocess
import numpy as np
from ase.calculators.calculator import Calculator, all_changes

class Calc(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, image_dir, label_name, species, command, **kwargs):
        Calculator.__init__(**kwargs)

        self.image_dir = image_dir
        self.label_name = label_name
        self.species = species
        self.command = command
