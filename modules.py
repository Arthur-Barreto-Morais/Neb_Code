#Módulos
import numpy as np

def Distance(atom1,atom2):
    return np.sqrt(sum((atom1[i] - atom2[i])**2 for i in range(len(atom1))))

def Angle(atom1,atom2,atom3):
    pass
