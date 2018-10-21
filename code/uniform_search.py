
# Pybel
from __future__ import print_function
import sys
import os
import math
import random
import pybel
from pybel import ob
import numpy as np 
from numpy.random import seed


# Seed for 5 rotor (200 iterations)   123,234,345,456
# Seed for 1-6 rotor (except four rotor)  123, 234, 345,456,321
# Seed for 4 rotor  1234,2345,3456,4567,321
seeds = [123,234,345,456]
ff = pybel._forcefields['mmff94']

def energy():
    return ff.Energy(False)

rl = ob.OBRotorList()

source = '/path/to/data/'
for dirpath, subdirs, files in os.walk(source):
    print(dirpath)
    #if os.path.exists(os.path.join(dirpath, "uniform.sdf")):continue
    for x in files:
        if x.endswith('base.sdf'):
            input_file = os.path.join(dirpath, x)
	    print(input_file)
            mol = pybel.readfile('sdf', input_file).next()
            rl.Setup(mol.OBMol)
            if rl.Size()<=3:
                K = 50
            elif np.logical_and(rl.Size()>=4, rl.Size()<=7):
                K = 200
            else:
                K = 250
                
            currentCoords = mol.OBMol.GetCoordinates()
            rotIterator = rl.BeginRotors()
            rotor = rl.BeginRotor(rotIterator)
            rotors = []
            angles = [] 
            # Get the rotors and angles
            while rotor is not None:
                angle = rotor.CalcTorsion(currentCoords)
                angles.append(angle)
                rotors.append(rotor)
                rotor = rl.NextRotor(rotIterator)

            for seed in seeds:
                min_energy = []
                min_angle = np.empty([K,rl.Size()])
                for j in range(K):
                    for i in range(len(angles)):
                        angle = random.uniform(0, 2.0 * math.pi)
                        angles[i] = angle
                        rotors[i].SetToAngle(currentCoords, angle)
                        mol.OBMol.SetCoordinates(currentCoords)
                    ff.SetCoordinates(mol.OBMol)
                    if j==0:
                        min_energy = energy()
                        min_angle[j,] = angles
                    elif energy()>=min_energy:
                        min_energy = min_energy
                        min_angle[j,] = min_angle[j-1,]
                    else:
                        min_energy = energy()
                        min_angle[j,] = angles


                for i in range(len(angles)):
                    angle = min_angle[K-1,i]
                    angles[i] = angle
                    rotors[i].SetToAngle(currentCoords, angle)
                    mol.OBMol.SetCoordinates(currentCoords)
            
                change_name = input_file.strip('base.sdf')
                new_name = change_name + 'uniform_{}.sdf'.format(seed)
                write_file = pybel.Outputfile('sdf', new_name, overwrite=True)
                write_file.write(mol)
	        write_file.close()
