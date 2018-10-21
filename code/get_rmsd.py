import numpy as np
import sys, os
import pandas as pd
import pybel
from pybel import ob
import openbabel
from scipy.stats import wilcoxon
import matplotlib.pyplot as plt
from itertools import chain

ff = pybel._forcefields['mmff94']
def energy():
    return ff.Energy(False)


def get_rmsd(source, file_in, seed):
    
    confab_rmsd_uniform = []
    confab_rmsd_EI = []
    confab_rmsd_LCB = []    
    
    bo_rmsd_confab = []
    bo_rmsd_uniform = []
    bo_rmsd_EI = []
    bo_rmsd_LCB = []

    bo_check_EI = []
    bo_check_LCB = []
    confab_check_EI = []
    confab_check_LCB = []

    bo_target = []
    confab_target = []

    for i in range(len(file_in)):
        if file_in.iloc[i,2] == "Yes":
            basenames = file_in.iloc[i,0] + '/' +file_in.iloc[i,1] + '.sdf'
            inputs = os.path.join(source,basenames)
            ref_mol = pybel.readfile("sdf", inputs).next()
            bo_target.append(file_in.iloc[i,0])

            # read EI_bayes
            EI_bases = file_in.iloc[i,0] + "/EI_bayes_{}.sdf".format(seed) 
            EI_input = os.path.join(source, EI_bases)
            EI_mol = pybel.readfile("sdf", EI_input).next()

            # read LCB_bayes
            LCB_bases = file_in.iloc[i,0] + "/LCB_bayes_{}.sdf".format(seed) 
            LCB_input = os.path.join(source, LCB_bases)
            LCB_mol = pybel.readfile("sdf", LCB_input).next()

            # read uniform
            uniform_bases = file_in.iloc[i,0] + "/uniform_{}.sdf".format(seed) 
            uniform_input = os.path.join(source, uniform_bases)
            uniform_mol = pybel.readfile("sdf", uniform_input).next()

            # read confab
            confab_bases = file_in.iloc[i,0] + "/confab.sdf" 
            confab_input = os.path.join(source, confab_bases)
            confab_mol = pybel.readfile("sdf", confab_input).next()

            # compute RMSD
            align = pybel.ob.OBAlign()
            align.SetRefMol(ref_mol.OBMol)
            
            # RMSD EI
            align.SetTargetMol(EI_mol.OBMol)
            align.Align()
            bo_rmsd_EI.append(align.GetRMSD())
            
            # RMSD LCB
            align.SetTargetMol(LCB_mol.OBMol)
            align.Align()
            bo_rmsd_LCB.append(align.GetRMSD())
            
            # RMSD Confab
            align.SetTargetMol(confab_mol.OBMol)
            align.Align()
            bo_rmsd_confab.append(align.GetRMSD())

            # RMSD Uniform
            align.SetTargetMol(uniform_mol.OBMol)
            align.Align()
            bo_rmsd_uniform.append(align.GetRMSD())
            
	
        else:
            basenames = file_in.iloc[i,0] + '/confab.sdf'
            inputs = os.path.join(source, basenames)
            ref_mol = pybel.readfile("sdf", inputs).next()
            confab_target.append(file_in.iloc[i,0])

            # read EI_bayes
            EI_bases = file_in.iloc[i,0] + "/EI_bayes_{}.sdf".format(seed) 
            EI_input = os.path.join(source, EI_bases)
            EI_mol = pybel.readfile("sdf", EI_input).next()

            # read LCB_bayes
            LCB_bases = file_in.iloc[i,0] + "/LCB_bayes_{}.sdf".format(seed) 
            LCB_input = os.path.join(source, LCB_bases)
            LCB_mol = pybel.readfile("sdf", LCB_input).next()

            # read uniform
            uniform_bases = file_in.iloc[i,0] + "/uniform_{}.sdf".format(seed) 
            uniform_input = os.path.join(source, uniform_bases)
            uniform_mol = pybel.readfile("sdf", uniform_input).next()

            # compute RMSD
            align = pybel.ob.OBAlign()
            align.SetRefMol(ref_mol.OBMol)
            
            # RMSD EI
            align.SetTargetMol(EI_mol.OBMol)
            align.Align()
            confab_rmsd_EI.append(align.GetRMSD())
            
            # RMSD LCB
            align.SetTargetMol(LCB_mol.OBMol)
            align.Align()
            confab_rmsd_LCB.append(align.GetRMSD())

            # RMSD Uniform
            align.SetTargetMol(uniform_mol.OBMol)
            align.Align()
            confab_rmsd_uniform.append(align.GetRMSD())

    bo_data = pd.DataFrame({"target": bo_target, "Uniform": bo_rmsd_uniform, "EI": bo_rmsd_EI, "LCB": bo_rmsd_LCB, "Confab": bo_rmsd_confab, "N_rot":1}, columns = ["target", "Uniform", "EI", "LCB", "Confab", "N_rot"])
    confab_data = pd.DataFrame({'target': confab_target, "Uniform": confab_rmsd_uniform, "EI": confab_rmsd_EI, "LCB": confab_rmsd_LCB, "N_rot":1}, columns = ["target","Uniform","EI","LCB", "N_rot"])   
    return  confab_data, bo_data


source = "/path/to/data/"
# {x}_lowest_point  x is the number of rotatable bond
file_in = pd.read_csv(os.path.dirname(os.getcwd()) + "/result/summary/one_lowest_point.txt") 

output = get_rmsd(source, file_in, 321)
output[0].to_csv(os.path.dirname(os.getcwd()) + '/result/rmsd/one_confab_rmsd_ref.txt', index=False)
output[1].to_csv(os.path.dirname(os.getcwd()) + '/result/rmsd/one_bo_rmsd_ref.txt', index=False)


