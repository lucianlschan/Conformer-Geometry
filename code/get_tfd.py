import numpy as np
import sys, os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import TorsionFingerprints as TFP

def get_tfd(source, file_in, seed):
    
    confab_tfd_uniform = []
    confab_tfd_EI = []
    confab_tfd_LCB = []    
    
    bo_tfd_confab = []
    bo_tfd_uniform = []
    bo_tfd_EI = []
    bo_tfd_LCB = []

    bo_check_EI = []
    bo_check_LCB = []
    confab_check_EI = []
    confab_check_LCB = []

    bo_target = []
    confab_target = []

    for i in range(len(file_in)):
        if file_in.iloc[i,2] == "Yes":
            basenames = file_in.iloc[i,0] + '/' + file_in.iloc[i,1] + '.sdf'
            inputs = os.path.join(source,basenames)
            ref_mol = Chem.SDMolSupplier(inputs)
            bo_target.append(file_in.iloc[i,0])

            # read EI_bayes
            EI_bases = file_in.iloc[i,0] + "/EI_bayes_{}.sdf".format(seed) 
            EI_input = os.path.join(source, EI_bases)
            EI_mol = Chem.SDMolSupplier(EI_input)

            # read LCB_bayes
            LCB_bases = file_in.iloc[i,0] + "/LCB_bayes_{}.sdf".format(seed) 
            LCB_input = os.path.join(source, LCB_bases)
            LCB_mol = Chem.SDMolSupplier(LCB_input)

            # read uniform
            uniform_bases = file_in.iloc[i,0] + "/uniform_{}.sdf".format(seed) 
            uniform_input = os.path.join(source, uniform_bases)
            uniform_mol = Chem.SDMolSupplier(uniform_input)

            # read confab
            confab_bases = file_in.iloc[i,0] + "/confab.sdf" 
            confab_input = os.path.join(source, confab_bases)
            confab_mol = Chem.SDMolSupplier(confab_input)

            
            bo_tfd_EI.append(TFP.GetTFDBetweenMolecules(EI_mol[0], ref_mol[0]))
            bo_tfd_LCB.append(TFP.GetTFDBetweenMolecules(LCB_mol[0], ref_mol[0]))
            bo_tfd_confab.append(TFP.GetTFDBetweenMolecules(confab_mol[0], ref_mol[0]))
            bo_tfd_uniform.append(TFP.GetTFDBetweenMolecules(uniform_mol[0], ref_mol[0]))
            
	
        else:
            basenames = file_in.iloc[i,0] + '/confab.sdf'
            inputs = os.path.join(source, basenames)
            ref_mol = Chem.SDMolSupplier(inputs)
            confab_target.append(file_in.iloc[i,0])

            # read EI_bayes
            EI_bases = file_in.iloc[i,0] + "/EI_bayes_{}.sdf".format(seed) 
            EI_input = os.path.join(source, EI_bases)
            EI_mol = Chem.SDMolSupplier(EI_input)

            # read LCB_bayes
            LCB_bases = file_in.iloc[i,0] + "/LCB_bayes_{}.sdf".format(seed) 
            LCB_input = os.path.join(source, LCB_bases)
            LCB_mol = Chem.SDMolSupplier(LCB_input)

            # read uniform
            uniform_bases = file_in.iloc[i,0] + "/uniform_{}.sdf".format(seed) 
            uniform_input = os.path.join(source, uniform_bases)
            uniform_mol = Chem.SDMolSupplier(uniform_input)


            confab_tfd_EI.append(TFP.GetTFDBetweenMolecules(EI_mol[0], ref_mol[0]))
            confab_tfd_LCB.append(TFP.GetTFDBetweenMolecules(LCB_mol[0], ref_mol[0]))
            confab_tfd_uniform.append(TFP.GetTFDBetweenMolecules(uniform_mol[0], ref_mol[0]))

    bo_data = pd.DataFrame({"target": bo_target, "Uniform": bo_tfd_uniform, "EI": bo_tfd_EI, "LCB": bo_tfd_LCB, "Confab": bo_tfd_confab, "N_rot":6}, columns = ["target", "Uniform", "EI", "LCB", "Confab","N_rot"])
    confab_data = pd.DataFrame({'target': confab_target, "Uniform": confab_tfd_uniform, "EI": confab_tfd_EI, "LCB": confab_tfd_LCB, "N_rot":6}, columns = ["target","Uniform","EI","LCB", "N_rot"])   
    return  confab_data, bo_data


source = "/path/to/data/"
file_in = pd.read_csv(os.path.dirname(os.getcwd()) + "/result/summary/six_lowest_point.txt")
output = get_tfd(source, file_in, 321)
output[0].to_csv(os.path.dirname(os.getcwd()) + '/result/tfd/six_confab_tfd_ref.txt', index=False)
output[1].to_csv(os.path.dirname(os.getcwd()) + '/result/tfd/six_bo_tfd_ref.txt', index=False)


