import numpy as np
import sys, os
import pandas as pd
import pybel
from pybel import ob
import openbabel


# Find reference conformation for each molecule 

ff = pybel._forcefields['mmff94']
def energy():
    return ff.Energy(False)


def analysis(file_in, seed):
    energy_uniform = []
    energy_200_uniform = []
    energy_EI = []
    energy_200_EI = []
    energy_LCB = []
    energy_200_LCB = []
    energy_confab = []
    path_record = [] 
    check_EI = []
    check_LCB= []	
    for paths, dirs, files in os.walk(file_in):
        if paths == file_in:
            continue
	path_record.append(os.path.basename(paths))
	print(paths)
        
        for x in files:
            if x.endswith('uniform_{}.sdf'.format(seed)):
                uniform_ob = pybel.readfile('sdf', os.path.join(paths, x)).next()
                ff.Setup(uniform_ob.OBMol)
                energy_uniform.append(energy())
            
            if x.endswith('EI_bayes_{}.sdf'.format(seed)):
                EI_ob = pybel.readfile('sdf', os.path.join(paths, x)).next()
                ff.Setup(EI_ob.OBMol)
                energy_EI.append(energy())

            if x.endswith('LCB_bayes_{}.sdf'.format(seed)):
                LCB_ob = pybel.readfile('sdf', os.path.join(paths, x)).next()
                ff.Setup(LCB_ob.OBMol)
                energy_LCB.append(energy())

            if x.endswith('uniform_200_{}.sdf'.format(seed)):
                uniform_ob = pybel.readfile('sdf', os.path.join(paths, x)).next()
                ff.Setup(uniform_ob.OBMol)
                energy_200_uniform.append(energy())

            if x.endswith('EI_bayes_200_{}.sdf'.format(seed)):
                EI_ob = pybel.readfile('sdf', os.path.join(paths, x)).next()
                ff.Setup(EI_ob.OBMol)
                energy_200_EI.append(energy())

            if x.endswith('LCB_bayes_200_{}.sdf'.format(seed)):
                LCB_ob = pybel.readfile('sdf', os.path.join(paths, x)).next()
                ff.Setup(LCB_ob.OBMol)
                energy_200_LCB.append(energy())

            if x.endswith('confab.sdf'):
                confabs_ob = pybel.readfile('sdf', os.path.join(paths, x)).next()
                ff.Setup(confabs_ob.OBMol)
                energy_confab.append(energy())

            if x.endswith("EI_evaluation_{}.txt".format(seed)):
                EI_eval = pd.read_table(os.path.join(paths, x), delimiter=" ")
                check_EI.append(len(EI_eval))

            if x.endswith("lcb_evaluation_{}.txt".format(seed)):
                LCB_eval = pd.read_table(os.path.join(paths, x), delimiter=" ")
                check_LCB.append(len(LCB_eval))

        
        if any(energy_200_EI) == True:

            EI_difference = np.array(energy_200_EI) - np.array(energy_confab)
            LCB_difference = np.array(energy_200_LCB) - np.array(energy_confab)
            uniform_difference = np.array(energy_200_uniform) - np.array(energy_confab)

        else:

            EI_difference = np.array(energy_EI) - np.array(energy_confab)
            LCB_difference = np.array(energy_LCB) - np.array(energy_confab)
            uniform_difference = np.array(energy_uniform) - np.array(energy_confab)


    data = pd.DataFrame({"Target": path_record, "Uniform": uniform_difference, "EI": EI_difference,"LCB": LCB_difference, "EI_iter": check_EI, "LCB_iter": check_LCB}, columns = ["Target","Uniform", "EI", "LCB", "EI_iter", "LCB_iter"])

    return data

# This commented block is used for molecules with 5 rotatable bonds.  
# We used the lowest energy conformation as reference point (consider outcome from both 100 iterations and 200 iterations)
"""
file_100_input = "/path/to/data/"
seeds_100 = [123,234,345,456,321]
results_100 = map(lambda x: analysis(file_100_input, x), seeds_100)
df_100 = [results_100[0].iloc[:,0:4], results_100[1].iloc[:,0:4], results_100[2].iloc[:,0:4], results_100[3].iloc[:,0:4], results_100[4].iloc[:,0:4]]
df_100_final = reduce(lambda left, right: pd.merge(left, right, on='Target'), df_100)
rename_seed = [["uniform_{}".format(x), "EI_bayes_{}".format(x), "LCB_bayes_{}".format(x)] for x in seeds_100]
renames = [item for sublist in rename_seed for item in sublist]
df_100_final.columns = ['Target'] + renames


file_input = '/path/to/data/' # for 200 iterations
seeds = [123, 234, 345, 456]
result = map(lambda x: analysis(file_input, x), seeds)
df = [result[0].iloc[:,0:4], result[1].iloc[:,0:4], result[2].iloc[:,0:4], result[3].iloc[:,0:4]]
df_final = reduce(lambda left, right: pd.merge(left, right, on='Target'), df)
rename_seed = [["uniform_200_{}".format(x), "EI_bayes_200_{}".format(x), "LCB_bayes_200_{}".format(x)] for x in seeds]
renames = [item for sublist in rename_seed for item in sublist]
df_final.columns = ['Target'] + renames

dataframe_combine = pd.merge(df_100_final, df_final, on='Target')
indxs = dataframe_combine.iloc[:,1:].idxmin(axis=1)
negative_molecules = []
minimum = dataframe_combine.iloc[:,1:].min(axis=1)
for item in minimum:
    if item <= 0:
        negative_molecules.append("Yes")
    else:
        negative_molecules.append("No")

dataframe = pd.DataFrame({"Target": df_final.iloc[:,0], "Lowest": indxs, "Best": negative_molecules, "N_rot": 5}, columns = ["Target", 'Lowest', 'Best', "N_rot"])

dataframe.to_csv(os.path.dirname(os.getcwd()) + "/result/summary/five_lowest_point_include_200_iterations.txt", index=False)

"""
# Check the lowest energy conformation 
file_input = '/path/to/data/'
seeds = [123, 234, 345, 456, 321]
result = map(lambda x: analysis(file_input, x), seeds)
df = [result[0].iloc[:,0:4], result[1].iloc[:,0:4], result[2].iloc[:,0:4], result[3].iloc[:,0:4], result[4].iloc[:,0:4]]
df_final = reduce(lambda left, right: pd.merge(left, right, on='Target'), df)
rename_seed = [["uniform_{}".format(x), "EI_bayes_{}".format(x), "LCB_bayes_{}".format(x)] for x in seeds]
renames = [item for sublist in rename_seed for item in sublist]
df_final.columns = ['Target'] + renames
indxs = df_final.iloc[:,1:].idxmin(axis=1)
negative_molecules = []
minimum = df_final.iloc[:,1:].min(axis=1)
for item in minimum:
    if item <= 0:
        negative_molecules.append("Yes")
    else:
        negative_molecules.append("No")

dataframe = pd.DataFrame({"Target": df_final.iloc[:,0], "Lowest": indxs, "Best": negative_molecules, "N_rot" : 1}, columns = ["Target", 'Lowest', 'Best', "N_rot"])



# Check for number of iterations in file (help remove early stopping molecules later)

df_check_EI_LCB = [result[0].iloc[:,[0,4,5]], result[1].iloc[:,[0,4,5]], result[2].iloc[:,[0,4,5]], result[3].iloc[:,[0,4,5]], result[4].iloc[:,[0,4,5]]]
df_check_final = reduce(lambda left, right: pd.merge(left, right, on='Target'), df_check_EI_LCB)
rename_check = [["EI_iter_{}".format(x), "LCB_iter_{}".format(x)] for x in seeds]
rename_check_1 = [item for sublist in rename_check for item in sublist]
df_check_final.columns = ["Target"] + rename_check_1
set_indx = df_check_final.iloc[:,1:].apply(set, axis=1)
set_len = [len(item) for item in set_indx]
early_stops = []
for item in set_len:
    if item == 1:
        early_stops.append("No")
    else:
        early_stops.append("Yes")
dataframe_check = pd.DataFrame({"Target": df_check_final.iloc[:,0], "Early_Stop": early_stops})


# Output file
dataframe.to_csv(os.path.dirname(os.getcwd()) + "/result/summary/one_lowest_point.txt", index=False)
dataframe_check.to_csv(os.path.dirname(os.getcwd()) + "/result/summary/one_iter_summary.txt", index=False)

