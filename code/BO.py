
import GPy
import GPyOpt
import numpy as np
import pandas as pd
import math
import sys
import os
import random
import pybel
from pybel import ob


# Load molecules sequentially 
source = '/path/to/data/'

# Change parameter file for different set of molecules
parameter_file = pd.read_table(os.dirname(os.getcwd()) + "/parameter_file/four_rotatable_bond.txt", delimiter=",", header=None)
	
ff = pybel._forcefields['mmff94']

def energy():
    return ff.Energy(False)

rl = ob.OBRotorList()

def calculate_energy(torsion_angles):
    output = []
    for j in range(torsion_angles.shape[0]):
        for i in range(torsion_angles.shape[1]):
            rotors[i].SetToAngle(currentCoords, torsion_angles[j,i])
            molecule.OBMol.SetCoordinates(currentCoords)
            ff.SetCoordinates(molecule.OBMol)
        output.append(ff.Energy())
    return np.array(output).reshape([torsion_angles.shape[0],1])

set_seed = [123, 234, 345, 456] 

for paths, dirs, files in os.walk(source):
    if paths == source:
        continue
    for x in files:
        para_loc = np.argwhere([parameter_file.iloc[i,0] == os.path.basename(paths) for i in range(len(parameter_file))]).item()
	if x.endswith("base.sdf"):
            # load molecules
            molecule = pybel.readfile('sdf', os.path.join(paths, x)).next()
            # Property of molecule
            rl.Setup(molecule.OBMol)
            ff.Setup(molecule.OBMol)
            number = rl.Size()
            print(number)
            currentCoords = molecule.OBMol.GetCoordinates()
            rotIterator = rl.BeginRotors()
            dimensions = rl.Size()
            period_parameter = parameter_file.iloc[para_loc,1:].tolist()
            for seeds in set_seed:
                np.random.seed(seeds)
                rotor = rl.BeginRotor(rotIterator)
                rotors = []
                angles = []
                while rotor is not None:
                    angle = rotor.CalcTorsion(currentCoords)
                    angles.append(angle)
                    rotors.append(rotor)
                    rotor = rl.NextRotor(rotIterator)

                # Bayesian Optimization
                # Set Budget
                if number<=3:
                    max_iter = 45
                elif np.logical_and(number>=4, number<8):
                    max_iter = 95
                else:
                    max_iter = 245

                design_domain = [{'name': 'torsion_{}'.format(x), 'type': 'continuous', 'domain': (0, 2*math.pi)} for x in range(number)]

                # ARD1=True : allow different period parameters
                periodic = GPy.kern.StdPeriodic(input_dim=dimensions, period = period_parameter, ARD1=True)
                rbfs = GPy.kern.RBF(input_dim=dimensions)
                locally_periodic = GPy.kern.Prod([periodic,rbfs])
                lperiodic = periodic*rbfs

                model = GPyOpt.models.GPModel(kernel=locally_periodic, exact_feval=True,optimize_restarts=10,verbose=False)
                objective = GPyOpt.core.task.SingleObjective(calculate_energy)
                space = GPyOpt.Design_space(design_domain)
                aquisition_optimizer = GPyOpt.optimization.AcquisitionOptimizer(space)
                initial_design = GPyOpt.experiment_design.initial_design('random', space, 5)

                acquisition_EI = GPyOpt.acquisitions.AcquisitionEI(model, space, optimizer=aquisition_optimizer)
                acquisition_lcb = GPyOpt.acquisitions.AcquisitionLCB(model, space, optimizer=aquisition_optimizer)

                evaluator_EI = GPyOpt.core.evaluators.Sequential(acquisition_EI)
                evaluator_LCB = GPyOpt.core.evaluators.Sequential(acquisition_lcb)

                EI_bo = GPyOpt.methods.ModularBayesianOptimization(model, space, objective, acquisition_EI, evaluator_EI, initial_design, normalize_Y=True)

                lcb_bo = GPyOpt.methods.ModularBayesianOptimization(model, space, objective, acquisition_lcb, evaluator_LCB, initial_design, normalize_Y=True)

                EI_bo.run_optimization(max_iter = max_iter, evaluations_file=os.path.join(paths,'EI_evaluation_{}.txt'.format(seeds)))

                lcb_bo.run_optimization(max_iter = max_iter, evaluations_file=os.path.join(paths, 'lcb_evaluation_{}.txt'.format(seeds)))


                for i in range(len(EI_bo.x_opt)):
                    angle = EI_bo.x_opt[i]
                    angles[i] = angle
                    rotors[i].SetToAngle(currentCoords, angle)
                    molecule.OBMol.SetCoordinates(currentCoords)

                write_EI_file = pybel.Outputfile('sdf', os.path.join(paths, 'EI_bayes_{}.sdf'.format(seeds)), overwrite=True)
                write_EI_file.write(molecule)
                write_EI_file.close()

		for i in range(len(lcb_bo.x_opt)):
		    angle = lcb_bo.x_opt[i]
                    angles[i] = angle
                    rotors[i].SetToAngle(currentCoords, angle)
                    molecule.OBMol.SetCoordinates(currentCoords)

                write_LCB_file = pybel.Outputfile('sdf', os.path.join(paths, "LCB_bayes_{}.sdf".format(seeds)), overwrite=True)
                write_LCB_file.write(molecule)
                write_LCB_file.close()

                # remove object after each loop
                del EI_bo, lcb_bo

