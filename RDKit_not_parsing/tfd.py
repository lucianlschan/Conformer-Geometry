import sys, os
from rdkit import Chem
from rdkit.Chem import TorsionFingerprints as TFP

dir_path = os.getcwd()

# load molecule
# confab.sdf is the lowest energy conformation found by confab
ei = Chem.SDMolSupplier(os.path.join(dir_path, 'EI_bayes_2345.sdf'))
lcb = Chem.SDMolSupplier(os.path.join(dir_path, 'LCB_bayes_2345.sdf'))
confab = Chem.SDMolSupplier(os.path.join(dir_path, 'confab.sdf'))
uniform = Chem.SDMolSupplier(os.path.join(dir_path, 'uniform_2345.sdf'))

# compute TFD

tfd_ei_confab = TFP.GetTFDBetweenMolecules(ei[0], confab[0])
tfd_lcb_confab = TFP.GetTFDBetweenMolecules(lcb[0], confab[0])
tfd_uniform_confab = TFP.GetTFDBetweenMolecules(uniform[0], confab[0])
