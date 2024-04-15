# Proteic Menger Curvature (PMC), Local Curvature (LC) and Local Flexibility (LF)

These 3 new metrics for disordered protein analysis were introduced in the following paper : 

"nP-collabs: Investigating counterion mediated bridges in the multiply phosphorylated tau-R2 repeat"
by Jules Marien [1,2], Chantal Prévost [1,2], and Sophie Sacquin-Mora [1,2]

[1] CNRS, Université de Paris, UPR 9080, Laboratoire de Biochimie Théorique, 13 rue Pierre et Marie Curie, 75005 Paris, France 

[2] Institut de Biologie Physico-Chimique-Fondation Edmond de Rotschild, PSL Research University, 13 rue Pierre et Marie Curie, 75005 Paris, France

# Description 

PMC is the application of the Menger curvature to the proteic backbone. Calculations are performed for each Cα($n$) where $n \in [s+1,N-s]$ with $s$ the spacing and $N$ the number of residues in the proteic chain from 1 to $N$. We advise to use $s$ = 2. PMCs present the great advantage to be an intrinsic property of each sampled conformation.

LC is the average of the PMC over every conformations for each residue.

LF is the standard deviation of the PMC over every conformations for each residue. 

LF provides an information on chain mobility for disordered proteic systems. It is philosophically in line with the RMSF but with the neat advantage of not requiring an alignement nor a reference structure. 


# Contents of interest

- A python module MODULE_Proteic_Menger_Curvature.py which provides the functions PMC, LC and LF to calculate the 3 metrics. 

- A python script Example_script_PMC_LC_LF.py which provides a workflow example for the calculation and plotting of PMCs, LCs and LFs

- Examples of plots to represent these metrics in the folder Plots.


# Package requirements
- MDAnalysis
- Numpy
- matplotlib (for plotting only)

# Minimal example for the calculation of PMCs

```
import numpy as np
import MDAnalysis as mda
import MODULE_Proteic_Menger_Curvature
    
result = MODULE_Proteic_Menger_Curvature.PMC('coordinates_file.pdb', 'trajectory.dcd', 'protein and name CA', spacing=2)
print(result)
np.savetxt("result.txt", result)
```

