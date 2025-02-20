# Launching a pCALVADOS simulation

The file pCALVADOS_residues_and_phosphoresidues.csv contains the CALVADOS 2 parameters developed by Tesei et al. for canonical amino acids, and the parameters derived by Perdikari et al for phosphoresidues

The notebook pCALVADOS_notebook.ipynb allows to run the simulation for p-state 2 of the Tau mutant described in the publication. 

Start by downloading these two files or pulling this directory



## How to set up the conda environment ? 

You will need to create a conda environment as follows :

```
conda create -n pCALVADOS_IDRLab python=3.10

conda activate pCALVADOS_IDRLab

conda install numpy matplotlib mdtraj ipykernel ipywidgets openmm=7.7 -c conda-forge --yes

pip install wget localcider==0.1.18
```


## How to modify the notebook ?

Change the NAME, SEQUENCE, Temperature, Ionic_strength, charged_N_terminal_amine, charged_C_terminal_carboxyl and charged_histidine to set your simulation details



## Output

Default is a 1ns equilibration and a 1µs production run, with 20000 frames produced for the production. You can modify nsteps, stride and eqsteps if you wish to modify this.
