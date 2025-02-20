# Design of the pipeline to convert the caorse-grained pCALVADOS simulation into a minimized all-atom trajectory

You will need to download the PULCHRA software [1] and VMD [2]


0. Use PULCHRA and the script_reconstruction_sequence_WT.py to rebuilt the all atom trajectory
The output is top_AA.pdb and traj_AA.dcd 

But the phosphoresidues are not phosphorylated ! So we need to rebuild them
Create the folder Trajectory_all_atoms and copy the top_AA.pdb and traj_AA.dcd files, and the topology_files folder

1. Rename the phosphoresidues 
The patching won't work if the residues do not have standard names, so we have to do on a copy of top_AA.pdb which we call top_AA_renamed.pdb for convenience:
HIS -> HSE
SEP -> SER
TPO -> THR

sed -i 's/HIS/HSE/g' top_AA_renamed.pdb 
sed -i 's/SEP/SER/g' top_AA_renamed.pdb 
sed -i 's/TPO/THR/g' top_AA_renamed.pdb 

We've chosen HSE quite arbitrarily since we could have chosen the other neutral form as well



2. Use the patch_phosphoresidues_and_add_hydrogens_V3.tcl with VMD
The script uses the CHARMM36m forcefield parameters to patch the phosphoresidues [3]
One must adapt the script to match the phosphoresidues to the sequence in Tau_sequences.txt !
The residues are numbered starting from 0 for some reason, we'll just play along 
Make sure to send the PDBs and PSFs to a subdirectory 


3. Create the constraint files
Use the SCRIPT_minimization_constraint_file.py to create the restraint files for the minimizations

5. Don't forget to save one PSF file before performing the minimization ! Take the modified_frame_0.psf and rename it as you please to keep as the topology file

5. Perform the minimization with the SCRIPT_minimization_*.sh files 
You can parallelize the calculation by running them at the same time, since each should run on a single CPU core

6. Load all the dcd files with the PSF file on VM using the command :

vmd saved_topology.psf minimization_frame_{?,??,???,????,?????}.dcd

7. Create a single DCD trajectory file with File -> Save Coordinates

8. Delete the rest of the individual DCDs

Here you go ! 



References

[1] Rotkiewicz, P, Skolnick J. 2008. Fast procedure for reconstruction of full-atom protein models from reduced representations. Journal of computational chemistry. 29(9):1460-5. https://sites.gatech.edu/cssb/pulchra/  https://github.com/euplotes/pulchra

[2] Humphrey, W., Dalke, A. and Schulten, K., "VMD - Visual Molecular Dynamics", J. Molec. Graphics, 1996, vol. 14, pp. 33-38. https://www.ks.uiuc.edu/Research/vmd/allversions/cite.html

[3] Huang J, Rauscher S, Nawrocki G, Ran T, Feig M, de Groot BL, Grubmüller H, MacKerell AD Jr. CHARMM36m: an improved force field for folded and intrinsically disordered proteins. Nat Methods. 2017 Jan;14(1):71-73. doi: 10.1038/nmeth.4067
