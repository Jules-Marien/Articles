# Github repository for the publication on pCALVADOS

This repository contains the additional contents described in the article :

"Hydrodynamic Radius Determination of Tau and AT8 Phosphorylated Tau Mutants: A Combined Simulation and Experimental Study"


by Cynthia Lohberger [1][∗], Jules Marien [4][∗], Clarisse Bridot [2], Chantal Prévost [4], Diane Allegro [1], Mario Tatoni [3], Isabelle Landrieu [2], Caroline Smet-Nocca [2], Sophie Sacquin-Mora [4][#] and Pascale Barbier [1][#]

[1] Aix Marseille Univ, CNRS, Institut de Neurophysiopathologie, Marseille, France
[2] Université de Lille, INSERM, Lille, France
[3] Aix Marseille Univ, CNRS, Bioénergétique et Ingénierie des Protéines, IMM, Marseille, France
[4] Université Paris Cité, CNRS, IBPC, Laboratoire de Biochimie Théorique
[∗] first-authors
[#] corresponding authors



The folder pCALVADOS/ contains a notebook to run a pCALVADOS simulation on a local conda environment and a CSV file containing the pCALVADOS parameters (resulting from the incorporation of the parameters of dianionic phosphoresidues by Perdikari et al [1] and the CALVADOS 2 parameters by Tesei et al [2])

The folder Scripts_CG_2_AA/ contains all the scripts necessary to convert a coarse-grained pCALVADOS simulation into a minimized all-atom simulation


All data, molecular dynamics simulations and scripts are available in the following Zenodo repository : 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14900596.svg)](https://doi.org/10.5281/zenodo.14900596)



References :

[1] Perdikari TM, Jovic N, Dignon GL, Kim YC, Fawzi NL, Mittal J. A predictive coarse-grained model for position-specific effects of post-translational modifications. Biophys J. 2021 Apr 6;120(7):1187-1197. doi: 10.1016/j.bpj.2021.01.034.

[2] Tesei G, Lindorff-Larsen K. Improved predictions of phase behaviour of intrinsically disordered proteins by tuning the interaction range. Open Res Eur. 2023 Jan 17;2:94. doi: 10.12688/openreseurope.14967.2.

