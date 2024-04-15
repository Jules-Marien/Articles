# Nonbonded parameters used for Figure SI-10

This repository contains the nonbonded parameters of the ions (Na+ and K+), the water oxygens and the terminal phosphate oxygens used to plots the theoretical nonbonded potentials in Figure SI-10.

The interaction potential is the sum of the Lennard-Jones and the Coulomb potentials following the GROMACS documentation : https://manual.gromacs.org/current/reference-manual/functions/nonbonded-interactions.html 

For the systems using a 4-points water model, the position of the partial charge $q$ was assumed to be 0.0159nm and 0.0155nm further away for A19 and A99 respectively.

The Lennard-Jones potential is of the form $4\epsilon_{ij} [ (\frac{\sigma_{ij}}{r})^{12} -  (\frac{\sigma_{ij}}{r})^{6} ]$. 

The Coulomb potential is of the form $f \frac{q_i q_j}{\epsilon_{r} r}$ where $f = \frac{1}{4\pi\epsilon_0} = {138.935 458}$ kJ.mol⁻1.nm.e⁻2 is the electric conversion factor.

$\epsilon_{r}$ is set to 1.

$\epsilon_{ij}$ and $\sigma_{ij}$ are calculated following Lorentz-Berthelot combination rules.
