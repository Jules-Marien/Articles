#!/bin/bash

for i in $(seq 16000 17999)
do
  namd2 --structure modified_frame_${i}.psf --coordinates modified_frame_${i}.pdb --outputName minimization_frame_${i} --consref modified_frame_${i}_restraints.ref --conskfile modified_frame_${i}_restraints.ref SCRIPT_minimization_input_file.inp > minimization_frame_${i}.log
  rm minimization_frame_${i}.coor
  rm minimization_frame_${i}.log
  rm minimization_frame_${i}.vel
  rm minimization_frame_${i}.xsc
  rm modified_frame_${i}.pdb 
  rm modified_frame_${i}.psf
  rm modified_frame_${i}_restraints.ref
  echo $i
done

