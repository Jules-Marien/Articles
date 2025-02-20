topology topology_files/top_all36_prot.rtf
topology topology_files/par_all36m_prot.prm
topology topology_files/top_all36_na.rtf
topology topology_files/par_all36_na.prm
topology topology_files/toppar_all36_prot_na_combined.str
 
pdbalias atom ILE CD1 CD 
pdbalias atom GLY OXT OT1 

set number_frames [molinfo top get numframes]



for { set i 0} {$i <= $number_frames} {incr i} {

set sel [atomselect top all frame $i]


$sel writepdb temp_pdb.pdb


puts "First step"
segment P1 { pdb temp_pdb.pdb
first NTER
last CTER } 

puts "Second step"
#patch SP2 P1:201
#patch THPB P1:204
#patch THPB P1:211

puts "Third step"
regenerate angles dihedrals 

pdbalias atom ILE CD1 CD 
pdbalias atom GLY OXT OT1

puts "Fourth step"
coordpdb temp_pdb.pdb P1

puts "Fifth step" 
guesscoord 

puts "Sixth step"

writepdb Trajectory/modified_frame_$i.pdb 
writepsf Trajectory/modified_frame_$i.psf
resetpsf
}



