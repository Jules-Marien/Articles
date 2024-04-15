import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda 

import MODULE_Proteic_Menger_Curvature





##### Calculate PMC matrices

print("Calculating PMC for replica 1...")

Curvatures_spacing_2_peptide_i_plus_5_replica1 = MODULE_Proteic_Menger_Curvature.PMC(coordinates_file_name = "Trajectories/Peptide_R_plus_5/Replica1/gro_peptide_R_plus_5_replica1_renamed_serines.gro",
                                          trajectory_file_name = "Trajectories/Peptide_R_plus_5/Replica1/production_peptide_R_plus_5_replica1_protein_ions.xtc", 
                                          selection_c_alphas = "protein and name CA",
                                          spacing = 2 )


print("Calculating PMC for replica 2...")

Curvatures_spacing_2_peptide_i_plus_5_replica2 = MODULE_Proteic_Menger_Curvature.PMC(coordinates_file_name = "Trajectories/Peptide_R_plus_5/Replica2/gro_peptide_R_plus_5_replica2_renamed_serines.gro",
                                          trajectory_file_name = "Trajectories/Peptide_R_plus_5/Replica2/production_peptide_R_plus_5_replica2_protein_ions.xtc", 
                                          selection_c_alphas = "protein and name CA",
                                          spacing = 2 )


print("Calculating PMC for replica 3...")

Curvatures_spacing_2_peptide_i_plus_5_replica3 = MODULE_Proteic_Menger_Curvature.PMC(coordinates_file_name = "Trajectories/Peptide_R_plus_5/Replica3/gro_peptide_R_plus_5_replica3_renamed_serines.gro",
                                          trajectory_file_name = "Trajectories/Peptide_R_plus_5/Replica3/production_peptide_R_plus_5_replica3_protein_ions.xtc", 
                                          selection_c_alphas = "protein and name CA",
                                          spacing = 2 )



print("Done ! Saving data...")





##### Save data

np.savetxt("Data_curvature/Curvatures_spacing_2_peptide_i_plus_5_replica1.txt", Curvatures_spacing_2_peptide_i_plus_5_replica1, header='Shape of the matrix : 52001 x 13')
np.savetxt("Data_curvature/Curvatures_spacing_2_peptide_i_plus_5_replica2.txt", Curvatures_spacing_2_peptide_i_plus_5_replica2, header='Shape of the matrix : 52001 x 13')
np.savetxt("Data_curvature/Curvatures_spacing_2_peptide_i_plus_5_replica3.txt", Curvatures_spacing_2_peptide_i_plus_5_replica3, header='Shape of the matrix : 52001 x 13')


print('Done ! Plot PMCs...')





##### Load curvatures if necessary
"""
Curvatures_spacing_2_peptide_i_plus_5_replica1 = MODULE_Proteic_Menger_Curvature.load_PMC("Data_curvature/Curvatures_spacing_2_peptide_i_plus_5_replica1.txt")[0]
Curvatures_spacing_2_peptide_i_plus_5_replica2 = MODULE_Proteic_Menger_Curvature.load_PMC("Data_curvature/Curvatures_spacing_2_peptide_i_plus_5_replica2.txt")[0]
Curvatures_spacing_2_peptide_i_plus_5_replica3 = MODULE_Proteic_Menger_Curvature.load_PMC("Data_curvature/Curvatures_spacing_2_peptide_i_plus_5_replica3.txt")[0]
"""





##### Plot PMCs

MODULE_Proteic_Menger_Curvature.timeplot_PMC_3_replicas(time =  np.arange(0,52001)/100,   #Time in ns
                                                        curvatures_array_replica1 = Curvatures_spacing_2_peptide_i_plus_5_replica1, 
                                                        curvatures_array_replica2 = Curvatures_spacing_2_peptide_i_plus_5_replica2,
                                                        curvatures_array_replica3 = Curvatures_spacing_2_peptide_i_plus_5_replica3,
                                                        file_name_residues = "list_residues_peptide_i_plus_5.txt",
                                                        size_figure_tuple = (22,5),
                                                        xlim_min = 0,
                                                        xlim_max = 520,
                                                        xlabel = "Time (ns)",
                                                        ylabel = "Residues of toy model R+5",
                                                        xticks=[0,50,100,150,200,250,300,350,400,450,500],
                                                        size_xlabel = 15, 
                                                        size_ylabel = 18,
                                                        size_xticks = 14, 
                                                        size_yticks = 18, 
                                                        v_min = 0, 
                                                        v_max = 0.35, 
                                                        wspace = 0.025,
                                                        cmap = "Blues",
                                                        label_colorbar = "Proteic Menger Curvature ($Å^{⁻1}$)",
                                                        fontsize_colorbar = 13,
                                                        file_name_figure = "Plots/timeplot_PMCs_peptide_i_plus_5.png", 
                                                        dpi = 300)


print('Done ! Concatenate the PMCs of each replica...')





##### Concatenate the 3 replicas 

Curvatures_spacing_2_peptide_i_plus_5 = np.concatenate( (Curvatures_spacing_2_peptide_i_plus_5_replica1,
                                       Curvatures_spacing_2_peptide_i_plus_5_replica2, 
                                       Curvatures_spacing_2_peptide_i_plus_5_replica3), axis = 0)


print('Done ! Calculate Local Curvatures...')




##### Calculate Local Curvatures

Local_curvatures_spacing_2_peptide_i_plus_5 = MODULE_Proteic_Menger_Curvature.LC(Curvatures_spacing_2_peptide_i_plus_5)


print('Done ! Calculate Local Flexibilities...')





##### Calculate Local Flexibilities

Local_flexibilities_spacing_2_peptide_i_plus_5 = MODULE_Proteic_Menger_Curvature.LF(Curvatures_spacing_2_peptide_i_plus_5)


print("Done ! Saving data...")





##### Save data

np.savetxt("Data_curvature/Local_curvatures_spacing_2_peptide_i_plus_5.txt", Local_curvatures_spacing_2_peptide_i_plus_5, header='Shape of the array : 13')
np.savetxt("Data_curvature/Local_flexibilities_spacing_2_peptide_i_plus_5.txt", Local_flexibilities_spacing_2_peptide_i_plus_5, header='Shape of the array : 13')


print('Done ! Plot Local Curvatures...')





##### Plot Local Curvatures 

MODULE_Proteic_Menger_Curvature.plot_local_curvature_i_plus_5(LC_i_plus_5 = Local_curvatures_spacing_2_peptide_i_plus_5,
                              file_name = "list_residues_peptide_i_plus_5.txt",
                              number_first_residue = 1,
                              number_last_residue = 13,
                              spacing = 2,
                              size_figure_tuple = (9, 6),
                              color_peptide_i_plus_5 = "#0091f7",
                              ylim_min = 0,
                              ylim_max = 0.25,
                              xlabel = "Residues of toy model R+5",
                              ylabel = "Local Curvature ($Å^{⁻1}$)",
                              size_xlabel = 18,
                              size_ylabel = 18,
                              size_xticks = 14, 
                              size_yticks = 14,
                              linewidth = 2,
                              size_scatter = 30,
                              alpha_line = 0.6,
                              DPI = 300,
                              name_figure = "Plots/Local_Curvatures.png")


print('Done ! Plot Local Curvatures...')





##### Plot Local Flexibilites 

# Call function
MODULE_Proteic_Menger_Curvature.plot_local_flexibility_i_plus_5(LF_i_plus_5 = Local_flexibilities_spacing_2_peptide_i_plus_5,
                              file_name = "list_residues_peptide_i_plus_5.txt",
                              number_first_residue = 1,
                              number_last_residue = 13,
                              spacing = 2,
                              size_figure_tuple = (9, 6),
                              color_peptide_i_plus_5 = "#0091f7",
                              ylim_min = 0,
                              ylim_max = 0.10,
                              xlabel = "Residues of toy model R+5",
                              ylabel = "Local Flexibility ($Å^{⁻1}$)",
                              size_xlabel = 18,
                              size_ylabel = 18,
                              size_xticks = 14, 
                              size_yticks = 14,
                              linewidth = 2,
                              size_scatter = 30,
                              alpha_line = 0.6,
                              DPI = 300,
                              name_figure = "Plots/Local_Flexibilities.png")

print('Done ! Script finished !')

