# MODULE_Proteic_Menger_Curvature.py

"""
Description:
    This module provides functions for calculating the Proteic Menger Curvature (PMC).
    It includes functions for calculating the Euclidean distance between 2 3D points and Menger curvature
    based on three 3D points. Finally, the function PMC computes the Menger curvature for
    residues in a trajectory given a certain spacing.
    Additional informations and practical usage can be obtained in the article :
    "nP-collabs: Investigating counterion mediated bridges in the multiply phosphorylated tau-R2 repeat"
    by Jules Marien [1,2], Chantal Prévost [1,2], and Sophie Sacquin-Mora [1,2]

    [1] CNRS, Université de Paris, UPR 9080, Laboratoire de Biochimie Théorique, 
        13 rue Pierre et Marie Curie, 75005 Paris, France 
    [2] Institut de Biologie Physico-Chimique-Fondation Edmond de Rotschild, PSL Research University, 
        13 rue Pierre et Marie Curie, 75005 Paris, France 

Authors:
    Jules Marien, Chantal Prévost, and Sophie Sacquin-Mora

Date:
    11th of April 2024

Dependencies:
    - numpy
    - MDAnalysis

Usage:
    To use the functions in this module, import them into your Python script:
    
    Example:

    import MODULE_Proteic_Menger_Curvature
    
    result = MODULE_Proteic_Menger_Curvature.PMC('coordinates_file.pdb', 'trajectory.dcd', 'protein and name CA', 2)
    print(result)
    np.savetxt("result.txt", result)

"""





import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patches as patches





def dist(array_1, array_2):
    """Calculate the Euclidean distance between two 3D points."""
    return np.sqrt((array_1[0] - array_2[0])**2 + (array_1[1] - array_2[1])**2 + (array_1[2] - array_2[2])**2)






def menger_curvature(coord_point_1, coord_point_2, coord_point_3):
    """
    Calculate the Menger curvature of three points.

    Menger curvature is defined as 4A / |x-y||y-z||z-x|, where A is the area
    enclosed by the triangle formed by the three points and x, y and z
    are the coordinates of the summits of the triangle.

    Args:
        coord_point_1 (3D numpy.ndarray): Coordinates of the first summit.
        coord_point_2 (3D numpy.ndarray): Coordinates of the second summit.
        coord_point_3 (3D numpy.ndarray): Coordinates of the third summit.

    Returns:
        float: Menger curvature value.
    """
    # Calculation of the norm of the sides
    dist_1_2 = dist(coord_point_1, coord_point_2)
    dist_2_3 = dist(coord_point_2, coord_point_3)
    dist_3_1 = dist(coord_point_3, coord_point_1)

    # Calculation of area A
    p = (dist_1_2 + dist_2_3 + dist_3_1) / 2
    A = np.sqrt(p * (p - dist_1_2) * (p - dist_2_3) * (p - dist_3_1))

    # Calculation of Menger curvature
    menger_curvature_val = (4 * A) / (dist_1_2 * dist_2_3 * dist_3_1)

    return menger_curvature_val





def PMC(coordinates_file_name, trajectory_file_name, selection_c_alphas, spacing):
    """
    Calculate the Proteic Menger Curvature for residues in [spacing+1:N-spacing] in a trajectory where N is the total number of residues (counting from 1 to N).

    Args:
        coordinates_file_name (str): Name of the coordinates file.
        trajectory_file_name (str): Name of the trajectory file.
        selection_c_alphas (str): Selection string for C-alpha atoms.
        spacing (int): Spacing between residues for curvature calculation.

    Returns:
        numpy.ndarray: Array of Menger curvature values over time per residue. 
    """
    u = mda.Universe(coordinates_file_name, trajectory_file_name)
    len_traj = len(u.trajectory)
    u.trajectory[0]  # set u to first frame

    # initialize array cumulative displacement
    curvature_over_time_per_residue = np.empty((len(u.trajectory), len(u.select_atoms(selection_c_alphas).positions[:, 0])))

    # Get number of residues
    c_alphas = u.select_atoms(selection_c_alphas)
    residues = c_alphas.resids
    number_residues = len(residues)

    for t in range(len_traj):
        # Update trajectory to each frame
        u.trajectory[t]

        # Get positions at time t
        c_alphas_positions_t = u.select_atoms(selection_c_alphas).positions

        # Loop over residues
        for resid in range(spacing, number_residues - spacing):
            curvature_over_time_per_residue[t, resid] = menger_curvature(
                c_alphas_positions_t[resid - spacing, :],
                c_alphas_positions_t[resid, :],
                c_alphas_positions_t[resid + spacing, :]
            )

    return curvature_over_time_per_residue






def load_PMC(file_curvatures):
    """
    Load PMC data from a file.

    Args:
        file_curvatures (str): Path to the file containing curvatures data.

    Returns:
        tuple: A tuple containing:
            - numpy.ndarray: PMC data.
            - int: Dimension of the curvature matrix along the x-axis.
            - int: Dimension of the curvature matrix along the y-axis.
    """
    # Get the shape of the matrix thanks to the header
    with open(file_curvatures, "r") as f:
        header = f.readline()

    dimension_x_curvatures = int(header.split()[-3])
    dimension_y_curvatures = int(header.split()[-1])

    # Load the data
    curvatures = np.loadtxt(file_curvatures)
    curvatures = curvatures.reshape(dimension_x_curvatures, dimension_y_curvatures)

    return curvatures, dimension_x_curvatures, dimension_y_curvatures





def LC(PMCs):
    """
    Calculate Local Curvatures.

    Args:
        PMCs (numpy.ndarray): Array containing the Protein Menger Curvatures per conformation per residue.

    Returns:
        numpy.ndarray: Array containing the Local Curvatures per residue.
    """
    #Calculate Local Curvatures
    Local_Curvatures = np.mean(PMCs, axis=0)

    return Local_Curvatures






def LF(PMCs):
    """
    Calculate Local Flexibilities.

    Args:
        PMCs (numpy.ndarray): Array containing the Protein Menger Curvatures per conformation per residue.

    Returns:
        numpy.ndarray: Array containing the Local Flexibilities per residue.
    """
    #Calculate Local Flexibilities
    Local_Flexibilities = np.std(PMCs, axis=0)

    return Local_Flexibilities





### Plotting functions for demonstration purposes


def timeplot_PMC_3_replicas(time, curvatures_array_replica1, curvatures_array_replica2, curvatures_array_replica3, file_name_residues, size_figure_tuple, xticks, xlim_min, xlim_max, xlabel, ylabel, size_xlabel, size_ylabel, size_xticks, size_yticks, v_min, v_max, wspace, cmap, label_colorbar, fontsize_colorbar, file_name_figure, dpi):
    """
    Plot PMC time evolution for 3 replicas.

    Args:
        time (array-like): Time array.
        curvatures_array_replica1 (array-like): Curvature array for replica 1.
        curvatures_array_replica2 (array-like): Curvature array for replica 2.
        curvatures_array_replica3 (array-like): Curvature array for replica 3.
        file_name_residues (str): Name of the file containing residue names.
        size_figure_tuple (tuple): Figure size as a tuple (width, height) in inches.
        xticks (array-like): Array of x-tick locations.
        xlim_min (float): Minimum value for x-axis limit.
        xlim_max (float): Maximum value for x-axis limit.
        xlabel (str): Label for the x-axis.
        ylabel (str): Label for the y-axis.
        size_xlabel (int): Font size for the x-axis label.
        size_ylabel (int): Font size for the y-axis label.
        size_xticks (int): Font size for the x-ticks.
        size_yticks (int): Font size for the y-ticks.
        v_min (float): Minimum value for color scale.
        v_max (float): Maximum value for color scale.
        wspace (float): Width reserved for space between subplots.
        cmap (str): Colormap name.
        label_colorbar (str): Label for the color bar.
        fontsize_colorbar (int): Font size for the color bar label.
        file_name_figure (str): Name of the output figure file.
        dpi (int): Dots per inch for the output figure.

    Returns:
        3 PMC time evolutions
    """
    # Get the names of residues
    name_list_residues = file_name_residues
    with open(name_list_residues) as list_residues_file:
        list_residues_raw = list_residues_file.readlines()

    list_residues = [elem.strip('\n') for elem in list_residues_raw]

    # Plot
    fig, ax = plt.subplots(1, 3, sharey=True, figsize=size_figure_tuple, gridspec_kw={'wspace': wspace})

    xlabel_fontsize = ylabel_fontsize = size_xlabel

    ax[0].set_xlabel(xlabel, fontsize=xlabel_fontsize)
    ax[1].set_xlabel(xlabel, fontsize=xlabel_fontsize)
    ax[2].set_xlabel(xlabel, fontsize=xlabel_fontsize)

    ax[0].set_ylabel(ylabel, fontsize=ylabel_fontsize)

    array_residus = np.linspace(0, len(list_residues), len(list_residues))

    ax[0].yaxis.set_major_locator(mticker.FixedLocator(array_residus))
    ax[0].yaxis.set_major_formatter(mticker.FixedFormatter(list_residues))

    ax[0].yaxis.get_label().set_fontsize(size_yticks)

    xticks_fontsize = size_xticks

    for a in ax:
        a.set_xticks(xticks)
        a.set_xticklabels(xticks, fontsize=xticks_fontsize)
        a.set_xlim([xlim_min, xlim_max])
        a.set_ylim([-0.5, len(list_residues) + 0.5])

    ax[0].pcolormesh(time, array_residus, np.transpose(curvatures_array_replica1), cmap=cmap, shading='nearest', vmin=v_min, vmax=v_max)
    ax[1].pcolormesh(time, array_residus, np.transpose(curvatures_array_replica2), cmap=cmap, shading='nearest', vmin=v_min, vmax=v_max)
    im = ax[2].pcolormesh(time, array_residus, np.transpose(curvatures_array_replica3), cmap=cmap, shading='nearest', vmin=v_min, vmax=v_max)

    cbar = fig.colorbar(im, ax=ax.ravel().tolist())

    cbar.set_label(label=label_colorbar, fontsize=fontsize_colorbar, weight='bold')

    plt.savefig(file_name_figure, dpi=dpi)













# Create custom functions for plotting of Local Curvatures and Local Flexibilities

def plot_local_curvature_i_plus_5(LC_i_plus_5, file_name, number_first_residue, number_last_residue, spacing, size_figure_tuple, color_peptide_i_plus_5, ylim_min, ylim_max, xlabel, ylabel, size_xlabel, size_ylabel, size_xticks, size_yticks, linewidth, size_scatter, alpha_line, DPI, name_figure):
    """
    Plot local curvature for i+5 residues.

    Args:
        LC_i_plus_5 (numpy.ndarray): Array containing local curvatures for i+5 residues.
        file_name (str): Path to the file containing residue names.
        number_first_residue (int): Index of the first residue.
        number_last_residue (int): Index of the last residue.
        spacing (int): Spacing for range of residues.
        size_figure_tuple (tuple): Size of the figure (width, height) in inches.
        color_peptide_i_plus_5 (str): Color for the i+5 peptide.
        ylim_min (float): Minimum value for y-axis limit.
        ylim_max (float): Maximum value for y-axis limit.
        xlabel (str): Label for the x-axis.
        ylabel (str): Label for the y-axis.
        size_xlabel (int): Font size for the x-axis label.
        size_ylabel (int): Font size for the y-axis label.
        size_xticks (int): Font size for the x-ticks.
        size_yticks (int): Font size for the y-ticks.
        linewidth (float): Width of the plot lines.
        size_scatter (int): Size of the scatter points.
        alpha_line (float): Transparency of the plot lines.
        DPI (int): Dots per inch for the output figure.
        name_figure (str): Name of the output figure.

    Returns:
        Figure of the Local Curvature per residue
    """
    # Open file and extract residue names
    with open(file_name) as list_names_file:
        list_names = [elem.strip('\n') for elem in list_names_file.readlines()]

    # Create a range to plot
    range_residues = range(number_first_residue, number_last_residue + 1)

    # Plot
    plt.figure(figsize=size_figure_tuple)
    ax = plt.axes()

    plt.xlim(number_first_residue - 0.5, number_last_residue + 0.5)
    plt.ylim(ylim_min, ylim_max)

    plt.xlabel(xlabel, fontsize=size_xlabel)
    plt.ylabel(ylabel, fontsize=size_ylabel)

    # Add phospho rectangles
    rect_phospho_1 = patches.Rectangle((3.5, 0), 1, ylim_max, color='red', alpha=0.3)
    rect_phospho_2 = patches.Rectangle((8.5, 0), 1, ylim_max, color='red', alpha=0.3)

    ax.add_patch(rect_phospho_1)
    ax.add_patch(rect_phospho_2)

    # Plot curvature data
    plt.scatter(range_residues[spacing:-spacing], LC_i_plus_5[spacing:-spacing], color=color_peptide_i_plus_5, s=size_scatter, label='R+5')
    plt.plot(range_residues[spacing:-spacing], LC_i_plus_5[spacing:-spacing], color=color_peptide_i_plus_5, linewidth=linewidth, alpha=alpha_line)

    ax.tick_params('x', labelrotation=90)

    ax.xaxis.set_major_locator(mticker.FixedLocator(range_residues))
    ax.xaxis.set_major_formatter(mticker.FixedFormatter(list_names))

    plt.xticks(fontsize=size_xticks)
    plt.yticks(fontsize=size_yticks)

    plt.legend()

    plt.tight_layout()

    plt.savefig(name_figure, dpi=DPI)





def plot_local_flexibility_i_plus_5(LF_i_plus_5, file_name, number_first_residue, number_last_residue, spacing, size_figure_tuple, color_peptide_i_plus_5, ylim_min, ylim_max, xlabel, ylabel, size_xlabel, size_ylabel, size_xticks, size_yticks, linewidth, size_scatter, alpha_line, DPI, name_figure):
    """
    Plot local flexibility for i+5 residues.

    Args:
        LF_i_plus_5 (numpy.ndarray): Array containing local flexibility for i+5 residues.
        file_name (str): Path to the file containing residue names.
        number_first_residue (int): Index of the first residue.
        number_last_residue (int): Index of the last residue.
        spacing (int): Spacing for range of residues.
        size_figure_tuple (tuple): Size of the figure (width, height) in inches.
        color_peptide_i_plus_5 (str): Color for the i+5 peptide.
        ylim_min (float): Minimum value for y-axis limit.
        ylim_max (float): Maximum value for y-axis limit.
        xlabel (str): Label for the x-axis.
        ylabel (str): Label for the y-axis.
        size_xlabel (int): Font size for the x-axis label.
        size_ylabel (int): Font size for the y-axis label.
        size_xticks (int): Font size for the x-ticks.
        size_yticks (int): Font size for the y-ticks.
        linewidth (float): Width of the plot lines.
        size_scatter (int): Size of the scatter points.
        alpha_line (float): Transparency of the plot lines.
        DPI (int): Dots per inch for the output figure.
        name_figure (str): Name of the output figure.

    Returns:
        Figure of the Local Flexibility per residue
    """
    # Open file and extract residue names
    with open(file_name) as list_names_file:
        list_names = [elem.strip('\n') for elem in list_names_file.readlines()]

    # Create a range to plot
    range_residues = range(number_first_residue, number_last_residue + 1)

    # Plot
    plt.figure(figsize=size_figure_tuple)
    ax = plt.axes()

    plt.xlim(number_first_residue - 0.5, number_last_residue + 0.5)
    plt.ylim(ylim_min, ylim_max)

    plt.xlabel(xlabel, fontsize=size_xlabel)
    plt.ylabel(ylabel, fontsize=size_ylabel)

    # Add phospho rectangles
    rect_phospho_1 = patches.Rectangle((3.5, 0), 1, ylim_max, color='red', alpha=0.3)
    rect_phospho_2 = patches.Rectangle((8.5, 0), 1, ylim_max, color='red', alpha=0.3)

    ax.add_patch(rect_phospho_1)
    ax.add_patch(rect_phospho_2)

    # Plot flexibility data
    plt.scatter(range_residues[spacing:-spacing], LF_i_plus_5[spacing:-spacing], color=color_peptide_i_plus_5, s=size_scatter, label='R+5')
    plt.plot(range_residues[spacing:-spacing], LF_i_plus_5[spacing:-spacing], color=color_peptide_i_plus_5, linewidth=linewidth, alpha=alpha_line)

    ax.tick_params('x', labelrotation=90)

    ax.xaxis.set_major_locator(mticker.FixedLocator(range_residues))
    ax.xaxis.set_major_formatter(mticker.FixedFormatter(list_names))

    plt.xticks(fontsize=size_xticks)
    plt.yticks(fontsize=size_yticks)

    plt.legend()

    plt.tight_layout()

    plt.savefig(name_figure, dpi=DPI)
