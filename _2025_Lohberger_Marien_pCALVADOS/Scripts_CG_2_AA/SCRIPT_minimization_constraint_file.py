#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 11:02:49 2023

@author: marien
"""



import numpy as np


nbre_lines_PDB_file = 6418
    



def PDB_ATOM_line_parser(line):
    """Takes an ATOM line as input and outputs its compenents based on the PDB format available at :
        https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html """
        
    ATOM = line[0:4]
    Space_1 = '  '
    Atom_number = line[6:11]
    Space_2 = ' '
    Atom_name = line[12:16]
    Alternate_location_indicator = line[16]
    Residue_name = line[17:20]
    Space_3 = ' '
    Chain_identifier = line[21]
    Residue_sequence_number = line[22:26]
    Code_insertion_residue = line[26]
    Space_4 = '   '
    X_coord = line[30:38]
    Y_coord = line[38:46]
    Z_coord = line[46:54]
    Occupancy = line[54:60]
    Temperature_factor = line[60:66]
    Space_5 = '      '
    Segment_identifier = line[72:76]
    Element_symbol = line[76:78]
    Charge = line[78:80]
    
    list_components = [ATOM,  #0
                       Space_1, #1
                       Atom_number, #2
                       Space_2, #3
                       Atom_name, #4
                       Alternate_location_indicator, #5
                       Residue_name, #6
                       Space_3, #7
                       Chain_identifier, #8
                       Residue_sequence_number, #9
                       Code_insertion_residue, #10
                       Space_4, #11
                       X_coord, #12
                       Y_coord, #13
                       Z_coord, #14
                       Occupancy, #15
                       Temperature_factor, #16
                       Space_5, #17
                       Segment_identifier, #18
                       Element_symbol, #19
                       Charge] #20
    
    return list_components
    



for file_number in range(0,20001):
    print(file_number)

    template = np.empty((nbre_lines_PDB_file),dtype=object)


    #Loading template into numpy array

    with open('modified_frame_'+ str(file_number) + '.pdb','r') as file:
        template = file.readlines()


    array_reorganization = np.empty((nbre_lines_PDB_file), dtype=object)  #Array which will contain the modified data 

    array_reorganization[0] = template[0] #Keeping the first and last lines intact
    array_reorganization[-1] = template[-1] #Keeping the first and last lines intact

    for i in range(1,nbre_lines_PDB_file-1):

        #print(i)

        line_elements = PDB_ATOM_line_parser(template[i])
    
        new_line_elements = line_elements.copy()
    
        new_atom_number = str(i)

        #Only put a restraint for the backbone 
        if new_line_elements[4] == " C  " or new_line_elements[4] == " CA " or new_line_elements[4] == " N  " :
            new_line_elements[15] = "  1.00" 

        else: 
            new_line_elements[15] = "  0.00"  
  
        array_reorganization[i] = ''.join(new_line_elements)




    #Create the output file

    with open('modified_frame_' + str(file_number) + '_restraints.ref','w') as outputfile:
        for i in range(len(array_reorganization)):
            outputfile.write(array_reorganization[i])



