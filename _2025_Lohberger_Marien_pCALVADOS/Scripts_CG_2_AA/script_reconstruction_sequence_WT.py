import wget
import os
import shutil
import numpy as np
import pandas as pd
import scipy.stats as scs
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression
import mdtraj as md
#from simtk import openmm, unit
#from simtk.openmm import app
import matplotlib as mpl
import matplotlib.pyplot as plt
#import BME as BME
from kneed import KneeLocator
#from google.colab import files
import matplotlib_inline.backend_inline
matplotlib_inline.backend_inline.set_matplotlib_formats('pdf', 'svg')
from ipywidgets import interactive
import ipywidgets as widgets
import subprocess
import sys


#USER 

SEQUENCE = "MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPGSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL"

NAME = "Simu_CALVADOS_sequence_WT"






#PROGRAM

wget.download('https://raw.githubusercontent.com/Jules-Marien/CALVADOS_phosphorylations/main/residues_and_phosphoresidues.csv')

residues = pd.read_csv('residues_and_phosphoresidues.csv')
residues = residues.set_index('one')


def backmapping(traj, dt):
    for i in np.arange(0,len(t_cg),dt):
        print("Backmapping step ", i)
        t_cg[int(i)].save_pdb('frame.pdb')
        subprocess.run(['./pulchra', 'frame.pdb'])
        if i == 0:
            traj_AA = md.load_pdb('frame.rebuilt.pdb')
            shutil.move('frame.rebuilt.pdb', 'top_AA.pdb')
        else:
            traj_AA += md.load_pdb('frame.rebuilt.pdb')
    traj_AA.save_dcd('traj_AA.dcd')
    #shutil.move('frame.rebuilt.pdb', 'top_AA.pdb')
    return traj_AA



def fix_topology(t,seq):
    cgtop = md.Topology()
    cgchain = cgtop.add_chain()
    for res in seq:
        cgres = cgtop.add_residue(res, cgchain)
        cgtop.add_atom('CA', element=md.element.carbon, residue=cgres)
    traj = md.Trajectory(t.xyz, cgtop, t.time, t.unitcell_lengths, t.unitcell_angles)
    traj = traj.superpose(traj, frame=0)
    return traj







SEQ3 = [residues.three[x] for x in SEQUENCE]
#t_cg = fix_topology(md.load_dcd('{:s}/traj.dcd'.format(NAME), top='{:s}/top.pdb'.format(NAME)), SEQ3)

t_cg = fix_topology(md.load_dcd('traj.dcd', top='top.pdb'), SEQ3)


#print('Subsampling trajectory prior to backmapping to all-atom. Taken 1 simulation frame every {}'.format(sub_f))
print('Backmapping to all-atom resolution...')


stride = 1

traj_AA = backmapping(t_cg, stride)





