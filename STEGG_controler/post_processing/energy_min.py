import os
import shutil
import numpy as np
import pandas as pd
from sys import stdout
import sys
import copy
from time import sleep

from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue

from openmm.app import *
from openmm import *
from openmm.unit import *

from pdbfixer import PDBFixer
from openmm.app import PDBFile
from biopandas.pdb import PandasPdb
import mdtraj

sys.path.insert(1,'/home/STEGG_controler/')
from utils import *

def minimizeConf(pdb_filename,device='CPU',with_solvent=False):

    # fix pdb
    print('fixing pdb file')
    fixer = PDBFixer(filename=pdb_filename)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    with open(pdb_filename[:-4]+'_fixed.pdb', 'w') as out_file:
        PDBFile.writeFile(fixer.topology, fixer.positions, out_file, keepIds=True)

    # Read PDB
    print('reading new pdb')
    pdb = PDBFile(pdb_filename[:-4]+'_fixed.pdb')
    top = pdb.getTopology()
    positions = np.array(pdb.positions)
    numAtoms = len(positions)
    positions = np.reshape(positions, (3*numAtoms,1))

    # Create forcefield
    print('creating forcefield, integrators, and simulation')
    if with_solvent:
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    else:
        forcefield = ForceField('amber14-all.xml')    
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.deleteWater()
    residues = modeller.addHydrogens(forcefield)
    if with_solvent:
        modeller.addSolvent(forcefield, padding=1.0*nanometer) # can keep explicit solvent and rigid waters
    # system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds) #cutoff refers to what molecules influence others (far away)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=CutoffNonPeriodic, constraints=HBonds)

    # Rest (Integrator + Platform)
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds) # normaly use 0.002*picoseconds for step size (to keep chains from drifting apart)
    # can also try including all of the TCR or covalently bindind cystines
    platform = Platform.getPlatformByName(device)

    # Create Simulation
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)

    print('minimizing energy')
    simulation.minimizeEnergy()
    energy_val = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    r = PDBReporter(pdb_filename[:-4]+'_minimized.pdb', 1)
    r.report(simulation, simulation.context.getState(getPositions=True, getEnergy=True))

    print('min_energy: ', energy_val)
    # print('running short simulation')
    # simulation.reporters.append(PDBReporter('output_conformations.pdb', 1000))
    # simulation.reporters.append(StateDataReporter('stats.txt', 1000, step=True,
    #         potentialEnergy=True, temperature=True))
    # simulation.step(10000)

    return energy_val
