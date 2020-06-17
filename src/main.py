import json
import os

import pandas as pd

from mcmc.DeterministicRotateDisplaceMove import \
    DeterministicRotateDisplaceMove
from openmmtools import cache, testsystems
from openmmtools.mcmc import GHMCMove, MCMCSampler
from openmmtools.states import SamplerState, ThermodynamicState
from potentialEnergyCalculator.InitialPotentialEnergyCalculator import \
    InitialPotentialEnergyCalculator
from potentialEnergyCalculator.TorsionNeglectingPotentialEnergyCalculator import \
    TorsionNeglectingPotentialEnergyCalculator
from simtk import openmm, unit
from simtk.openmm import *
from simtk.openmm.app import *

# setup our system
# Two possible ways:
# a) use predifined system
# alanine = testsystems.AlanineDipeptideVacuum()
# system = alanine.system
# positions = alanine.positions

# b) create & load our own
pdb = PDBFile(os.path.dirname(os.path.realpath(__file__)) +
              '/input/alanine-dipeptide-nowater.pdb')
forcefield = ForceField('amber10.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
                                 nonbondedCutoff=1*unit.nanometer, constraints=HBonds)
thermodynamic_state = ThermodynamicState(
    system=system, temperature=298.15*unit.kelvin)
positions = pdb.getPositions()

# finally, start with sampling
sampler_state = SamplerState(positions=positions)

# TODO: decide on nr. of timesteps to go before recalculation
V_calculator = TorsionNeglectingPotentialEnergyCalculator(timesteps=5)
move = DeterministicRotateDisplaceMove(
    displacement_sigma=5.0*unit.nanometer, potential_energy_calculator=V_calculator)
sampler = MCMCSampler(thermodynamic_state, sampler_state, move=move)

# run the sampler
sampler.minimize()
sampler.run(n_iterations=400)  # TODO: increase iterations for final run

# write the results
move.write_to_file()

# get the results of the last/final iterations
# this is an example on how to do it in the actual, needed logger
results = sampler.sampler_state
velocities = results.velocities.__getstate__()
positions = results.positions.__getstate__()
k_energy = results.kinetic_energy.__getstate__()

data = {
    'velocities [' + str(velocities['unit'].get_name()) + ']': velocities['_value'].tolist(),
    'positions [' + str(positions['unit'].get_name()) + ']': positions['_value'].tolist(),
    'kinetic_energy [' + str(k_energy['unit'].get_name()) + ']': k_energy['_value']
}

out_filename = os.path.dirname(
    os.path.realpath(__file__)) + "/../out/results.csv"
f = open(out_filename, "w+")
json.dump(data, fp=f)
f.close()

# TODO:
# - read results.(potential_energy, positions) from file written by CustomEnergiedMetropolizedMove,
# - calculate the two angles we are interested in
# - plot them vs. the energies
