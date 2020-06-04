import json
import os

import pandas as pd

from mcmc.TorsionNeglectingMove import TorsionNeglectingMove
from openmmtools import cache, testsystems
from openmmtools.mcmc import GHMCMove, MCMCSampler
from openmmtools.states import SamplerState, ThermodynamicState
from simtk import unit

# setup our system
alanine = testsystems.AlanineDipeptideVacuum()
thermodynamic_state = ThermodynamicState(
    system=alanine.system, temperature=298.15*unit.kelvin)
sampler_state = SamplerState(positions=alanine.positions)

# TODO: decide on nr. of timesteps to go before recalculation
move = TorsionNeglectingMove(timesteps=5)
sampler = MCMCSampler(thermodynamic_state, sampler_state, move=move)

# run the sampler
sampler.minimize()
sampler.run(n_iterations=2)  # TODO: increase iterations for final run

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
