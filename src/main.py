from openmmtools import cache, testsystems
from openmmtools.states import SamplerState, ThermodynamicState
from simtk import unit

alanine = testsystems.AlanineDipeptideVacuum()
thermodynamic_state = ThermodynamicState(
    system=alanine.system, temperature=298.15*unit.kelvin)
sampler_state = SamplerState(positions=alanine.positions)

# TODO: replace with our own Move implementation
# see also: https://github.com/choderalab/openmmtools/blob/d4097719d49437c87ad8ca7136d74791fe762dcf/openmmtools/mcmc.py#L770
ghmc_move = GHMCMove(timestep=1.0*unit.femtosecond, n_steps=50)
sampler = MCMCSampler(thermodynamic_state, sampler_state, move=ghmc_move)

# run the sampler
sampler.minimize()
sampler.run(n_iterations=2) # TODO: increase iterations for final run

# get the results
results = sampler.sampler_state

# TODO: 
# - write results.(potential_energy, positions) to file, 
# - calculate the two angles we are interested in
# - plot them vs. the energies
