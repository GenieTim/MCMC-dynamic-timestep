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

# TODO: replace with our own Move implementation
# see also: https://github.com/choderalab/openmmtools/blob/d4097719d49437c87ad8ca7136d74791fe762dcf/openmmtools/mcmc.py#L770
move = TorsionNeglectingMove(timesteps=5)
sampler = MCMCSampler(thermodynamic_state, sampler_state, move=move)

# run the sampler
sampler.minimize()
sampler.run(n_iterations=2)  # TODO: increase iterations for final run

# get the results
results = sampler.sampler_state
print(results)

# TODO:
# - write results.(potential_energy, positions) to file,
# - calculate the two angles we are interested in
# - plot them vs. the energies
