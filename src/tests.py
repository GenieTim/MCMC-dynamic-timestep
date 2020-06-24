import copy
import os
import time

import simtk.testInstallation
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

simtk.testInstallation.main()

# setup our system
alanine = testsystems.AlanineDipeptideVacuum()
thermodynamic_state = ThermodynamicState(
    system=alanine.system, temperature=298.15*unit.kelvin)
sampler_state = SamplerState(positions=alanine.positions)

# run sampler multiple times with different parameters
# in order to check for differences due to different timesteps
# (runtime can only be compared for same acceptance rate)
f = open("../out/test_results.txt", mode='w+')
accept_state = [True, False]
for accept in accept_state:
    text = "Accepting anyway: {}".format(accept)
    print(text)
    print(text, file = f)
    text = "Timestepsize\tAccepted\tProposed\tTime [s]"
    print(text)
    print(text, file = f)
    # run different nr. of timesteps to go before recalculation
    for n_steps in range(1, 10):
        V_calculator = TorsionNeglectingPotentialEnergyCalculator(
            timesteps=n_steps)
        move = DeterministicRotateDisplaceMove(
            displacement_sigma=5.0*unit.nanometer, atom_subset = 1, potential_energy_calculator=V_calculator, accept_anyway=accept, rng_seed = 42)
        sampler = MCMCSampler(copy.deepcopy(thermodynamic_state),
                              copy.deepcopy(sampler_state), move=move)
        sampler.minimize()
        # take time
        start = time.time()
        # TODO: increase iterations for final run
        sampler.run(n_iterations=400)
        end = time.time()

        text = "{}\t{}\t{}\t{}".format(n_steps, move.n_accepted,
                                      move.n_proposed, end - start)

        print(text)
        print(text, file = f)

        # write the results
        # only interesting if they are not all the same
        if (not accept):
            out_filename = os.path.dirname(os.path.realpath(
                __file__)) + "/../out/test_results_" + str(n_steps) + ".json"
            move.write_to_file(out_filename=out_filename)

# test with different number of atoms to move
# and different displacement sigma

print("Testing different atom numbers to move, accepting anyways: False")
print("Atmom Number\tAccepted\tProposed\tTime [s]")

for a_number in range(1, 10):
    V_calculator = TorsionNeglectingPotentialEnergyCalculator(timesteps=1)
    move = DeterministicRotateDisplaceMove(
        displacement_sigma=5.0*unit.nanometer, atom_subset = a_number, potential_energy_calculator=V_calculator, accept_anyway=False, rng_seed = 42)
    sampler = MCMCSampler(copy.deepcopy(thermodynamic_state), copy.deepcopy(sampler_state), move=move)
    sampler.minimize()
    # take time
    start = time.time()
    sampler.run(n_iterations=400)
    end = time.time()

    text = "{}\t{}\t{}\t{}".format(a_number, move.n_accepted,
                                      move.n_proposed, end - start)
    print(text)

print("Testing different displacement distance to move, accepting anyways: False. 1 Atom")
print("Sigma [nm]\tAccepted\tProposed\tTime [s]")

for sigma_factor in range(1, 11):
    V_calculator = TorsionNeglectingPotentialEnergyCalculator(timesteps=1)
    move = DeterministicRotateDisplaceMove(
        displacement_sigma=0.5*unit.nanometer*sigma_factor, atom_subset = 1, potential_energy_calculator=V_calculator, accept_anyway=False, rng_seed = 42)
    sampler = MCMCSampler(copy.deepcopy(thermodynamic_state), copy.deepcopy(sampler_state), move=move)
    sampler.minimize()
    # take time
    start = time.time()
    sampler.run(n_iterations=400)
    end = time.time()

    text = "{}\t{}\t{}\t{}".format(sigma_factor * 0.5, move.n_accepted,
                                      move.n_proposed, end - start)
    print(text)