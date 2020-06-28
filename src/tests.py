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
from simtk.openmm import *
from simtk.openmm.app import *

# simtk.testInstallation.main()

# setup our system
pdb = PDBFile(os.path.dirname(os.path.realpath(__file__)) +
              '/input/alanine-dipeptide-nowater-noneq.pdb')
forcefield = ForceField('amber10.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
                                 nonbondedCutoff=1*unit.nanometer, constraints=HBonds)
thermodynamic_state = ThermodynamicState(
    system=system, temperature=298.15*unit.kelvin)
positions = pdb.getPositions()

# finally, start with sampling
sampler_state = SamplerState(positions=positions)

# run sampler multiple times with different parameters
# in order to check for differences due to different timesteps
# (runtime can only be compared for same acceptance rate)
headerText = "Timestepsize\tAccepted\tProposed\tTime [s]"
accept_state = [True, False]
for accept in accept_state:
    text = "Accepting anyway: {}".format(accept)
    print(headerText)
    f = open("out/test_results_perf_" + str(accept) + ".tsv", mode='w+')
    print(headerText, file=f)
    # run different nr. of timesteps to go before recalculation
    for n_steps in range(1, 10):
        V_calculator = TorsionNeglectingPotentialEnergyCalculator(
            timesteps=n_steps
        )
        move = DeterministicRotateDisplaceMove(
            displacement_sigma=0.01*unit.nanometer,
            atom_subset=1,
            potential_energy_calculator=V_calculator,
            accept_anyway=accept,
            rng_seed=42
        )
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
        print(text, file=f)

        # write the results
        # only interesting if they are not all the same
        if (not accept):
            out_filename = os.path.dirname(os.path.realpath(
                __file__)) + "/../out/test_results_" + str(accept) + "_" + str(n_steps) + ".json"
            move.write_to_file(out_filename=out_filename)
