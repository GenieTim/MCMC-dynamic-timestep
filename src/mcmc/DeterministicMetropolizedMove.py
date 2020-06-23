# =============================================================================
# METROPOLIZED MOVE BASE CLASS
# =============================================================================
# As adapted from https://github.com/choderalab/openmmtools/blob/d4097719d49437c87ad8ca7136d74791fe762dcf/openmmtools/mcmc.py#L770
import abc
import copy
import json
import os

import numpy as np
from numpy.random import default_rng, choice

from openmmtools import cache, integrators, utils
from openmmtools.mcmc import MCMCMove
from openmmtools.utils import SubhookedABCMeta, Timer


class DeterministicMetropolizedMove(MCMCMove):
    """A base class for metropolized moves with a custom pot energy calculator
    and .
    This class is intended to be inherited by MCMCMoves that needs to
    accept or reject a proposed move with a Metropolis criterion. Only
    the proposal needs to be specified by subclasses through the method
    _propose_positions(), the way to calculate potential energy through 
    _calculate_potential_energy.

    Parameters
    ----------
    accept_anyway : bool
        Indicate that you don't care about the result of the potential energy calculation, 
        you want to accept all moves. This is useful in case you want to compare two methods,
        to get an overall sum of how many moves are accepted.
    atom_subset : slice or list of int or int, optional
        If specified, the move is applied only to those atoms specified by these
        indices. If int, move is applied to this many random atoms.
        If None, the move is applied to all atoms (default is None).
    context_cache : openmmtools.cache.ContextCache, optional
        The ContextCache to use for Context creation. If None, the global cache
        openmmtools.cache.global_context_cache is used (default is None).
    potential_energy_calculator : an object with a function 
        _calculate_potential_energy(self, thermodynamic_state, context)

    Attributes
    ----------
    n_accepted : int
        The number of proposals accepted.
    n_proposed : int
        The total number of attempted moves.
    accept_anyway
    atom_subset
    context_cache
    potential_energy_calculator

    Examples
    --------
    >>> from simtk import unit
    >>> from openmmtools import testsystems, states
    >>> class AddOneVector(CustomEnergiedMetropolizedMove):
    ...     def __init__(self, **kwargs):
    ...         super(AddOneVector, self).__init__(**kwargs)
    ...     def _propose_positions(self, initial_positions):
    ...         print('Propose new positions')
    ...         displacement = unit.Quantity(np.array([1.0, 1.0, 1.0]), initial_positions.unit)
    ...         return initial_positions + displacement
    ...
    >>> alanine = testsystems.AlanineDipeptideVacuum()
    >>> sampler_state = states.SamplerState(alanine.positions)
    >>> thermodynamic_state = states.ThermodynamicState(alanine.system, 300*unit.kelvin)
    >>> move = AddOneVector(atom_subset=list(range(sampler_state.n_particles)))
    >>> move.apply(thermodynamic_state, sampler_state)
    Propose new positions
    >>> move.n_accepted
    1
    >>> move.n_proposed
    1
    """

    def __init__(self, potential_energy_calculator=None, atom_subset=None, context_cache=None, accept_anyway=False, rng_seed=42):
        self.n_accepted = 0
        self.n_proposed = 0
        self.accept_anyway = accept_anyway
        self.potential_energy_calculator = potential_energy_calculator
        self.atom_subset = atom_subset
        self.context_cache = context_cache
        self.state_history = []
        # use "deterministic" random generator
        self.rng = default_rng(rng_seed)

    @property
    def statistics(self):
        """The acceptance statistics as a dictionary."""
        return dict(n_accepted=self.n_accepted, n_proposed=self.n_proposed)

    @statistics.setter
    def statistics(self, value):
        self.n_accepted = value['n_accepted']
        self.n_proposed = value['n_proposed']

    # MARK: this is the function to adapt to our ideas
    def _calculate_potential_energy(self, thermodynamic_state, context):
        return self.potential_energy_calculator._calculate_potential_energy(thermodynamic_state, context)

    def apply(self, thermodynamic_state, sampler_state):
        """Apply a metropolized move to the sampler state.
        Total number of acceptances and proposed move are updated.
        Parameters
        ----------
        thermodynamic_state : openmmtools.states.ThermodynamicState
           The thermodynamic state to use to apply the move.
        sampler_state : openmmtools.states.SamplerState
           The initial sampler state to apply the move to. This is modified.
        """
        timer = Timer()
        benchmark_id = 'Applying {}'.format(self.__class__.__name__)
        timer.start(benchmark_id)

        # Check if we have to use the global cache.
        if self.context_cache is None:
            context_cache = cache.global_context_cache
        else:
            context_cache = self.context_cache

        # Create context, any integrator works.
        context, unused_integrator = context_cache.get_context(
            thermodynamic_state)

        # Compute initial energy. We don't need to set velocities to compute the potential.
        sampler_state.apply_to_context(context, ignore_velocities=True)
        initial_energy = self._calculate_potential_energy(
            thermodynamic_state, context)

        # Handle default and weird cases for atom_subset.
        if self.atom_subset is None:
            atom_subset = slice(None)
        elif isinstance(self.atom_subset, int): # we've been given a number of atoms to move randomly
            if(self.atom_subset < 2):  
                #the second case, wanting too many atoms, will cause random.choice to raise an error anyway
                raise ValueError("We tried to pick too few random atoms.")
            atom_subset = self.rng.choice(sampler_state.n_particles, self.atom_subset, replace = False)
            #choses self.atom_subset numbers from the range 0 to n_particles, without replacements
        #I am not sure we need to care about edge cases like this, commented out right now because it caused issues
        #elif not isinstance(self.atom_subset, slice) and len(self.atom_subset) == 1:
            # Slice so that initial_positions (below) will have a 2D shape.
        #    atom_subset = slice(self.atom_subset[0], self.atom_subset[0]+1)
        else:
            atom_subset = self.atom_subset

        # Store initial positions of the atoms that are moved.
        # We'll use this also to recover in case the move is rejected.
        if isinstance(atom_subset, slice):
            # Numpy array when sliced return a view, they are not copied.
            initial_positions = copy.deepcopy(
                sampler_state.positions[atom_subset])
        else:
            # This automatically creates a copy.
            # If random particles were picked, then atom_subset is a list of integers. 
            initial_positions = sampler_state.positions[atom_subset]

        # Propose perturbed positions. Modifying the reference changes the sampler state.
        proposed_positions = self._propose_positions(initial_positions)

        # Compute the energy of the proposed positions.
        sampler_state.positions[atom_subset] = proposed_positions
        sampler_state.apply_to_context(context, ignore_velocities=True)
        proposed_energy = self._calculate_potential_energy(
            thermodynamic_state, context)

        # Accept or reject with Metropolis criteria.
        delta_energy = proposed_energy - initial_energy
        if (not np.isnan(proposed_energy) and
                (delta_energy <= 0.0 or self.rng.random() < np.exp(-delta_energy))):
            self.n_accepted += 1
        elif (not self.accept_anyway):
            # Restore original positions.
            sampler_state.positions[atom_subset] = initial_positions
        self.n_proposed += 1

        # Print timing information.
        timer.stop(benchmark_id)
        # timer.report_timing()
        # remember state so we can write to file later in order to create trajectory
        # note that this uses very much memory – ok for our small molecule,
        # but definintely not price worthy
        self.state_history.append({
            'accepted': self.n_accepted,
            'proposed': self.n_proposed,
            'positions': {
                'value': sampler_state.positions.__getstate__()['_value'].tolist(),
                'unit': str(sampler_state.positions.__getstate__()['unit'].get_name())
            },
            # 'velocities': { # in this MC, velocities by sampler state will always be 0. they are not relevant for the angles anyway
            #     'value': sampler_state.velocities.__getstate__()['_value'].tolist(),
            #     'unit': str(sampler_state.velocities.__getstate__()['unit'].get_name())
            # },
            'initial_energy': initial_energy,
            'proposed_energy': proposed_energy,
            'unmodified_energy': thermodynamic_state.reduced_potential(context),
            'current_energy': {
                'unit': str(sampler_state.kinetic_energy.__getstate__()['unit'].get_name()) + ']',
                'value': copy.deepcopy(sampler_state.kinetic_energy.__getstate__()['_value'])
            },
            'timing': copy.deepcopy(timer.report_timing())
        })

    def write_to_file(self, out_filename=None):
        if (out_filename is None):
            out_filename = os.path.dirname(os.path.realpath(
                __file__)) + "/../../out/results_" + str(self.__class__.__name__) + ".json"
        output_file = open(out_filename, "w+")
        json.dump(self.state_history, fp=output_file)
        output_file.close()

    def __getstate__(self):
        if self.context_cache is None:
            context_cache_serialized = None
        else:
            context_cache_serialized = utils.serialize(self.context_cache)
        serialization = dict(atom_subset=self.atom_subset,
                             context_cache=context_cache_serialized)
        serialization.update(self.statistics)
        return serialization

    def __setstate__(self, serialization):
        self.atom_subset = serialization['atom_subset']
        if serialization['context_cache'] is None:
            self.context_cache = None
        else:
            self.context_cache = utils.deserialize(
                serialization['context_cache'])
        self.statistics = serialization

    @abc.abstractmethod
    def _propose_positions(self, positions):
        """Return new proposed positions.
        These method must be implemented in subclasses.
        Parameters
        ----------
        positions : nx3 numpy.ndarray
            The original positions of the subset of atoms that these move
            applied to.
        Returns
        -------
        proposed_positions : nx3 numpy.ndarray
            The new proposed positions.
        """
        pass
