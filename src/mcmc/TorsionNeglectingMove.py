import numpy as np

from simtk import openmm, unit

from .RotateDisplaceMove import RotateDisplaceMove


# inspired by https://github.com/choderalab/openmmtools/blob/d4097719d49437c87ad8ca7136d74791fe762dcf/openmmtools/mcmc.py#L1757
# and https://github.com/choderalab/openmmtools/blob/d4097719d49437c87ad8ca7136d74791fe762dcf/openmmtools/mcmc.py#L1680
class TorsionNeglectingMove(RotateDisplaceMove):
    """A metropolized move that randomly displaces & rotates a subset of atoms.
    Also, it calculates the energies with the torsional potential evaluated lazily
    only sometimes after a few timesteps.
    Parameters
    ----------
    timesteps : int
        number of timesteps to do before updating the chached torsional component 
        of the potential energy
    displacement_sigma : simtk.unit.Quantity
        The standard deviation of the normal distribution used to propose the
        random displacement (units of length, default is 1.0*nanometer).
    atom_subset : slice or list of int, optional
        If specified, the move is applied only to those atoms specified by these
        indices. If None, the move is applied to all atoms (default is None).
    context_cache : openmmtools.cache.ContextCache, optional
        The ContextCache to use for Context creation. If None, the global cache
        openmmtools.cache.global_context_cache is used (default is None).
    Attributes
    ----------
    n_accepted : int
        The number of proposals accepted.
    n_proposed : int
        The total number of attempted moves.
    cached_torsional_term : float
        The last value of the torsional component of the potential energy
    timestep : int
        The number of moves taken
    displacement_sigma
    timesteps
    atom_subset
    context_cache
    See Also
    --------
    TotateDisplaceMove
    """

    def __init__(self, timesteps=5, **kwargs):
        super(TorsionNeglectingMove, self).__init__(**kwargs)
        self.timestep = 0
        self.timesteps = timesteps
        self.cached_torsional_term = 0

    def __getstate__(self):
        serialization = super(TorsionNeglectingMove, self).__getstate__()
        serialization['timesteps'] = self.timesteps
        serialization['timestep'] = self.timestep
        return serialization

    def __setstate__(self, serialization):
        super(TorsionNeglectingMove, self).__setstate__(serialization)
        self.timesteps = serialization['timesteps']
        self.timestep = serialization['timestep']

    def _calculate_potential_energy(self, thermodynamic_state, context):
        self.timestep += 1
        # TODO: calculate potential energy / replace following line
        return thermodynamic_state.reduced_potential(context)
        # if (timestep % timesteps == 0): calculate new torsional term
        #   self.cached_torsional_term = 
        # 
        # anyway calculate all other terms
        # sum them up (+ self.cached_torsional_term) and return
