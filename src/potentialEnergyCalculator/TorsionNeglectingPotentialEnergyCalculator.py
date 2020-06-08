class TorsionNeglectingPotentialEnergyCalculator(object):
    """A calculator of the energies with the torsional potential evaluated lazily
    only sometimes after a few timesteps.

    Parameters
    ----------
    timesteps : int
        The number of timesteps to do before updating the chached torsional component 
        of the potential energy

    Attributes
    ----------
    cached_torsional_term : float
        The last value of the torsional component of the potential energy
    timestep : int
        The number of moves taken
    timesteps : int 
        The number of timesteps to do before updating the chached torsional component 
        of the potential energy 
    """

    def __init__(self, timesteps=5, **kwargs):
        super(TorsionNeglectingPotentialEnergyCalculator, self).__init__(**kwargs)
        self.timestep = 0
        self.timesteps = timesteps
        self.cached_torsional_term = 0

    def __getstate__(self):
        serialization = super(TorsionNeglectingPotentialEnergyCalculator, self).__getstate__()
        serialization['timesteps'] = self.timesteps
        serialization['timestep'] = self.timestep
        return serialization

    def __setstate__(self, serialization):
        super(TorsionNeglectingPotentialEnergyCalculator, self).__setstate__(serialization)
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
