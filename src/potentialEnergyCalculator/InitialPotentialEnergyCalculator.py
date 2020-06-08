class InitialPotentialEnergyCalculator(object):
    """A calculator of the energies.
    """

    def __init__(self, **kwargs):
        super(InitialPotentialEnergyCalculator, self).__init__(**kwargs)

    def __getstate__(self):
        serialization = super(InitialPotentialEnergyCalculator, self).__getstate__()
        return serialization

    def __setstate__(self, serialization):
        super(InitialPotentialEnergyCalculator, self).__setstate__(serialization)

    def _calculate_potential_energy(self, thermodynamic_state, context):
        # calculate potential energy
        return thermodynamic_state.reduced_potential(context)
