from simtk import unit


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
        serialization = super(
            TorsionNeglectingPotentialEnergyCalculator, self).__getstate__()
        serialization['timesteps'] = self.timesteps
        serialization['timestep'] = self.timestep
        return serialization

    def __setstate__(self, serialization):
        super(TorsionNeglectingPotentialEnergyCalculator,
              self).__setstate__(serialization)
        self.timesteps = serialization['timesteps']
        self.timestep = serialization['timestep']

    def _calculate_potential_energy(self, thermodynamic_state, context):
        sys = context.getSystem()

        # this might should be done in the init
        # if we were still inheriting from metropolizedMove
        # but we are not and don't care about performance anyway (Python, you know...), so
        for i in range(0, sys.getNumForces()):
            force = sys.getForce(i)
            className = force.__class__.__name__

            # print("Got force: {} with group: {}".format(
            #     className, force.getForceGroup()))

            if ("TorsionForce" in className):
                force.setForceGroup(2)
            else:
                force.setForceGroup(1)
        
        # see also: https://openmmtools.readthedocs.io/en/latest/_modules/openmmtools/states.html#ThermodynamicState.reduced_potential
        # calculate all forces but torsional
        n_particles = context.getSystem().getNumParticles()
        openmm_state = context.getState(getEnergy=True, groups=1)
        potential_energy = openmm_state.getPotentialEnergy()
        volume = openmm_state.getPeriodicBoxVolume()

        # then, calculate also the torsional energy
        calculateTorsionalForce = self.timestep % self.timesteps == 0
        if (calculateTorsionalForce):
            openmm_state = context.getState(
                getEnergy=True, groups=2)
            self.cached_torsional_term = openmm_state.getPotentialEnergy()

        # print("V: {}, T: {}".format(potential_energy, self.cached_torsional_term))
        potential_energy += self.cached_torsional_term
        self.timestep += 1

        # Convert potential energy into reduced potential.
        # see also: https://openmmtools.readthedocs.io/en/latest/_modules/openmmtools/states.html#ThermodynamicState._compute_reduced_potential
        beta = 1.0 / (unit.BOLTZMANN_CONSTANT_kB *
                      thermodynamic_state.temperature)
        reduced_potential = potential_energy / unit.AVOGADRO_CONSTANT_NA
        if thermodynamic_state.pressure is not None:
            reduced_potential += thermodynamic_state.pressure * volume
        return beta * reduced_potential
