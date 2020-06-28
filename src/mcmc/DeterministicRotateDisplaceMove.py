import numpy as np
from numpy.random import default_rng

from simtk import openmm, unit

from .DeterministicMetropolizedMove import DeterministicMetropolizedMove


# inspired by https://github.com/choderalab/openmmtools/blob/d4097719d49437c87ad8ca7136d74791fe762dcf/openmmtools/mcmc.py#L1757
# and https://github.com/choderalab/openmmtools/blob/d4097719d49437c87ad8ca7136d74791fe762dcf/openmmtools/mcmc.py#L1680
class DeterministicRotateDisplaceMove(DeterministicMetropolizedMove):
    """A metropolized move that deterministic-randomly displaces & rotates a subset of atoms.
    Parameters
    ----------
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
    displacement_sigma
    atom_subset
    context_cache
    See Also
    --------
    MetropolizedMove
    """

    # still to be overriden, since super().__init__ expects class name too
    def __init__(self, displacement_sigma=1.0*unit.nanometer, **kwargs):
        super(DeterministicRotateDisplaceMove, self).__init__(**kwargs)
        self.displacement_sigma = displacement_sigma
        self.rotate_rng = default_rng(5)

    def displace_positions(self, positions, displacement_sigma=1.0*unit.nanometer):
        """Return the positions after applying a random displacement to them.
        Parameters
        ----------
        positions : nx3 numpy.ndarray simtk.unit.Quantity
            The positions to displace.
        displacement_sigma : simtk.unit.Quantity
            The standard deviation of the normal distribution used to propose
            the random displacement (units of length, default is 1.0*nanometer).
        Returns
        -------
        rotated_positions : nx3 numpy.ndarray simtk.unit.Quantity
            The displaced positions.
        """
        positions_unit = positions.unit
        unitless_displacement_sigma = displacement_sigma / positions_unit
        displacement_vector = unit.Quantity(self.rotate_rng.standard_normal(3) * unitless_displacement_sigma,
                                            positions_unit)
        # print("Displacement vector: {}".format(displacement_vector))
        return positions + displacement_vector

    def rotate_positions(self, positions):
        """Return the positions after applying a random rotation to them.
        Parameters
        ----------
        positions : nx3 numpy.ndarray simtk.unit.Quantity
            The positions to rotate.
        Returns
        -------
        rotated_positions : nx3 numpy.ndarray simtk.unit.Quantity
            The rotated positions.
        """
        positions_unit = positions.unit
        x_initial = positions / positions_unit

        # Compute center of geometry of atoms to rotate.
        x_initial_mean = x_initial.mean(0)

        # Generate a random rotation matrix.
        rotation_matrix = self.generate_random_rotation_matrix()

        # Apply rotation.
        x_proposed = (rotation_matrix * np.matrix(x_initial -
                                                  x_initial_mean).T).T + x_initial_mean
        return unit.Quantity(x_proposed, positions_unit)

    def generate_random_rotation_matrix(self):
        """Return a random 3x3 rotation matrix.
        Returns
        -------
        Rq : 3x3 numpy.ndarray
            The random rotation matrix.
        """
        q = self._generate_uniform_quaternion()
        return self._rotation_matrix_from_quaternion(q)

    def _rotation_matrix_from_quaternion(self, q):
        """Compute a 3x3 rotation matrix from a given quaternion (4-vector).
        Parameters
        ----------
        q : 1x4 numpy.ndarray
            Quaterion (need not be normalized, zero norm OK).
        Returns
        -------
        Rq : 3x3 numpy.ndarray
            Orthogonal rotation matrix corresponding to quaternion q.
        Examples
        --------
        >>> q = np.array([0.1, 0.2, 0.3, -0.4])
        >>> Rq = MCRotationMove._rotation_matrix_from_quaternion(q)
        References
        ----------
        [1] http://en.wikipedia.org/wiki/Rotation_matrix#Quaternion
        """

        w, x, y, z = q
        Nq = (q**2).sum()  # Squared norm.
        if Nq > 0.0:
            s = 2.0 / Nq
        else:
            s = 0.0

        X = x*s
        Y = y*s
        Z = z*s
        wX = w*X
        wY = w*Y
        wZ = w*Z
        xX = x*X
        xY = x*Y
        xZ = x*Z
        yY = y*Y
        yZ = y*Z
        zZ = z*Z

        Rq = np.matrix([[1.0-(yY+zZ),     xY-wZ,          xZ+wY],
                        [xY+wZ,        1.0-(xX+zZ),       yZ-wX],
                        [xZ-wY,           yZ+wX,    1.0-(xX+yY)]])

        return Rq

    def _generate_uniform_quaternion(self):
        """Generate a uniform normalized quaternion 4-vector.
        References
        ----------
        [1] K. Shoemake. Uniform random rotations. In D. Kirk, editor,
        Graphics Gems III, pages 124-132. Academic, New York, 1992.
        [2] Described briefly here: http://planning.cs.uiuc.edu/node198.html
        Examples
        --------
        >>> q = MCRotationMove._generate_uniform_quaternion()
        """
        u = self.rotate_rng.random(3)
        q = np.array([np.sqrt(1-u[0])*np.sin(2*np.pi*u[1]),
                      np.sqrt(1-u[0])*np.cos(2*np.pi*u[1]),
                      np.sqrt(u[0])*np.sin(2*np.pi*u[2]),
                      np.sqrt(u[0])*np.cos(2*np.pi*u[2])])
        return q

    def __getstate__(self):
        serialization = super(
            DeterministicRotateDisplaceMove, self).__getstate__()
        serialization['displacement_sigma'] = self.displacement_sigma
        return serialization

    def __setstate__(self, serialization):
        super(DeterministicRotateDisplaceMove,
              self).__setstate__(serialization)
        self.displacement_sigma = serialization['displacement_sigma']

    def _propose_positions(self, initial_positions):
        """Implement MetropolizedMove._propose_positions for apply()."""
        # print("Proposing for {} positions".format(len(initial_positions)))
        displaced_positions = self.displace_positions(
            initial_positions, self.displacement_sigma)
        return displaced_positions
        # rotated_positions = self.rotate_positions(initial_positions)
        # return rotated_positions
