"""Random Number Generator Class."""


class RNG(object):

    """Pseudo-Random Number Generator Class."""

    N = 624 # size of the state vector
    
    def __init__(self, s=5489):
        """Default constructor for RNG."""

        # state vector
        self.mt = [0]*RNG.N
        self.mti = RNG.N+1

        # init the generator
        self.seed(s)

    def seed(self, s):
        """Initialize the prng with a seed."""
        pass

    def random_float(self):
        """Generate a random float."""
        pass

    def random_uint(self):
        """Generate a random unsigned integer."""
        pass
    
