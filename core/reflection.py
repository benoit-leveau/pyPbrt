"""BSDF Classes."""


class BSDFSample(object):

    """BSDFSample class."""

    def __init__(self, up0, up1, ucomp):
        """Default Constructor for BSDFSample."""
        self.u_dir = [up0, up1]
        self.u_component = ucomp

    @classmethod
    def from_rng(cls, rng):
        """Constructs a random sample."""
        return cls(rng.random_float(),
                   rng.random_float(),
                   rng.random_float())

    def __str__(self):
        """Return a string describing the sample."""
        return "BSDFSample (u_dir=[%f, %f], u_component=%f)" % (self.u_dir[0],
                                                                self.u_dir[1],
                                                                self.u_component)
        
                   
