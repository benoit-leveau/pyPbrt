"""Material Class."""

from abc import ABCMeta, abstractmethod


class Material(object):

    """Interface Class for all Materials."""

    def __init__(self):
        """Default constructor for Material."""
        pass

    @abstractmethod
    def get_bsdf(self, dg_geom, dg_shading):
        """Get the BSDF at intersection."""
        pass

    def get_bssrdf(self, dg_geom, dg_shading):
        """Get the BSDF at intersection."""
        return None

    def bump(self, d, dg_geom, dg_shading):
        """Retrieve bumped intersection on material."""
        raise NotImplementedError()

    @abstractmethod
    def __str__(self):
        """Return a string describing the material."""
        pass
