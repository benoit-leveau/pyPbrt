"""Renderer Class."""

from abc import ABCMeta, abstractmethod


class Renderer(metaclass=ABCMeta):

    """Renderer Interface."""

    def __init__(self):
        """Default constructor."""
        pass

    @abstractmethod
    def render(self, scene):
        """Render a scene."""
        pass

    @abstractmethod
    def Li(self, scene, ray, sample, rng, isect=None, T=None):
        """Compute the incident radiance along a given ray."""
        pass

    @abstractmethod
    def transmittance(self, scene, ray, sample, rng):
        """Compute the attenuation by volumetric scattering along a ray."""
        pass
