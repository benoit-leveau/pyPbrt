"""Interface Classes for Integrators."""


from abc import ABCMeta, abstractmethod

from core.geometry import abs_dot
from core.reflection import BSDFSample


class Integrator(object):
    
    """Interface for Integrator Classes."""
    
    def preprocess(self, scene, camera, renderer):
        """Optional preprocess for the integrator."""
        pass

    def request_samples(self, sampler, sample, scene):
        """Request samples to the sampler."""
        pass


class SurfaceIntegrator(Integrator):

    """SurfaceIntegrator Interface."""

    __metaclass__ = ABCMeta

    @abstractmethod
    def Li(self, scene, renderer, ray, intersection, sample, rng):
        """Return light emitted at the intersection."""
        pass


def specular_reflect(ray, bsdf, rng, intersection, renderer, scene, sample):
    """Trace a ray for specular reflection."""
    wo = -ray.d
    p = bsdf.dg_shading.p
    n = bsdf.dg_shading.nn

    # retrieve the value, vector and probability of a sampled direction
    f, wi, pdf = bsdf.sample_f(wo, BSDFSample.from_rng(rng),
                               BxDFType(BSDF_REFLECTION | BSDF_SPECULAR))

    L = 0.0
    if pdf>0.0 and (not f.is_black()) and abs_dot(wi, n) != 0.0:
        # compute ray differential rd for specular reflection
        Li = renderer.Li(scene, rd, sample, rng)
        L = f * Li * abs_dot(wi, n) / pdf
    return L


def specular_transmit(ray, bsdf, rng, intersection, renderer, scene, sample):
    """Trace a ray for specular transmission."""
    wo = -ray.d
    p = bsdf.dg_shading.p
    n = bsdf.dg_shading.nn

    # retrieve the value, vector and probability of a sampled direction
    f, wi, pdf = bsdf.sample_f(wo, BSDFSample.from_rng(rng),
                               BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR))

    L = 0.0
    if pdf>0.0 and (not f.is_black()) and abs_dot(wi, n) != 0.0:
        # compute ray differential rd for specular reflection
        Li = renderer.Li(scene, rd, sample, rng)
        L = f * Li * abs_dot(wi, n) / pdf
    return L
