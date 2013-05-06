"""Interface Classes for Integrators."""


from core.geometry import abs_dot
from core.reflection import BSDFSample


class Integrator(object):

    """Interface for Integrator Classes."""

    def preprocess(self):
        pass

    def request_samples(self):
        pass


class SurfaceIntegrator(Integrator):

    """SurfaceIntegrator Interface."""

    def Li(self, scene, renderer, ray, isect, sample, rng):
        pass


def specular_reflect(ray, bsdf, rng, isect, renderer, scene, sample):
    """Trace a ray for specular reflection."""
    wo = -ray.d * wi ###
    p = bsdf.dg_shading.p
    n = bsdf.dg_shading.nn

    # retrieve the value, vector and probability of a sampled direction
    f, wi, pdf = bsdf.sample_f(wo, BSDFSample(rng),
                               BxDFType(BSDF_REFLECTION | BSDF_SPECULAR))

    L = 0.0
    if pdf>0.0 and (not f.is_black()) and abs_dot(wi, n) != 0.0:
        # compute ray differential rd for specular reflection
        Li = renderer.Li(scene, rd, sample, rng)
        L = f * Li * abs_dot(wi, n) / pdf
    return L


def specular_transmit(ray, bsdf, rng, isect, renderer, scene, sample):
    """Trace a ray for specular transmission."""
    wo = -ray.d * wi ###
    p = bsdf.dg_shading.p
    n = bsdf.dg_shading.nn

    # retrieve the value, vector and probability of a sampled direction
    f, wi, pdf = bsdf.sample_f(wo, BSDFSample(rng),
                               BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR))

    L = 0.0
    if pdf>0.0 and (not f.is_black()) and abs_dot(wi, n) != 0.0:
        # compute ray differential rd for specular reflection
        Li = renderer.Li(scene, rd, sample, rng)
        L = f * Li * abs_dot(wi, n) / pdf
    return L
