"""WhittedIntegrator Class."""


from core.spectrum import Spectrum
from core.light import LightSample
from core.geometry import abs_dot
from core.integrator import SurfaceIntegrator


class WhittedIntegrator(SurfaceIntegrator):

    """WhittedIntegrator."""

    def __init__(self, max_depth):
        """Default constructor for WhittedIntegrator."""
        self.max_depth = max_depth

    def Li(self, scene, renderer, ray, intersection, sample, rng):
        """Computes the radiance along a ray."""
        return Spectrum(1.0)
    
        L = Spectrum(0.0)

        # evaluate BSDF at hit point
        bsdf = intersection.get_bsdf(ray)
        
        # initialize common variables for Whitted integrator
        p = bsdf.dg_shading.p
        n = bsdf.dg_shading.nn
        wo = -ray.d

        # compute emitted light if ray hit an area light source
        L += intersection.Le(wo)
        
        # add contribution of each light source
        for light in self.scene.lights:
            light_sample = LightSample.from_rng(rng)
            Li, wi, pdf, visibility = light.sample_L(p,
                                                     intersection.ray_epsilon,
                                                     light_sample,
                                                     ray.time)
            if Li.is_black() or pdf == 0.0:
                continue

            f = bsdf.f(wo, wi)
            if (not f.is_black()) and visibility.unoccluded(scene):
                L+= f * Li * abs_dot(wi, n) * visibility.transmittance(scene,
                                                                       renderer,
                                                                       sample,
                                                                       rng) / pdf
        if ray.depth+1 < self.max_depth:
            # trace rays for specular reflection and refraction
            L += specular_reflect(ray, bsdf, rng, intersection,
                                       renderer, scene, sample)

            L+= specular_transmit(ray, bsdf, rng, intersection,
                                       renderer, scene, sample)

        return L

    def __str__(self):
        """Return a string describing the integrator."""
        return "WhittedIntegrator (max_depth='%d')" % self.max_depth
