"""Class for the Uber Material."""

from core.material import Material


def UberMaterial(Material):

    """UberMaterial Class."""

    def __init__(self, kd, ks, kr, kt, rough, opacity, eta, bump):
        """Default constructor for UberMaterial."""
        self.Kd = kd
        self.Ks = ks
        self.Kr = kr
        self.Kt = kt
        self.roughness = rough
        self.opacity = opacity
        self.eta = eta
        self.bump_map = bump

    def get_bsdf(self, dg_geom, dg_shading):
        """Get the BSDF at intersection."""
        if self.bump_map:
            dgs = self.bump(self.bump_map, dg_geom, dg_shading)
        else:
            dgs = dg_shading

        bsdf = BSDF(dgs, dg_geom.nn)
        op = self.opacity.evaluate(dgs).clamp()
        if op != Spectrum(1.0):
            tr = SpecularTransmission(-op + Spectrum(1.0), 1.0, 1.0)
            bsdf.add(tr)

        kd = op * self.Kd.evaluate(dgs).clamp()
        if not kd.is_black():
            diff = Lambertian(kd)
            bsdf.add(diff)

        e = self.eta.evaluate(dgs)
        ks = op * self.Ks.evaluate(dgs).clamp()
        if not ks.is_black():
            fresnel = FresnelDielectric(e, 1.0)
            rough = self.roughness.evaluate(dgs)
            spec = Microfacet(ks, fresnel, Blinn(1.0/rough))
            bsdf.add(spec)

        kr = op * self.Kr.evaluate(dgs).clamp()
        if not kr.is_black():
            fresnel = FresnelDielectric(e, 1.0)
            bsdf.add(SpecularReflection(kr, fresnel))

        kt = op * self.Kt.evaluate(dgs).clamp()
        if not kt.is_black():
            bsdf.add(SpecularTransmission(kt, e, 1.0))

        return bsdf

    def get_bssrdf(self, dg_geom, dg_shading):
        """Get the BSDF at intersection."""
        pass

    def __str__(self):
        """Return a string describing the material."""
        return "UberMaterial ( )"

        
