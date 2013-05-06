"""Classes for Spectrum/Color."""


import math


class CoefficientSpectrum(object):

    """Spectrum Class."""

    def __init__(self, n_samples, v=0.0):
        """Default constructor for CoefficientSpectrum."""
        self.c = [v] * n_samples

    def has_NaNs(self):
        """Return True if any of the coefficient is NaN."""
        for c in self.c:
            if math.isnan(c):
                return True
        return False

    def __str__(self):
        """Return a string describing the spectrum."""
        return "CoefficientSpectrum (%s)" % ", ".join([str(c) for c in self.c])


class RGBSpectrum(CoefficientSpectrum):

    """Specialization of Spectrum Class for 3 values (R, G, B)."""

    def __init__(self, v=0.0):
        """Default constructor for RGBSpectrum."""
        super(RGBSpectrum, self).__init__(n_samples=3, v=v)

    def __str__(self):
        """Return a string describing the rgb color."""
        return "RGBSpectrum (r=%f g=%f b=%f)" % (self.c[0],
                                                 self.c[1],
                                                 self.c[2])
