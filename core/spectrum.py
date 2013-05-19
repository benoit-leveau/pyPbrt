"""Classes for Spectrum/Color."""


import math


class CoefficientSpectrum(object):

    """Spectrum Class."""

    def __init__(self, n_samples, v=0.0):
        """Default constructor for CoefficientSpectrum."""
        self.n_samples = n_samples
        self.c = [v] * n_samples

    @classmethod
    def from_spectrum(cls, s2):
        """Copy constructor."""
        ret = cls(s2.n_samples)
        for i in xrange(s2.n_samples):
            ret.c[i] = s2.c[i]
        return ret
        
    def has_nan(self):
        """Return True if any of the coefficient is NaN."""
        for c in self.c:
            if math.isnan(c):
                return True
        return False

    def __add__(self, s2):
        """Overload the addition operator."""
        ret = type(self).from_spectrum(self)
        for i in xrange(self.n_samples):
            ret.c[i] += s2.c[i]
        return ret

    def __sub__(self, s2):
        """Overload the subtraction operator."""
        ret = type(self).from_spectrum(self)
        for i in xrange(self.n_samples):
            ret.c[i] -= s2.c[i]
        return ret

    def __mul__(self, s2):
        """Overload the multiplication operator."""
        ret = type(self).from_spectrum(self)
        if isinstance(s2, CoefficientSpectrum):
            for i in xrange(self.n_samples):
                ret.c[i] *= s2.c[i]
        else:
            for i in xrange(self.n_samples):
                ret.c[i] *= s2            
        return ret

    def __rmul__(self, s2):
        """Overload the right multiplication operator."""
        ret = type(self).from_spectrum(self)
        if isinstance(s2, CoefficientSpectrum):
            for i in xrange(self.n_samples):
                ret.c[i] *= s2.c[i]
        else:
            for i in xrange(self.n_samples):
                ret.c[i] *= s2            
        return ret

    def __div__(self, s2):
        """Overload the division operator."""
        ret = type(self).from_spectrum(self)
        if isinstance(s2, CoefficientSpectrum):
            for i in xrange(self.n_samples):
                ret.c[i] /= s2.c[i]
        else:
            for i in xrange(self.n_samples):
                ret.c[i] /= s2            
        return ret

    def __str__(self):
        """Return a string describing the spectrum."""
        return "CoefficientSpectrum (%s)" % ", ".join([str(c) for c in self.c])


class RGBSpectrum(CoefficientSpectrum):

    """Specialization of Spectrum Class for 3 values (R, G, B)."""
    
    y_weight = [0.212671, 0.715160, 0.072169]

    def __init__(self, v=0.0):
        """Default constructor for RGBSpectrum."""
        super(RGBSpectrum, self).__init__(n_samples=3, v=v)

    def y(self):
        """Return luminance."""
        return self.c[0] * RGBSpectrum.y_weight[0] + \
               self.c[1] * RGBSpectrum.y_weight[1] + \
               self.c[2] * RGBSpectrum.y_weight[2]

    def to_xyz(self):
        """Convert to XYZ."""
        return rgb_to_xyz(self.c[0], self.c[1], self.c[2])

    def to_rgb(self):
        """Convert to RGB."""
        return self.c[0], self.c[1], self.c[2]
    
    def __str__(self):
        """Return a string describing the rgb color."""
        return "RGBSpectrum (r=%f g=%f b=%f)" % (self.c[0],
                                                 self.c[1],
                                                 self.c[2])


# default Spectrum is RGBSpectrum
Spectrum = RGBSpectrum


def rgb_to_xyz(r, g, b):
    """Convert from RGB to XYZ."""
    x = 0.412453*r + 0.357580*g + 0.180423*b
    y = 0.212671*r + 0.715160*g + 0.072169*b
    z = 0.019334*r + 0.119193*g + 0.950227*b
    return x, y, z


def xyz_to_rgb(x, y, z):
    """Convert from XYZ to RGB."""
    r =  3.240479*x - 1.537150*y - 0.498535*z
    g = -0.969256*x + 1.875991*y + 0.041556*z
    b =  0.055648*x - 0.204043*y + 1.057311*z
    return r, g, b
