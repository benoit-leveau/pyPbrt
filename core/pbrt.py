"""Helper Functions."""

import math


def round_up_pow_2(v):
    """Round up specified integer to the closest (upper) power of 2."""
    v -= 1
    v |= v >> 1
    v |= v >> 2
    v |= v >> 4
    v |= v >> 8
    v |= v >> 16
    return v+1


def lerp(t, v1, v2):
    """Linear Interpolation between v1 and v2."""
    return (1.0-t)*v1 + t*v2


def eq(a, b, epsilon=1e-10):
    """Compares two float values."""
    return a==b or abs(a-b)<epsilon


def clamp(value, low, high):
    """Return the clamped value in the specified range."""
    if value < low:
        return low
    elif value > high:
        return high
    return value


def quadratic(A, B, C):
    # Find quadratic discriminant
    discrim = B*B - 4.0*A*C
    if (discrim < 0.0):
        return False, float('inf'), float('inf')
    root_discrim = math.sqrt(discrim)

    # Compute quadratic _t_ values
    if (B < 0):
        q = -0.5 * (B-root_discrim)
    else:
        q = -0.5 * (B+root_discrim)

    if A != 0.0:
        t0 = q / A
    else:
        t0 = float('inf')
    if q != 0.0:
        t1 = C / q
    else:
        t1 = float('inf')    
    if (t0 > t1):
        t0, t1 = t1, t0
    return True, t0, t1

def round_to_int(value):
    """Round to closest integer."""
    # might be problematic, see online discussions about round errors
    return int(value+0.5)

