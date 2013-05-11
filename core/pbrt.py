"""Helper Functions."""


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
