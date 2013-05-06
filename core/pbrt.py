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
