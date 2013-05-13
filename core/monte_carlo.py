"""Monte Carlo Classes and Functions."""

import math


def concentric_sample_disk(u1, u2):
    """Generate sample on a disk."""
    # Map uniform random numbers to $[-1,1]^2$
    sx = 2 * u1 - 1
    sy = 2 * u2 - 1

    # Map square to $(r,\theta)$

    # Handle degeneracy at the origin
    if (sx == 0.0 and sy == 0.0):
        return 0.0, 0.0

    if (sx >= -sy):
        if (sx > sy):
            # Handle first region of disk
            r = sx
            if (sy > 0.0):
                theta = sy/r
            else:
                theta = 8.0 + sy/r
        else:
            # Handle second region of disk
            r = sy
            theta = 2.0 - sx/r;
    else:
        if (sx <= sy):
            # Handle third region of disk
            r = -sx
            theta = 4.0 - sy/r
        else:
            # Handle fourth region of disk
            r = -sy
            theta = 6.0 + sx/r
    theta *= math.pi / 4.0
    dx = r * math.cos(theta)
    dy = r * math.sin(theta)
    return dx, dy
