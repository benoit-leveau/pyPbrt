"""Light Classes."""


class Light(object):

    """Base Light Class."""

    def __init__(self):
        pass


class LightSample(object):

    """LightSample Class."""

    def __init__(self):
        pass

    @classmethod
    def from_rng(cls, rng):
        ret = cls()
        ret.u_pos[0] = rng.random_float()
        ret.u_pos[1] = rng.random_float()
        ret.u_component = rng.random_float()
