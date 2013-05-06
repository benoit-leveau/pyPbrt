"""Functions and Classes for Parallelization of pyPbrt."""


import multiprocessing


def num_system_cores():
    # if pbrt_options.n_cores > 0:
    #    return pbrt_options.n_cores
    return multiprocessing.cpu_count()
