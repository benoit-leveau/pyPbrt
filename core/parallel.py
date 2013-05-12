"""Functions and Classes for Parallelization of pyPbrt."""


import multiprocessing
from abc import ABCMeta, abstractmethod


class Task(object):

    """Interface class for tasks."""

    __metaclass__ = ABCMeta

    @abstractmethod
    def run(self):
        """Run the task."""
        pass


def num_system_cores():
    # if pbrt_options.n_cores > 0:
    #    return pbrt_options.n_cores
    return multiprocessing.cpu_count()
