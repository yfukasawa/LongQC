# coding=utf-8
import numpy as np

from .distribution import Distribution


class ExponentialDistribution(Distribution):
    """Exponential distribution with parameter (lambda)."""

    def __init__(self, lmbda):
        self.lmbda = lmbda

    def log_density(self, data):
        assert(len(data.shape) == 1), "Expect 1D data!"
        return np.log(self.lmbda) - self.lmbda * data

    def estimate_parameters(self, data, weights):
        assert(len(data.shape) == 1), "Expect 1D data!"
        self.lmbda = np.sum(weights) / np.sum(weights * data)

    def __repr__(self):
        return "Exp[λ={lmbda:.4g}]".format(lmbda=self.lmbda)
