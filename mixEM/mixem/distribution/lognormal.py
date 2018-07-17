# coding=utf-8
import numpy as np
import scipy.stats

from .distribution import Distribution


class LogNormalDistribution(Distribution):
    """Univariate log-normal distribution with parameters (mu, sigma)."""

    def __init__(self, mu, sigma):
        self.mu = mu
        self.sigma = sigma

    def log_density(self, data):
        assert(len(data.shape) == 1), "Expect 1D data!"

        return - (np.log(data) - self.mu) ** 2 / (2 * self.sigma ** 2) - np.log(self.sigma) - 0.5 * np.log(2 * np.pi) - np.log(data)

    def estimate_parameters(self, data, weights):
        assert(len(data.shape) == 1), "Expect 1D data!"

        wsum = np.sum(weights)

        self.mu = np.sum(weights * np.log(data)) / wsum
        self.sigma = np.sqrt(np.sum(weights * (np.log(data) - self.mu) ** 2) / wsum)

    def __repr__(self):
        return "LogNormal[μ={mu:.4g}, σ={sigma:.4g}]".format(mu=self.mu, sigma=self.sigma)


    def get_mu(self):
        return self.mu

    def get_sigma(self):
        return self.sigma
