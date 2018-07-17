# coding=utf-8
import numpy as np
import scipy.stats

from .distribution import Distribution


class NormalDistribution(Distribution):
    """Univariate normal distribution with parameters (mu, sigma)."""

    def __init__(self, mu, sigma):
        self.mu = mu
        self.sigma = sigma

    def log_density(self, data):
        assert(len(data.shape) == 1), "Expect 1D data!"

        return - (data - self.mu) ** 2 / (2 * self.sigma ** 2) - np.log(self.sigma) - 0.5 * np.log(2 * np.pi)

    def estimate_parameters(self, data, weights):
        assert(len(data.shape) == 1), "Expect 1D data!"

        wsum = np.sum(weights)

        self.mu = np.sum(weights * data) / wsum
        self.sigma = np.sqrt(np.sum(weights * (data - self.mu) ** 2) / wsum)

    def __repr__(self):
        return "Norm[μ={mu:.4g}, σ={sigma:.4g}]".format(mu=self.mu, sigma=self.sigma)

    def get_mu(self):
        return self.mu

    def get_sigma(self):
        return self.sigma

class MultivariateNormalDistribution(Distribution):
    """Multivariate normal distribution with parameters (mu, Sigma)."""

    def __init__(self, mu, sigma):
        mu = np.array(mu)
        sigma = np.array(sigma)

        assert len(mu.shape) == 1, "Expect mu to be 1D vector!"
        assert len(sigma.shape) == 2, "Expect sigma to be 2D matrix!"

        assert sigma.shape[0] == sigma.shape[1], "Expect sigma to be a square matrix!"

        self.mu = mu
        self.sigma = sigma

    def log_density(self, data):
        return scipy.stats.multivariate_normal.logpdf(data, self.mu, self.sigma)

    def estimate_parameters(self, data, weights):
        self.mu = np.sum(data * weights[:, np.newaxis], axis=0) / np.sum(weights)

        center_x = data - self.mu[np.newaxis, :]

        # sigma = (np.diag(weights) @ center_x).T @ center_x / np.sum(weights)
        self.sigma = np.dot(
            np.dot(
                np.diag(weights),
                center_x
            ).T,
            center_x
        ) / np.sum(weights)

    def __repr__(self):
        po = np.get_printoptions()

        np.set_printoptions(precision=3)

        try:
            result = "MultiNorm[μ={mu}, σ={sigma}]".format(mu=self.mu, sigma=str(self.sigma).replace("\n", ","))
        finally:
            np.set_printoptions(**po)

        return result

    def get_mu(self):
        return self.mu

    def get_sigma(self):
        return self.sigma
