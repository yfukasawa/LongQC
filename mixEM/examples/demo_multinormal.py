#!/usr/bin/env python

import numpy as np
import mixem
from mixem.distribution import MultivariateNormalDistribution


def generate_data():
    dist_params = [
        (np.array([4]), np.diag([1])),
        (np.array([1]), np.diag([0.5]))
    ]

    weights = [0.3, 0.7]

    n_data = 5000
    data = np.zeros((n_data, 1))
    for i in range(n_data):
        dpi = np.random.choice(range(len(dist_params)), p=weights)
        dp = dist_params[dpi]
        data[i] = np.random.multivariate_normal(dp[0], dp[1])

    return data


def recover(data):

    mu = np.mean(data)
    sigma = np.var(data)

    init_params = [
        (np.array([mu + 0.1]), np.diag([sigma])),
        (np.array([mu - 0.1]), np.diag([sigma]))
    ]

    weight, distributions, ll = mixem.em(data, [MultivariateNormalDistribution(mu, sigma) for mu, sigma in init_params])

    print(weight, distributions, ll)


if __name__ == '__main__':
    data = generate_data()
    recover(data)
