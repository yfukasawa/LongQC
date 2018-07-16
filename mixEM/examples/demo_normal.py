#!/usr/bin/env python

import numpy as np
import mixem
from mixem.distribution import NormalDistribution


def generate_data():
    dist_params = [
        (4, 1),
        (1, 0.5)
    ]

    weights = [0.3, 0.7]

    n_data = 5000
    data = np.zeros((n_data,))
    for i in range(n_data):
        dpi = np.random.choice(range(len(dist_params)), p=weights)
        mu, sigma = dist_params[dpi]
        data[i] = np.random.normal(loc=mu, scale=sigma)

    return data


def recover(data):

    mu = np.mean(data)
    sigma = np.var(data)

    init_params = [
        (mu + 0.1, sigma),
        (mu - 0.1, sigma)
    ]

    weight, distributions, ll = mixem.em(data, [NormalDistribution(mu, sigma) for mu_sigma in init_params])

    print(weight, distributions, ll)


if __name__ == '__main__':
    data = generate_data()
    recover(data)
