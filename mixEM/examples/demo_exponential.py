#!/usr/bin/env python

import numpy as np
import mixem
from mixem.distribution import ExponentialDistribution


def generate_data():
    dist_params = [1, 10]
    weights = [0.4, 0.6]

    n_data = 10000
    data = np.zeros((n_data))
    for i in range(n_data):
        dpi = np.random.choice(range(len(dist_params)), p=weights)
        dp = dist_params[dpi]
        data[i] = np.random.exponential(scale=1.0 / dp)

    return data


def recover(data):

    init_params = np.random.choice(data, size=2)

    weight, distributions, ll = mixem.em(data, [ExponentialDistribution(l) for l in init_params])

    print(weight, distributions, ll)


if __name__ == '__main__':
    data = generate_data()
    recover(data)
