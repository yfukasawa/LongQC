#!/usr/bin/env python

import numpy as np
import mixem
from mixem.distribution import GeometricDistribution


def generate_data():
    dist_params = [0.1, 0.3]
    weights = [0.4, 0.6]

    n_data = 10000
    data = np.zeros((n_data))
    for i in range(n_data):
        dpi = np.random.choice(range(len(dist_params)), p=weights)
        dp = dist_params[dpi]
        data[i] = np.random.geometric(p=dp)

    return data


def recover(data):

    weight, distributions, ll = mixem.em(data, [
        GeometricDistribution(0.8),
        GeometricDistribution(0.1),
    ])

    print(weight, distributions, ll)


if __name__ == '__main__':
    data = generate_data()
    recover(data)
