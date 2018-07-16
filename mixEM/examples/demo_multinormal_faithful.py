#!/usr/bin/env python

import os

import numpy as np
import pandas as pd

import mixem
from mixem.distribution import MultivariateNormalDistribution


def main():

    data = pd.read_csv(os.path.join(os.path.dirname(__file__), "faithful.csv"))

    data = np.array(data)

    init_params = [
        (np.array((2, 50)), np.identity(2)),
        (np.array((4, 80)), np.identity(2)),
    ]

    weight, distributions, ll = mixem.em(data, [MultivariateNormalDistribution(mu, sigma) for mu, sigma in init_params], initial_weights=[0.3, 0.7])

    print(weight, distributions, ll)


if __name__ == '__main__':
    main()
