import numpy as np


def probability(data, weights, distributions):
    """Compute the probability for data of the mixture density model given by weights and a list of distributions"""

    assert len(weights) == len(distributions), "Need matching number of weights and distributions!"

    if not hasattr(data, '__len__'):
        data = [data]

    data = np.array(data)

    n_dist = len(distributions)

    densities = np.empty((data.shape[0], n_dist))
    for d in range(n_dist):
        densities[:, d] = np.exp(distributions[d].log_density(data))

    return np.sum(weights[np.newaxis, :] * densities, axis=1)
