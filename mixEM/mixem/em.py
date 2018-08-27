import numpy as np

from .progress import logged_simple_progress


def em(data, distributions, initial_weights=None, max_iterations=100, tol=1e-15, tol_iters=10, progress_callback=logged_simple_progress):
    """Fit a mixture of probability distributions using the Expectation-Maximization (EM) algorithm.

    :param data: The data to fit the distributions for. Can be an array-like or a :class:`numpy.ndarray`
    :type data: numpy.ndarray

    :param distributions: The list of distributions to fit to the data.
    :type distributions: list of :class:`mixem.distribution.Distribution`

    :param initial_weights:  Inital weights for the distributions. Must be the same size as distributions. If None, will use uniform initial weights for all distributions.
    :type initial_weights: numpy.ndarray

    :param max_iterations:  The maximum number of iterations to compute for.
    :type max_iterations: int

    :param tol: The minimum relative increase in log-likelihood after tol_iters iterations
    :type tol: float

    :param tol_iters: The number of iterations to go back in comparing log-likelihood change
    :type tol_iters: int

    :param progress_callback: A function to call to report progress after every iteration.
    :type progress_callback: function or None

    :rtype: tuple (weights, distributitons, log_likelihood)
    """

    n_distr = len(distributions)
    n_data = data.shape[0]

    if initial_weights is not None:
        weight = np.array(initial_weights)
    else:
        weight = np.ones((n_distr,))

    last_ll = np.zeros((tol_iters, ))
    resp = np.empty((n_data, n_distr))
    log_density = np.empty((n_data, n_distr))

    iteration = 0
    while True:
        # E-step #######

        # compute responsibilities
        for d in range(n_distr):
            log_density[:, d] = distributions[d].log_density(data)

        # normalize responsibilities of distributions so they sum up to one for example
        resp = weight[np.newaxis, :] * np.exp(log_density)
        resp /= np.sum(resp, axis=1)[:, np.newaxis]

        log_likelihood = np.sum(resp * log_density)

        # M-step #######
        for d in range(n_distr):
            distributions[d].estimate_parameters(data, resp[:, d])

        weight = np.mean(resp, axis=0)

        if progress_callback:
            progress_callback(iteration, weight, distributions, log_likelihood)

        # Convergence check #######
        if np.isnan(log_likelihood):
            last_ll[0] = log_likelihood
            break

        if iteration >= tol_iters and (last_ll[-1] - log_likelihood) / last_ll[-1] <= tol:
            last_ll[0] = log_likelihood
            break

        if iteration >= max_iterations:
            last_ll[0] = log_likelihood
            break

        # store value of current iteration in last_ll[0]
        # and shift older values to the right
        last_ll[1:] = last_ll[:-1]
        last_ll[0] = log_likelihood

        iteration += 1

    return weight, distributions, last_ll[0]
