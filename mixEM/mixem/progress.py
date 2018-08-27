# coding=utf-8
from logging import getLogger
logger = getLogger(__name__)

def simple_progress(iteration, weights, distributions, log_likelihood):
    """A simple default progress callback to use with mixem.em"""

    print("iteration {iteration:4d} (log-likelihood={log_likelihood:.5e}): p(x|Φ) = {formatted_distributions}".format(
        iteration=iteration,
        log_likelihood=log_likelihood,
        formatted_distributions=" + ".join("{w:.3g}*{d}".format(w=w, d=d) for w, d in zip(weights, distributions))
    ))

def logged_simple_progress(iteration, weights, distributions, log_likelihood):

    logger.info("iteration {iteration:4d} (log-likelihood={log_likelihood:.5e}): p(x|Φ) = {formatted_distributions}".format(
        iteration=iteration,
        log_likelihood=log_likelihood,
        formatted_distributions=" + ".join("{w:.3g}*{d}".format(w=w, d=d) for w, d in zip(weights, distributions))
    ))
