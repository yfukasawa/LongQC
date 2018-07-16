Getting started
===============

Design
------

From :math:`N` points of data :math:`X=(\vec x_1, \ldots, \vec x_N)`, mix'EM will let you fit a probability distribution :math:`P(x|\phi)` from a mixture of :math:`D` probability distributions with parameters :math:`\phi_d` and mixing weights :math:`w_d`:

.. math::
    P(x|\phi) = \sum_{d=1}^D w_d P(x|\phi_d), \; \sum_{d=1}^D w_d = 1

To fit the data to your user-defined model, you pass both the data :math:`X` and a :py:class:`list` of :py:class:`mixem.distribution.Distribution` objects with initial parameters to the :py:func:`mixem.em` function. Mix'EM will optimize the model parameters to best fit the provided data and return the mixing weights :math:`w_d` and distributions with parameters :math:`\phi_d` plus the final log-likelihood:

.. code-block:: python

    weights, distributions, log_likelihood = mixem.em(my_data, [
        mixem.distribution.NormalDistribution(mu=1, sigma=1),
        mixem.distribution.ExponentialDistribution(lmbda=1),
    ], initial_weights=[0.5, 0.5])

    print(("Final model: {w0} * Norm[mu={mu}, sigma={sigma}] +" 
           " {w1} * Exp[lambda={lmbda}]").format(
        w0=weights[0], w1=weights[1],
        mu=distributions[0].mu, sigma=distributions[0].sigma,
        lmbda=distributions[1].lmbda
    ))


Examples
--------

A good walkthrough of all major functions of mix'EM is given in the :doc:`Old Faithful example<examples/old_faithful>`.

If you have a particular set of distributions that you'd like to fit, you can check out the examples in the `examples <https://github.com/sseemayer/mixem/tree/master/examples>`_ subdirectory of the mixem project.

