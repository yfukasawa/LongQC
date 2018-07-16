import abc


class Distribution(object):
    """
    Base class for a mixEM probability distribution.

    To define your own new distribution, all methods of this class will have to be implemented.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def log_density(self, data):
        """Compute the log-probability density :math:`\log P(x|\phi)`

        :param data: The data :math:`x` to compute a probability density for. A :math:`(N \\times D)` :class:`numpy.ndarray` where N is the number of examples and D is the dimensionality of the data
        :type data: numpy.ndarray

        :returns: The log-probability for observing the data, given the probability distribution's parameters
        :rtype: float
        """
        raise NotImplementedError("Need to implement density calculation!")

    @abc.abstractmethod
    def estimate_parameters(self, data, weights):
        """Estimate the probabilities' parameters using weighted maximum-likelihood estimation and update parameters in-place.

        :param data: The data :math:`x` to estimate parameters for. A :math:`(N \\times D)` :class:`numpy.ndarray` where N is the number of examples and D is the dimensionality of the data
        :type data: numpy.ndarray

        :param weights: The weights :math:`\gamma` for individual data points. A N-element :class:`numpy.ndarray` where N is the number of examples.

        Choose those parameters :math:`\phi` that maximize the weighted log-likelihood function:

        .. math::
            ll_\gamma(x|\phi) = \sum_{n=1}^N \gamma_{n} \log [P(x|\phi)]

        Generally, this will involve differentiating the log-likelihood function for all parameters.
        You can set the derivative of the gradient to 0 and try to solve for the parameter to find
        a closed-form solution for the maximum-likelihood estimate or use a numerical optimizer to
        find the maximum-likelihood parameter.

        Once parameter estimates are found, update the attributes *in place*.
        """
        raise NotImplementedError("Need to implement parameter estimation!")

    @abc.abstractmethod
    def __repr__(self):
        """Create a string representation of the probability distribution"""
        raise NotImplementedError("Need to implement string representation!")
