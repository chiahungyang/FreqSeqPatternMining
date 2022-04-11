"""To compute support thresholds of significance of DNA sequential patterns."""


import numpy as np
from scipy.special import gammaincc
from scipy.optimize import root


class RootFailure(Exception):
    """Exception that an attempt to find numerical root fails."""
    pass


def significant_support(size, _len, conf):
    """
    Return the support threshold for a pattern of given length such that the
    occurence indicates a significantly rare event in a null model where
    patterns are drawn uniformly from all possibilities conditioned on length.

    Parameters
    ----------
    size: int
        Size of hypothetical samples from the null model.

    _len: int
        Given length of sequential patterns.

    conf: float
        Confidence level of significance.

    Return
    ------
    thrd: int
        Support threshold of significance, that is the minimum support at which
        the cumulative distribution function is greater than the confidence
        level.

    Note
    ----
    The probability distribution of the null model is exact binomial. For
    simplicity computing the cumulative distribution, it is assumed that the
    sample size is large and the distribution is approximately Poisson, whose
    cumulative form is a regularized Gamma function.

    """
    prob = 0.25 ** _len  # Probability that any specific pattern is realized
    mean = size * prob

    def cdf(x):
        return gammaincc(x, mean)

    # Solve for the threshold such that the cumulative function equates to the
    # confidence level
    def func(x):
        return cdf(x) - conf
    guess = mean + np.sqrt(mean)
    sol = root(func, guess)

    # Return the threshold if root finding was successful
    # Otherwise raise a RootFailure exception
    if sol.success:
        thrd = np.ceil(sol.x[0])
        return thrd
    else:
        raise RootFailure


def thresholds(num, min_len, lwr_bd, conf):
    """
    Return support thresholds of significance for patterns from a given minimum
    length up to the one such that the support threshold is not less than a
    pre-specified lower bound.

    Parameters
    ----------
    num: int
        Number of nucleotides (characters) in the genomic data.

    min_len: int
        Minimum pattern length of interests.

    lwr_bd: int
        Pre-specified lower bound of support thresholds.

    conf: float
        Confidence level of significance.

    Return
    ------
    len_thrd: dict
        Dictinoary mapping from pattern length to the corresponding support
        thresholds.

    """
    len_thrd = dict()
    _len = min_len
    reach_bd = False  # Indicator of whether the lower bound is reached
    while not reach_bd:
        thrd = significant_support(num - _len, _len, conf)
        if thrd >= lwr_bd:
            len_thrd[_len] = thrd
        else:
            reach_bd = True
        _len += 1

    return len_thrd


def main():
    """Empty main function."""
    return


if __name__ == '__main__':
    main()
