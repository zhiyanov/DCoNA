"""Statistics related to pearson/spearman correlations

The module allows to compute pearson/spearman correlations,
their quantiles, mean values and standard deviations.

Also it can be used to compute the probability that
one correltion differs from another one no more then
a given threshold
"""


import numpy as np
import scipy.stats
from ..utils import \
    bound, bootstrap_sample


PARTITION_SIZE = 10**3

EPS1 = 1e-5
EPS2 = 1e-2
EPS3 = 1e-5


# Spearmanr block
def spearmanr_std(rs, size):
    """Computes standard deviation of spearman correlation

    Parameters
    ----------
    rs: numerical value or numpy.array
        Sperman correlation(s).
    size: numerical value or numpy.array
        Size(s) of sample(s) that was(were) used to compute
        rs correlation(s).
    Returns
    -------
    numerical value or numpy.array respectively to the input
        Assuming that "normalized" spearman correlation has
        a distribution t(df=size - 2),
        it returns standard deviation (std) of arctanh(rs).

        Normlizaion: "(arctanh(x) - arctanh(rs)) / std" has
        "t(df=size - 2)" distribution, where x is a random variable
        distributed as the correlation with an actual value rs.
    """
    return np.sqrt((1 + rs**2 / 2) / (size - 3))

def spearmanr_mean(source, target):
    """Computes actual value of spearman correlation

    Parameters
    ----------
    source: list
        First sample.
    target: list
        Second sample.
    Returns
    -------
    numerical value
        Spearman correlation between source and target.
    """
    return scipy.stats.spearmanr(source, target)[0]

def spearmanr_cdf(quantiles, rs, size, ss=None):
    """Computes cdf of spearman correlation distribution

    Parameters
    ----------
    quantiles: numerical value or numpy.array
    rs: numerical value or numpy.array
        Sperman correlation(s).
    size: numerical value or numpy.array
        Size(s) of sample(s) that were used to compute
        "rs".
    ss: numerical value or numpy.array
        Std of "arctanh(rs)" computed by "spearman_std".
    Returns
    -------
    numerical value or numpy.array respectively to "quantiles"
        Assuming that "normalized" spearman correlation has
        a distribution t(df=size - 2),
        it computes quantiles of spearman correlation.

        Normlizaion: "(arctanh(x) - arctanh(rs)) / ss" has
        "t(df=size - 2)" distribution, where x is a random variable
        distributed as the correlation with an actual value rs.
    """
    if not ss:
        ss = spearmanr_std(rs, size)
    return scipy.stats.t.cdf(
        (np.arctanh(quantiles) - np.arctanh(rs)) / ss,
        size - 2
    )

def spearmanr_ppf(proba, rs, size):
    """Computes quantiles of spearman correlation distribution

    Parameters
    ----------
    proba: numerical value or numpy.array
    rs: numerical value or numpy.array
        Sperman correlation(s).
    size: numerical value or numpy.array
        Size(s) of sample(s) that were used to compute
        "rs".
    Returns
    -------
    Quantiles corrsponding to probas
    """
    ss = spearmanr_std(rs, size)  
    z_a = scipy.stats.norm.ppf(proba)

    return np.tanh(np.arctanh(rs) + ss * z_a)

def spearmanr_conf_interval(confidence, rs, size):
    """Computes confidence interval of spearman correlation

    Parameters
    ----------
    confidence: numerical value or numpy.array
    rs: numerical value or numpy.array
        Sperman correlation(s).
    size: numerical value or numpy.array
        Size(s) of sample(s) that were used to compute
        "rs".
    Returns
    -------
    tuple of numerical value or numpy.array respectively to "confidence" and "rs"
    """
    return (spearmanr_ppf((1 - confidence) / 2, rs, size),
            spearmanr_ppf((1 + confidence) / 2, rs, size))
           
# Pearsonr block
def pearsonr_std(rs, size):
    """Computes standard deviation of pearson correlation

    Parameters
    ----------
    rs: numerical value or numpy.array
        Pearson correlation(s).
    size: numerical value or numpy.array
        Size(s) of sample(s) that were used to compute
        "rs" correlations.
    Returns
    -------
    numerical value or numpy.array respectively to the input
        Assuming that a "normalized" pearson correlation has
        a distribution t(df=size - 2),
        it returns standard deviation (std) of "arctanh(rs)".

        Normlizaion: "(arctanh(x) - arctanh(rs)) / std" has
        normal distribution, where "x" is a random variable
        distributed as the correlation with an actual value "rs".
    """
    return 1 / np.sqrt(size - 3)

def pearsonr_mean(source, target):
    """Computes actual value of pearson correlation

    Parameters
    ----------
    source: numpy.array
        First sample.
    target: numpy.array
        Second sample.
    Returns
    -------
    numerical value
        Pearson correlation between "source" and "target".
    """
    return scipy.stats.pearsonr(source, target)[0]

def pearsonr_cdf(quantiles, rs, size, ss=None):
    """Computes cdf of pearson correlation distribution

    Parameters
    ----------
    quantiles: numerical value or numpy.array
    rs: numerical value or numpy.array
        Pearson correlation(s).
    size: numerical value or numpy.array
        Size(s) of sample(s) that were used to compute
        "rs" correlation(s).
    ss: numerical value or numpy.array
        Std(s) computed by "pearsonr_std".
    Returns
    -------
    numerical value or numpy.array respectively to "quantiles"
        Assuming that "normalized" pearson correlation has
        normal distribution, it computes quantiles of
        pearson correlation.

        Normlizaion: "(arctanh(x) - arctanh(rs)) / ss" has
        normal distribution, where x is a random variable
        distributed as the correlation with an actual value rs.
    """
    if not ss:
        ss = pearsonr_std(rs, size)
    return scipy.stats.norm.cdf(
        (np.arctanh(quantiles) - np.arctanh(rs)) / ss
    )

def pearsonr_ppf(proba, rs, size):
    """Computes quantiles of pearson correlation distribution

    Parameters
    ----------
    proba: numerical value or numpy.array
    rs: numerical value or numpy.array
        Sperman correlation(s).
    size: numerical value or numpy.array
        Size(s) of sample(s) that were used to compute
        "rs".
    Returns
    -------
    Quantiles corrsponding to probas
    """
    ss = pearsonr_std(rs, size)  
    z_a = scipy.stats.norm.ppf(proba)

    return np.tanh(np.arctanh(rs) + ss * z_a)

def pearsonr_conf_interval(confidence, rs, size):
    """Computes confidence interval of pearson correlation

    Parameters
    ----------
    confidence: numerical value or numpy.array
    rs: numerical value or numpy.array
        Sperman correlation(s).
    size: numerical value or numpy.array
        Size(s) of sample(s) that were used to compute
        "rs".
    Returns
    -------
    tuple of numerical value or numpy.array respectively to "confidence" and "rs"
    """
    return (pearsonr_ppf((1 - confidence) / 2, rs, size),
            pearsonr_ppf((1 + confidence) / 2, rs, size))


# Correlation difference block
def correlation_diff_cdf(
    first_source, first_target,
    second_source, second_target,
    delta,
    correlation="spearman",
    alternative="two-sided",
    method="analytic"
):
    """Computes the probability that differnce of
    correlations is lower than
    the given threshold "delta"

    Parameters
    ----------
    first_source: numpy.array
    first_target: numpy.array
        These arguments are used to compute
        the first correlation and have equal size(shape).
    second_source: numpy.array
    second_target: numpy.array
        These arguments are used to compute
        the second correlation and have equal size(shape).
    delta: numerical value or numpy.array
    correlation: "spearman" (default), "pearson"
        Type of the correlation which is used
        in the probability computation.
    alternative: "two-sided" (default), "less", "greater"
        Computes the probability of the following events:
        "two-sided" "|rs1 - rs2| <= delta",
        "less" "rs1 <= rs2 + delta",
        "greater" "rs1 > rs2 + delta".
    method: "analytic" (default), "bootstrap"
        Defines the way of the probability computation.
    Returns
    -------
    numerical value or numpy.array of "delta" shape
        Computes the probability of "alternative"
    """
    if (correlation == "spearman"):
        correlation_std = spearmanr_std
        correlation_mean = spearmanr_mean
        correlation_cdf = spearmanr_cdf
    elif (correlation == "pearson"):
        correlation_std = pearsonr_std
        correlation_mean = pearsonr_mean
        correlation_cdf = pearsonr_cdf

    if (isinstance(delta, float)):
        delta = np.array([delta])
    else:
        delta = np.array(delta)

    first_target = np.array(first_target)
    second_target = np.array(second_target)

    proba = np.zeros(len(first_target))

    if (method == "bootstrap"):
        first_rs_targets = []
        second_rs_targets = []
        for i in range(len(first_target)):
            first_rs_targets.append([
                r for r in bootstrap_sample(
                    first_source,
                    first_target[i],
                    statistic=correlation_mean
                )
            ])

        for i in range(len(second_target)):
            second_rs_targets.append([
                r for r in bootstrap_sample(
                    second_source,
                    second_target[i],
                    statistic=correlation_mean
                )
            ])

        first_rs_targets = np.array(first_rs_targets)
        second_rs_targets = np.array(second_rs_targets)

        if (alternative == "greater"):
            first_rs_targets, second_rs_targets =\
                second_rs_targets, first_rs_targets
            delta = -delta
            alternative = "less"

        return diff_bootstrap_cdf(
            first_rs_targets,
            second_rs_targets,
            delta,
            alternative=alternative
        )

    elif (method == "analytic"):
        first_rs = []
        second_rs = []

        for i in range(len(first_target)):
            rs = correlation_mean(
                first_source, first_target[i]
            )
            first_rs.append(rs)

        for i in range(len(second_target)):
            rs = correlation_mean(
                second_source, second_target[i]
            )
            second_rs.append(rs)

        first_rs = np.array(first_rs)
        second_rs = np.array(second_rs)

        if (alternative == "greater"):
            first_rs, second_rs =\
                second_rs, first_rs
            delta = -delta
            alternative = "less"

        return correlation_diff_analytic_cdf(
            first_rs, len(first_source),
            second_rs, len(second_source),
            delta,
            correlation=correlation,
            alternative=alternative
        )

    return None

def correlation_diff_analytic_cdf(
    first_rs, first_size,
    second_rs, second_size,
    delta,
    correlation="spearman",
    alternative="two-sided"
):
    """Computes the probability that the differnce
    of correlations is lower than
    the given threshold "delta"

    Parameters
    ----------
    first_rs: numerical value or numpy.array
        Correlation(s).
    first_size: numerical value or numpy.array
        Size(s) of sample(s) that were used
        to compute correlation(s) "first_rs".
    second_rs: numerical value or numpy.array
        Correlations.
    second_size: numerical value or numpy.array
        Size(s) of sample(s) that were used
        to compute correlation(s) "second_rs".
    delta: numerical value or numpy.array
    correlation: "spearman" (default), "pearson"
        Type of the correlation which is used
        in the probability computation.
    alternative: "two-sided" (default), "less"
        Computes the probability of the following events:
        "two-sided" "|rs1 - rs2| <= delta",
        "less" "rs1 <= rs2 + delta".
    Returns
    -------
    numerical value or numpy.array of "delta" shape
        Computes the probability of "alternative".
    """
    proba = np.zeros(len(delta))

    indexes_1 = np.logical_and(np.abs(first_rs) < 1 - EPS2,
        np.abs(second_rs) < 1 - EPS2)
    indexes_2 = np.logical_xor(np.abs(first_rs) >= 1 - EPS2,
        np.abs(second_rs) >= 1 - EPS2)
    indexes_3 = np.logical_and(np.abs(first_rs) >= 1 - EPS2,
        np.abs(second_rs) >= 1 - EPS2)

    if (np.sum(indexes_1) > 0):
        proba[indexes_1] = _correlation_diff_analytic_cdf_1(
            first_rs[indexes_1], first_size,
            second_rs[indexes_1], second_size,
            delta[indexes_1],
            correlation=correlation,
            alternative=alternative
        )

    if (np.sum(indexes_2) > 0):
        proba[indexes_2] = _correlation_diff_analytic_cdf_2(
            first_rs[indexes_2], first_size,
            second_rs[indexes_2], second_size,
            delta[indexes_2],
            correlation=correlation,
            alternative=alternative
        )

    if (np.sum(indexes_3) > 0):
        proba[indexes_3] = _correlation_diff_analytic_cdf_3(
            first_rs[indexes_3], first_size,
            second_rs[indexes_3], second_size,
            delta[indexes_3],
            correlation=correlation,
            alternative=alternative
        )

    return proba

def _correlation_diff_analytic_cdf_1(
    first_rs, first_size,
    second_rs, second_size,
    delta,
    correlation="spearman",
    alternative="two-sided"
):
    """
    Similar to "correlation_diff_analytic_cdf". It
    is used in the case of both "first_rs" and "second_rs"
    are separated from 1 and -1.
    """
    # check that all correlations are not equal to 1
    assert(np.sum(np.logical_or(np.abs(first_rs) >= 1 - EPS2,
        np.abs(second_rs) >= 1 - EPS2)) == 0)

    if (correlation == "spearman"):
        correlation_std = spearmanr_std
        correlation_mean = spearmanr_mean
        correlation_cdf = spearmanr_cdf
    elif (correlation == "pearson"):
        correlation_std = pearsonr_std
        correlation_mean = pearsonr_mean
        correlation_cdf = pearsonr_cdf

    delta = delta.reshape((-1, 1))

    quantiles = np.linspace(-1, 1, PARTITION_SIZE).reshape((1, -1))
    quantiles_pd = quantiles + delta
    quantiles_md = quantiles - delta

    quantiles = bound(quantiles, -1 + EPS1, 1 - EPS1)
    quantiles_pd = bound(quantiles_pd, -1 + EPS1, 1 - EPS1)
    quantiles_md = bound(quantiles_md, -1 + EPS1, 1 - EPS1)

    first_ss = correlation_std(first_rs, first_size)
    second_ss = correlation_std(second_rs, second_size)

    first_rs = first_rs.reshape((-1, 1))
    first_ss = first_ss.reshape((-1, 1))
    second_rs = second_rs.reshape((-1, 1))
    second_ss = second_ss.reshape((-1, 1))

    fd_pd = correlation_cdf(quantiles_pd, first_rs, first_size, ss=first_ss)
    fd_md = correlation_cdf(quantiles_md, first_rs, first_size, ss=first_ss)
    sd = correlation_cdf(quantiles, second_rs, second_size, ss=second_ss)

    if (alternative == "less"):
        return np.sum(
            fd_pd[:, :-1] * (sd[:, 1:] - sd[:, :-1]),
            axis=1
        )

    return np.sum(
        (fd_pd[:, :-1] - fd_md[:, :-1]) * (sd[:, 1:] - sd[:, :-1]),
        axis=1
    )

def _correlation_diff_analytic_cdf_2(
    first_rs, first_size,
    second_rs, second_size,
    delta,
    correlation="spearman",
    alternative="two-sided"
):
    """
    Similar to "correlation_diff_analytic_cdf". It
    is used in the case: only one of "first_rs" and "second_rs"
    is near to -1 or 1.
    """
    # check that all correlations are equal to 1
    assert(np.sum(np.logical_xor(np.abs(first_rs) >= 1 - EPS2,
        np.abs(second_rs) >= 1 - EPS2)) == len(delta))

    if (correlation == "spearman"):
        correlation_std = spearmanr_std
        correlation_mean = spearmanr_mean
        correlation_cdf = spearmanr_cdf
    elif (correlation == "pearson"):
        correlation_std = pearsonr_std
        correlation_mean = pearsonr_mean
        correlation_cdf = pearsonr_cdf

    if (alternative == "less"):
        proba = np.zeros(len(delta))

        first_ss = correlation_std(first_rs, first_size)
        second_ss = correlation_std(second_rs, second_size)

        f_indexes = (np.abs(first_rs) >= 1 - EPS2)
        s_indexes = (np.abs(second_rs) >= 1 - EPS2)

        proba[f_indexes] = 1 - correlation_cdf(
            bound(first_rs[f_indexes] - delta[f_indexes], -1 + EPS1, 1 - EPS1),
            second_rs[f_indexes], second_size, ss=second_ss[f_indexes]
        )

        proba[s_indexes] = correlation_cdf(
            bound(second_rs[s_indexes] + delta[s_indexes], -1 + EPS1, 1 - EPS1),
            first_rs[s_indexes], first_size, ss=first_ss[s_indexes]
        )

        return proba

    indexes = (np.abs(first_rs) >= 1 - EPS2)
    first_rs[indexes], second_rs[indexes] =\
        second_rs[indexes], first_rs[indexes]

    first_ss = correaltion_std(first_rs, first_size)

    proba = np.zeros(len(delta))

    p_indexes = (second_rs >= 1 - EPS2)
    m_indexes = (second_rs <= -1 + EPS2)

    proba[p_indexes] = 1 - correlation_cdf(
        bound(1 - delta[p_indexes], -1 + EPS1, 1 - EPS1),
        first_rs[p_indexes], first_size, ss=first_ss[p_indexes]
    )

    proba[m_indexes] = correlation_cdf(
        bound(-1 + delta[m_indexes], -1 + EPS1, 1 - EPS1),
        first_rs[m_indexes], first_size, ss=first_ss[m_indexes]
    )

    return proba

def _correlation_diff_analytic_cdf_3(
    first_rs, first_size,
    second_rs, second_size,
    delta,
    correlation="spearman",
    alternative="two-sided"
):
    """
    Similar to "correlation_diff_analytic_cdf". It
    is used in the case of both "first_rs"" and "second_rs"
    are close to 1 or -1.
    """
    # check that all correlations are equal to 1
    assert(np.sum(np.logical_and(np.abs(first_rs) < 1 - EPS2,
        np.abs(second_rs) < 1 - EPS2)) == 0)

    proba = np.zeros(len(delta))
    if (alternative == "less"):
        proba[first_rs <= second_rs + delta] = 1
        return proba

    proba[np.abs(first_rs - second_rs) <= delta] = 1
    return proba

def diff_bootstrap_cdf(
    first_sample, second_sample,
    delta,
    alternative="two-sided"
):
    """Computes the probability that the differnce
    "first_sample - second_sample" is lower than
    treshold "delta"

    Parameters
    ----------
        first_sample: numpy.ndarray
        second_sample: numpy.ndarray
            The previous argements have equal shape.
            Their rows are samples (for example,
            received by bootstrap method) of equal size.
        alternative: "two-sided" (default), "less"
            Computes the probability of the following events:
            "two-sided" "|first_sample - second_sample| <= delta",
            "less" "first_sample <= second_sample + delta".
    Returns
    -------
    numerical value or numpy.array of delta.shape
        Computes the probability of "alternative"
    """
    delta = delta.reshape((-1, 1))
    if (alternative == "less"):
        return np.sum(
            first_sample <= second_sample + delta,
            axis=1
        ) / len(first_sample[0])

    return np.sum(
        np.abs(first_sample - second_sample) <= delta,
        axis=1
    ) / len(first_sample[0])
