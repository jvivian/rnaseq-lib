from scipy import stats

DISTRIBUTIONS = [
    stats.alpha,
    stats.arcsine,
    stats.betaprime,
    stats.burr,
    stats.chi,
    stats.cosine,
    stats.dweibull,
    stats.expon,
    stats.exponweib,
    stats.f,
    stats.fisk,
    stats.foldnorm,
    stats.genpareto,
    stats.genexpon,
    stats.gausshyper,
    stats.gengamma,
    stats.gilbrat,
    stats.gumbel_r,
    stats.halfcauchy,
    stats.halfnorm,
    stats.hypsecant,
    stats.invgauss,
    stats.johnsonsb,
    stats.ksone,
    stats.laplace,
    stats.levy_l,
    stats.logistic,
    stats.loglaplace,
    stats.lomax,
    stats.mielke,
    stats.ncx2,
    stats.nct,
    stats.pareto,
    stats.powerlaw,
    stats.powernorm,
    stats.reciprocal,
    stats.rice,
    stats.semicircular,
    stats.triang,
    stats.truncnorm,
    stats.uniform,
    stats.vonmises_line,
    stats.weibull_min]


def name_from_dist(dist_func):
    """
    Derives string name from scipy dist function

    :param func dist_func: Scipy distribution function
    :return: Name of function
    :rtype: str
    """
    return str(dist_func).split()[0].split('.')[-1][:-4]
