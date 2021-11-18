"""
This package includes tools for statistical tests.
"""

from statsmodels.stats import proportion


def one_proportions_z_test(prop1, n1, expected):

    stat, pval = proportion.proportions_ztest(
            int(prop1),
            int(n1),
            value=expected,
            alternative='two-sided',
            prop_var=False)

    print('{0:0.1000f}'.format(pval))
    return stat, pval
