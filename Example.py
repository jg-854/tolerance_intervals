"""
@author: jeevi
description: an example showing the use of the tolerances.py module for both one sided and double sided data.
"""

from numpy import genfromtxt
import tolerances

conf = 0.8  # this is the confidence that beta*100 % of a sample will be within a certain interval
coverage = 0.95  # this is the coverage, this value will change

# example showing two sided tolerance interval

data = genfromtxt("two_sided_data.csv", delimiter=',')
updated_beta, upper, lower, _, _, _ = tolerances.double_sided_tolerance(conf, coverage, data)

print("Over {beta_val} % of the population covered, with {confidence_val} % confidence that values lie between {"
      "lower_val} and {upper_val}".format(beta_val=updated_beta.round(2) * 100, confidence_val=conf * 100,
                                          lower_val=lower.round(4), upper_val=upper.round(4)))

# example showing one sided tolerance interval

data = genfromtxt("one_sided_data.csv", delimiter=',')
updated_beta, upper, _, _, _ = tolerances.one_sided_tolerance(conf, coverage, data)

print("Over {beta}% of the population covered, with {confidence}% confidence that values are under {upper_val}".format(
    beta=updated_beta.round(2) * 100, confidence=conf * 100, upper_val=upper.round(4)))


