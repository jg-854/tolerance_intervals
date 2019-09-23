"""
@author: jeevi
description: A module that will calculate tolerance intervals based on coverage and confidence for a given dataset.
"""
from scipy.stats import chi2
from scipy.stats import norm
import math
import numpy as np
from numpy import percentile


def find_nearest(array, value):
    # nearest neighbour interpolation to find the input value which maps to the CDF's output value
    idx = (np.abs(array - value)).argmin()
    return idx


def double_sided_tolerance(alpha, beta, dataset, bootstrap_iterations=100):
    norm_val = (beta + 1) / 2
    num_samples = len(dataset)
    dof = num_samples - 1  # degrees of freedom, or number of residuals

    # these are the empirical mean and std from the given sample, as opposed to the bootstrapped values seen later
    sample_mu = dataset.mean()
    sample_std = dataset.std()

    z = norm.ppf(norm_val)
    chi_sq = chi2.ppf(1 - alpha, dof)

    k = z * math.sqrt((dof * (1 + 1 / num_samples)) / chi_sq)
    w = math.sqrt(1 + (num_samples - 3 - chi2.ppf(1 - alpha, dof)) /
                  (2 * (num_samples + 1) * (num_samples + 1)))

    # w is very close to 1 but we will still use the correction just because it doesn't add much computational time
    k_corrected = w * k

    # these are the tolerance intervals ASSUMING our data is normally distributed
    lower_bound = sample_mu - (k_corrected * sample_std)
    upper_bound = sample_mu + (k_corrected * sample_std)

    dataset.sort()  # this forms an empirical cdf of the original dataset
    p = np.linspace(0, 0.995, num_samples)  # this generates the probability (y axis) for the cdf
    d = []
    for i in range(bootstrap_iterations):  # value is arbitrary: the higher the better the accuracy
        # the mean and std are calculated for the bootstrapped sample

        bsample = []
        for i in range(num_samples):
            bsample.append(np.random.choice(dataset, replace=True))
        bsample = np.asarray(bsample)
        bmu = bsample.mean()
        bstd = bsample.std()
        bsample.sort()

        f_sam_upper = p[find_nearest(bsample, bmu + (k_corrected * bstd))]
        f_sam_lower = p[find_nearest(bsample, bmu - (k_corrected * bstd))]
        f_emp_upper = p[find_nearest(dataset, bmu + (k_corrected * bstd))]
        f_emp_lower = p[find_nearest(dataset, bmu - (k_corrected * bstd))]

        db = math.pow(num_samples, 0.5) * (f_sam_upper - f_sam_lower - (f_emp_upper - f_emp_lower))
        d.append(db)

    pcnt = percentile(d, alpha * 100)
    updated_beta = p[find_nearest(dataset, upper_bound)] - p[find_nearest(dataset, lower_bound)] - pcnt / math.sqrt(
        num_samples)

    return updated_beta, lower_bound, upper_bound, sample_mu, sample_std, num_samples


def one_sided_tolerance(alpha, beta, dataset, bootstrap_iterations=100):
    num_samples = len(dataset)

    # these are the empirical mean and std from the given sample, as opposed to the bootstrapped values seen later
    sample_mu = dataset.mean()
    sample_std = dataset.std()

    zp = norm.ppf(beta)
    za = norm.ppf(alpha)
    a = 1 - 0.5 * za * za / (num_samples - 1)
    b = zp * zp - za * za / num_samples

    k = (zp + math.pow(zp * zp - a * b, 0.5)) / a

    # these are the tolerance intervals ASSUMING our data is normally distributed
    upper_bound = sample_mu + (k * sample_std)

    dataset.sort()  # this forms an empirical cdf of the original dataset

    p = np.linspace(0, 0.995, num_samples)  # this generates the probability (y axis) for the cdf

    d = []
    for i in range(bootstrap_iterations):  # value is arbitrary: the higher the better the accuracy
        # the mean and std are calculated for the bootstrapped sample
        bsample = []
        for i in range(num_samples):
            bsample.append(np.random.choice(dataset, replace=True))
        bsample = np.asarray(bsample)
        bmu = bsample.mean()
        bstd = bsample.std()
        bsample.sort()

        f_sam_upper = p[find_nearest(bsample, bmu + (k * bstd))]
        f_emp_upper = p[find_nearest(dataset, bmu + (k * bstd))]

        db = math.pow(num_samples, 0.5) * (f_sam_upper - f_emp_upper)
        d.append(db)

    pcnt = percentile(d, alpha * 100)
    updated_beta = p[find_nearest(dataset, upper_bound)] - pcnt / math.sqrt(num_samples)

    return updated_beta, upper_bound, sample_mu, sample_std, num_samples
