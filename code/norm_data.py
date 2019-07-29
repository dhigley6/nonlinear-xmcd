"""Calculating normalized/weighted data

Daniel Higley
"""

import numpy as np

def norm_disc_scan(bin_data, sigs, norms):
    """Return normalized binned data for a scan with discrete steps
    """
    points = np.unique(bin_data)
    space_between_points_min = np.amin(np.diff(points))
    bin_edges_except_last = points-space_between_points_min/2
    bin_edge_last = points[-1]+space_between_points_min/2
    bin_edges = np.append(bin_edges_except_last, bin_edge_last)
    norm_data = norm_trace(bin_data, sigs, norms, bin_edges)
    norm_data['points'] = points
    return norm_data

def norm_trace(bin_data, sigs, norms, bins=100):
    """Return normalized binned data
    """
    # We calculate the normalized data of each bin as sum(signal)/sum(norm).
    # This is equivalent to calculating the weighted mean of signal, with
    # weights of norm, as done below:
    data = weight_trace(bin_data, sigs, norms, bins)
    return data

def weight_trace(bin_data, weighted_sigs, weights, bins=100):
    """Return weighted mean and standard deviation of binned input data
    """
    sigs = weighted_sigs/weights
    if not hasattr(bins, '__getitem__'):
        bin_start = np.amin(bin_data)
        bin_end = np.amax(bin_data)
        bins = np.linspace(bin_start, bin_end, bins)
    bin_centers = (bins[1:]+bins[:-1])/2
    bin_counts, unused = np.histogram(bin_data, bins=bins)
    weight_sig_sums, unused = np.histogram(bin_data, weights=weighted_sigs,
                                           bins=bins)
    weight_sums, unused = np.histogram(bin_data, weights=weights, 
                                       bins=bins)
    weight_sig_mean = weight_sig_sums/weight_sums
    weight_sig_squared_sums, unused = np.histogram(bin_data,
                                                   weights=sigs**2*weights,
                                                   bins=bins)
    weight_sig_squared_mean = weight_sig_squared_sums/weight_sums
    weight_sig_var = (weight_sig_squared_mean-weight_sig_mean**2)
    # weight_sig_std is the weighted standard deviation and is
    # computed following the formula on the Wikipedia entry for
    # mean square weighted deviation
    weight_sig_std = np.sqrt(weight_sig_var)
    weight_sig_std = weight_sig_std/np.sqrt(bin_counts)
    data = {'bin_edges': bins,
            'bin_centers': bin_centers,
            'norm_sig': weight_sig_mean,
            'err': weight_sig_std,
            'norm_sums': weight_sums,
            'bin_counts': bin_counts}
    return data
