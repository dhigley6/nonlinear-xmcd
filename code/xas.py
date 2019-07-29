"""Module for normalizing/calibrating XAS spectra
"""

import numpy as np

import matplotlib.pyplot as plt
from scipy.optimize import leastsq

def get_bg_from_ref(phot_sam, spec_sam, phot_ref, spec_ref, valid_phots):
    # interpolate reference spectrum at photon energy locations of sample
    # spectrum:
    ref_interp = np.interp(phot_sam, phot_ref, spec_ref)
    sam_minus_ref = spec_sam-ref_interp
    phot_in_range = []
    for phot in phot_sam:
        curr_phot_in_range = []
        for valid_phot in valid_phots:
            curr_phot_in_range.append((phot >= valid_phot[0]) & (phot <= valid_phot[1]))
        curr_phot_in_range = np.any(curr_phot_in_range)
        phot_in_range.append(curr_phot_in_range)
    phot_reduced = phot_sam[np.array(phot_in_range)]
    sam_minus_ref_reduced = sam_minus_ref[np.array(phot_in_range)]
    lin_fit = np.polyfit(phot_reduced, sam_minus_ref_reduced, deg=1)
    return lin_fit

def norm_to_no_sample(phot_sam, spec_sam,
                      phot_no_sam, spec_no_sam, fit_deg=1):
    """Normalize the sample spectrum to a reference no sample spectrum
    """
    # Fit the reference no sample spectrum with NaNs removed
    not_nan_idx = ~np.isnan(spec_no_sam)
    phot_no_sam_no_nans = phot_no_sam[not_nan_idx]
    spec_no_sam_no_nans = spec_no_sam[not_nan_idx]
    no_sam_fit = np.polyfit(phot_no_sam_no_nans, spec_no_sam_no_nans, deg=fit_deg)
    # Divide the sample spectrum by the fit of the reference spectrum to
    # get the normalized spectrum
    plt.figure()
    plt.plot(phot_no_sam_no_nans, spec_no_sam_no_nans)
    spec_no_sam_fitted = np.polyval(no_sam_fit, phot_sam)
    spec_normalized = spec_sam-spec_no_sam_fitted
    return spec_normalized

def sub_preedge(phot, raw_xas, phot_bounds_preedge):
    """Subtract a linear fit to the pre-edge region
    """
    preedge_poly = fit_preedge(phot, raw_xas, phot_bounds_preedge)
    proc_xas = raw_xas-preedge_poly(phot)
    return proc_xas

def fit_preedge(phot, raw_xas, phot_bounds_preedge):
    """Return a linear fit to the pre-edge region
    """
    preedge_start = np.amin(phot_bounds_preedge)
    preedge_end = np.amax(phot_bounds_preedge)
    preedge_region = (phot >= preedge_start) & (phot <= preedge_end)
    preedge_phot = phot[preedge_region]
    preedge_xas = raw_xas[preedge_region]
    preedge_phot_nonans = preedge_phot[~np.isnan(preedge_xas)]
    preedge_xas_nonans = preedge_xas[~np.isnan(preedge_xas)]
    preedge_fit = np.polyfit(preedge_phot_nonans, preedge_xas_nonans, deg=1)
    preedge_poly = np.poly1d(preedge_fit)
    return preedge_poly

def edgejump_norm(phot, raw_xas, phot_bounds_preedge, phot_bounds_postedge,
                  edge_jump=1):
    """Subtract a linear fit to the preedge region and set the edge jump
    """
    preedge_start = np.amin(phot_bounds_preedge)
    preedge_end = np.amax(phot_bounds_preedge)
    preedge_region = (phot >= preedge_start) & (phot <= preedge_end)
    postedge_start = np.amin(phot_bounds_postedge)
    postedge_end = np.amax(phot_bounds_postedge)
    postedge_region = (phot >= postedge_start) & (phot <= postedge_end)
    preedge_phot = phot[preedge_region]
    preedge_xas = raw_xas[preedge_region]
    preedge_fit = np.polyfit(preedge_phot, preedge_xas, deg=1)
    preedge_poly = np.poly1d(preedge_fit)
    proc_xas = raw_xas-preedge_poly(phot)
    postedge_avg = np.mean(proc_xas[postedge_region])
    #proc_xas = proc_xas*edge_jump/postedge_avg
    return proc_xas

def fit_abs(phot, spec, cal_phot, cal_spec,
            phot_scale_guess=1, phot_shift_guess=0,
            spec_scale_guess=1, spec_shift_guess=0,
            spec_slope_guess=0):
    """Fit one spectrum to another (e.g., for calibration)
    If one sets a guess parameter to None, then it will not be used
    to fit and will simply be set to its default value.
    """
    if phot_scale_guess is None:
        phot_scale_guess = 1
        do_phot_scale = False
    else:
        do_phot_scale = True
    if phot_shift_guess is None:
        phot_shift_guess = 0
        do_phot_shift = False
    else:
        do_phot_shift = True
    if spec_scale_guess is None:
        spec_scale_guess = 1
        do_spec_scale = False
    else:
        do_spec_scale = True
    if spec_shift_guess is None:
        spec_shift_guess = 0
        do_spec_shift = False
    else:
        do_spec_shift = True
    if spec_slope_guess is None:
        spec_slope_guess = 0
        do_spec_slope = False
    else:
        do_spec_slope = True

    def trans_spec_diffs(fit_coef):
        fit_coef = {'phot_scale': fit_coef[0],
                    'phot_shift': fit_coef[1],
                    'spec_scale': fit_coef[2],
                    'spec_shift': fit_coef[3],
                    'spec_slope': fit_coef[4]}
        if not do_phot_scale:
            fit_coef['phot_scale'] = 1
        if not do_phot_shift:
            fit_coef['phot_shift'] = 0
        if not do_spec_scale:
            fit_coef['spec_scale'] = 1
        if not do_spec_shift:
            fit_coef['spec_shift'] = 0
        if not do_spec_slope:
            fit_coef['spec_slope'] = 0
        scaled_phot, scaled_spec = get_trans_spec(phot, spec, fit_coef)
        diffs = np.interp(scaled_phot, cal_phot, cal_spec)-scaled_spec
        return diffs

    start_params = (phot_scale_guess, phot_shift_guess, 
                    spec_scale_guess, spec_shift_guess, 
                    spec_slope_guess)
    fit_coef, pcov = leastsq(trans_spec_diffs, start_params)
    
    fit_coef_dict = {'phot_scale': fit_coef[0],
                     'phot_shift': fit_coef[1],
                     'spec_scale': fit_coef[2],
                     'spec_shift': fit_coef[3],
                     'spec_slope': fit_coef[4]}

    # Make a plot to see how well the fitting worked:
    f, axs = plt.subplots(2, 1, sharex=True)
    axs[0].plot(cal_phot, cal_spec)
    scaled_phot, scaled_spec = get_trans_spec(phot, spec, fit_coef_dict)
    axs[0].plot(scaled_phot, scaled_spec)
    #axs[1].plot(scaled_phot, spec-scaled_spec)
    axs[1].set_ylabel('Difference')
    axs[1].plot(scaled_phot, scaled_spec-np.interp(scaled_phot, cal_phot, cal_spec))
    return fit_coef_dict

def get_trans_spec(phot, spec, fit_coef):
    transformed_phot = phot*fit_coef['phot_scale']-fit_coef['phot_shift']
    transformed_spec = spec*fit_coef['spec_scale']-fit_coef['spec_shift']
    transformed_spec = transformed_spec-fit_coef['spec_slope']*transformed_phot
    return transformed_phot, transformed_spec

#def get_xas_scan(data, bins=100, equal_norms=True):
#    if equal_norms is True:
#        xas_scan = norm_data.norm_trace_equal_norms(data['delay'], data['andor'], data['norm'], bins=bins)
#    else:
