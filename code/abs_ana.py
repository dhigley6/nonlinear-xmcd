"""Processing of 2016 LCLS nonlinear XMCD experiment

01/2016
Daniel Higley
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats

import norm_data
import abs_hdf2rec
Pd_ABS = 0.315    # Calculated from CXRO
import als_process2019
ALS_SPEC = als_process2019.get_als_spectra()

# fluence values obtained from abs_hdf2rec are pulse energy/(FWHM area).
# The below factor can be used to convert these to average fluence values
# (one half of the peak fluence)
FWHM_FLUENCE2AVERAGE_FLUENCE = 0.4412712003053032

def get_spectra(data, bins=100):
    """Return XMCD data for input data set
    """

    def calc_diff(spec, ref_spec):
        """Calculate the difference between input spectrum and reference spectrum
        """
        ref_interped = np.interp(spec['phot'], ref_spec['phot'], ref_spec['spec'])
        diff_spec = spec['spec']-ref_interped
        out_of_range = (spec['phot'] > np.amax(ref_spec['phot'])) | (spec['phot'] < np.amin(ref_spec['phot']))
        diff_spec[out_of_range] = np.nan
        return diff_spec

    data_p = data[data['magnet_dir'] == 1]
    data_m = data[data['magnet_dir'] == -1]
    spec_p = norm_data.norm_trace(data_p['phot'], data_p['fluence_out'], data_p['fluence_in'], bins=bins)
    spec_p['spec'] = -1*np.log(spec_p['norm_sig'])
    spec_m = norm_data.norm_trace(data_m['phot'], data_m['fluence_out'], data_m['fluence_in'], bins=bins)
    spec_m['spec'] = -1*np.log(spec_m['norm_sig'])
    spec_p['fluence_in'] = spec_p['norm_sums']/spec_p['bin_counts']
    spec_m['fluence_in'] = spec_m['norm_sums']/spec_m['bin_counts']
    abs_frac_plus_als_phot = 1-np.exp(-1*ALS_SPEC['plus']['spec'])
    abs_frac_plus = np.interp(spec_p['bin_centers'], ALS_SPEC['plus']['phot'], abs_frac_plus_als_phot)
    spec_p['abs_fluence'] = spec_p['fluence_in']*abs_frac_plus
    abs_frac_minus_als_phot = 1-np.exp(-1*ALS_SPEC['minus']['spec'])
    abs_frac_minus = np.interp(spec_m['bin_centers'], ALS_SPEC['minus']['phot'], abs_frac_minus_als_phot)
    spec_m['abs_fluence'] = spec_m['fluence_in']*abs_frac_minus

    specs = {}
    specs['plus'] = spec_p
    specs['minus'] = spec_m
    xas = (specs['minus']['spec']+specs['plus']['spec'])/2
    specs['xas'] = {'bin_centers': spec_m['bin_centers'],
                    'bin_edges': spec_m['bin_edges'],
                    'spec': xas}
    xmcd = (specs['plus']['spec']-specs['minus']['spec'])
    specs['xmcd'] = {'bin_centers': spec_m['bin_centers'],
                     'bin_edges': spec_m['bin_edges'],
                     'spec': xmcd}
    for key in specs:
        if key != 'xmcd':
            specs[key]['spec'] = specs[key]['spec']+Pd_ABS
        specs[key]['phot'] = specs[key]['bin_centers']
        specs[key]['diff'] = calc_diff(specs[key], ALS_SPEC[key])
    specs['fluence'] = np.sum(data['fluence_in']**2)/np.sum(data['fluence_in'])
    specs['absorbed_fluence'] = np.sum(data['absorbed_fluence']*data['fluence_in'])/np.sum(data['fluence_in'])
    specs['absorbed_fluence2'] = np.mean(spec_p['abs_fluence']+spec_m['abs_fluence'])/2
    specs['fluence'] = specs['fluence']*FWHM_FLUENCE2AVERAGE_FLUENCE
    specs['absorbed_fluence'] = specs['absorbed_fluence']*FWHM_FLUENCE2AVERAGE_FLUENCE
    specs['absorbed_fluence2'] = specs['absorbed_fluence2']*FWHM_FLUENCE2AVERAGE_FLUENCE
    return specs

def get_data(runs, sample_filter=True, mcp_filter=None):
    """Return data filtered by specifications
    """
    data = abs_hdf2rec.get_rec_array(runs)
    if sample_filter is not None:
        data = data[data['intact_sample'] == sample_filter]
    if mcp_filter is not None:
        mcp_low_perc = mcp_filter[0]
        mcp_high_perc = mcp_filter[1]
        mcp_low_val = scipy.stats.scoreatpercentile(data['mcp'], mcp_low_perc)
        mcp_high_val = scipy.stats.scoreatpercentile(data['mcp'], mcp_high_perc)
        mcp_high_enough = data['mcp'] >= mcp_low_val
        mcp_low_enough = data['mcp'] < mcp_high_val
        data = data[(mcp_high_enough) & (mcp_low_enough)]
    return data