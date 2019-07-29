"""Scripts for figuring out the calibration of LK30 spectra
"""

import numpy as np
import matplotlib.pyplot as plt

import abs_ana
import abs_hdf2rec
import als_process2019
import run_cals
import abs_get_processed_data
import xas

import calibrate_signals
#### Get top point fluence cutoff percentile for the data sets
signal_calibrations = calibrate_signals.get_calibrations()
HIGH_CUTOFF = signal_calibrations['fluence_cutoff']['high_fluence']
LOW_CUTOFF = signal_calibrations['fluence_cutoff']['low_fluence']
LOWEST_CUTOFF = signal_calibrations['fluence_cutoff']['lowest_fluence']
####

#### Bottom point fluence cutoff percentile for using data with better stats
CUTOFF_BOTTOM = 30
####

ALS_SPEC = als_process2019.get_als_spectra()

def do_cals():
    """Perform all calibrations
    """
    cal_coef = cal_to_whole_spectrum([124, 125], phot_max=782, mcp_filter=(CUTOFF_BOTTOM, LOWEST_CUTOFF))
    runs_to_change = [105, 108, 109, 120, 123, 124, 126, 129]
    for run in runs_to_change:
        run_cals.change_run_pars(run, cal_coef)
    cal_coef = cal_away_from_edge([126, 127], mcp_filter=(CUTOFF_BOTTOM, LOW_CUTOFF))
    run_cals.change_run_pars(126, cal_coef)
    #cal_coef = cal_away_from_edge([109, 111, 117, 118], sample_filter=True, mcp_filter=(CUTOFF_BOTTOM, 99.9), fluence_filter=(30, 65))
    cal_coef = cal_away_from_edge([109, 111, 117, 118], sample_filter=True, mcp_filter=(CUTOFF_BOTTOM, 99.9), fluence_filter=(30, 65))
    run_cals.change_run_pars(108, cal_coef)
    run_cals.change_run_pars(109, cal_coef)
    #cal_coef = cal_to_whole_spectrum([105], phot_max=782)
    #run_cals.change_run_pars(105, cal_coef)

def cal_to_whole_spectrum(runs=[126, 127], phot_max=10000, mcp_filter=(CUTOFF_BOTTOM, 95)):
    lcls_xas = get_lcls_xas(runs, phot_max, mcp_filter)
    run_pars = abs_hdf2rec.get_run_params(runs[0])
    als = ALS_SPEC['xas']
    print(''.join(['\nFitting runs ', str(runs), ' to ALS:']))
    cal_coef = cal_spec(lcls_xas, als, run_pars, spec_scale_guess=1)
    #run_cals.change_run_pars(runs[0], cal_coef)
    plt.title('Fitting for runs '+str(runs)+' to ALS')
    return cal_coef

def cal_away_from_edge(runs=[126, 127], sample_filter=None, mcp_filter=(CUTOFF_BOTTOM, 95), fluence_filter=None):
    lcls_xas = get_lcls_xas(runs, sample_filter=sample_filter, mcp_filter=mcp_filter, fluence_filter=fluence_filter)
    phot = lcls_xas['phot']
    phot_below_edge = (phot < 773.5) & (phot < 774.5)
    phot_above_edge = (phot > 781) & (phot < 782)
    phot_in_range = phot_below_edge | phot_above_edge
    lcls_xas['spec'] = lcls_xas['spec'][phot_in_range]
    lcls_xas['phot'] = lcls_xas['phot'][phot_in_range]
    run_pars = abs_hdf2rec.get_run_params(runs[0])
    als = ALS_SPEC['xas']
    print(''.join(['\nFitting runs ', str(runs), ' to ALS:']))
    cal_coef = cal_spec(lcls_xas, als, run_pars, phot_scale_guess=None, phot_shift_guess=None, spec_scale_guess=None)
    #run_cals.change_run_pars(runs[0], cal_coef)
    plt.title('Fitting for runs '+str(runs)+' to ALS')
    return cal_coef
    

def get_lcls_xas(runs, phot_max=10000, mcp_filter=(CUTOFF_BOTTOM, 95), sample_filter=None, fluence_filter=None):
    lcls_data = abs_ana.get_data(runs, sample_filter, mcp_filter=mcp_filter)
    if fluence_filter is not None:
        lcls_data = lcls_data[lcls_data['fluence_in'] > fluence_filter[0]]
        lcls_data = lcls_data[lcls_data['fluence_in'] < fluence_filter[1]]
    lcls_spec = abs_ana.get_spectra(lcls_data, bins=50)
    lcls_xas = lcls_spec['xas']
    # Take out only LCLS photon energies where we know there
    # are no artifacts
    good = (lcls_xas['phot'] < phot_max) & (~np.isnan(lcls_xas['spec']))
    lcls_xas['spec'] = lcls_xas['spec'][good]
    lcls_xas['phot'] = lcls_xas['phot'][good]
    return lcls_xas

def cal_spec(spec, fit_spec, run_pars, **kwargs):
    fit_coef = xas.fit_abs(spec['phot'], spec['spec'], fit_spec['phot'], fit_spec['spec'], **kwargs)
    # Print fit coefficients:
    print(''.join(['Photon energy scaling: ', str(fit_coef['phot_scale'])]))
    print(''.join(['Photon energy shift: ', str(fit_coef['phot_shift'])]))
    print(''.join(['Spectrum scaling: ', str(fit_coef['spec_scale'])]))
    print(''.join(['Spectrum shift: ', str(fit_coef['spec_shift'])]))
    print(''.join(['Spectrum slope: ', str(fit_coef['spec_slope'])]))
    # Calculate how to change calibration coefficients based on fit
    # coefficients
    cal_coef = {'phot_scale': fit_coef['phot_scale']*run_pars['phot_scale'],
                'phot_shift': run_pars['phot_shift']*fit_coef['phot_scale']+fit_coef['phot_shift'],
                'norm_offset': run_pars['norm_offset']+fit_coef['spec_shift']+fit_coef['spec_slope']*778,
                'norm_slope': run_pars['norm_slope']+fit_coef['spec_slope']}
    # Print calibration coefficients
    print('***Calculated Calibration Coefficients***')
    print(''.join(['phot_scale ', str(cal_coef['phot_scale'])]))
    print(''.join(['phot_shift ', str(cal_coef['phot_shift'])]))
    print(''.join(['norm_offset ', str(cal_coef['norm_offset'])]))
    print(''.join(['norm_slope ', str(cal_coef['norm_slope'])]))
    return cal_coef


