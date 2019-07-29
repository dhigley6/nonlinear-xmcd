"""Calibrate signals recorded in LK30 experiment
"""

import numpy as np
import matplotlib.pyplot as plt
import h5py
import pickle
import os
from scipy.stats import percentileofscore
from scipy.stats import scoreatpercentile

import run_cals

#### Define runs containing different fluence data sets
HIGH_FLUENCE_RUNS = [109, 111, 117, 118]
LOW_FLUENCE_RUNS = [126, 127]
LOWEST_FLUENCE_RUNS = [124, 125]
####

# Path to the file with dark (no X-ray) data
DARK_DATA_PATH = os.path.join(os.path.dirname(__file__), '../data/preprocessed/run112allevts.h5')
# Path to data directory
DATA_PATH = os.path.join(os.path.dirname(__file__), '../data/preprocessed/')
# File for dark calibrations:
DARK_CALIBRATION_FILE = os.path.join(os.path.dirname(__file__), '../cals/dark_calibrations.p')
# File for andor saturation calibrations:
SATURATION_CALIBRATION_FILE = os.path.join(os.path.dirname(__file__), '../cals/saturation_calibrations.p')
# File for MCP nonlinearity calibrations:
MCP_CALIBRATION_FILE = os.path.join(os.path.dirname(__file__), '../cals/mcp_calibrations.p')

def do_calibrations():
    do_dark_calibration()
    do_andor_saturation_calibration()
    do_mcp_nonlinearity_calibration()

def get_calibrations():
    dark_calibration = _get_dark_calibration()
    andor_calibration = _get_andor_calibration()
    mcp_calibration = _get_mcp_nonlinearity_calibration()
    return {'dark': dark_calibration,
            'fluence_cutoff': andor_calibration,
            'mcp': mcp_calibration}

def _get_dark_calibration():
    pickle_out = open(DARK_CALIBRATION_FILE, 'rb')
    dark_calibrations = pickle.load(pickle_out)
    pickle_out.close()
    return dark_calibrations
    
def _get_andor_calibration():
    pickle_out = open(SATURATION_CALIBRATION_FILE, 'rb')
    andor_calibrations = pickle.load(pickle_out)
    pickle_out.close()
    return andor_calibrations
    
def _get_mcp_nonlinearity_calibration():
    pickle_out = open(MCP_CALIBRATION_FILE, 'rb')
    mcp_calibrations = pickle.load(pickle_out)
    pickle_out.close()
    return mcp_calibrations

def do_dark_calibration():
    f = h5py.File(DARK_DATA_PATH, 'r')
    dark_andor_mean = np.mean(np.array(f['Andor']['signal'])-np.array(f['Andor']['reference']))
    dark_mcp_mean = np.mean(f['Acqiris2']['acq'])
    dark_cals = {'andor': dark_andor_mean,
                 'mcp': dark_mcp_mean}
    pickle_on = open(DARK_CALIBRATION_FILE, 'wb')
    pickle.dump(dark_cals, pickle_on)
    pickle_on.close()
    
def do_andor_saturation_calibration():
    cutoff_percentile_lowest_fluence = _calculate_percentile_cutoff(LOWEST_FLUENCE_RUNS)
    cutoff_percentile_low_fluence = _calculate_percentile_cutoff(LOW_FLUENCE_RUNS)
    cutoff_percentile_high_fluence = _calculate_percentile_cutoff(HIGH_FLUENCE_RUNS)
    saturation_cals = {'lowest_fluence': cutoff_percentile_lowest_fluence,
                       'low_fluence': cutoff_percentile_low_fluence,
                       'high_fluence': cutoff_percentile_high_fluence}
    pickle_on = open(SATURATION_CALIBRATION_FILE, 'wb')
    pickle.dump(saturation_cals, pickle_on)
    pickle_on.close()

def _calculate_percentile_cutoff(run_numbers):
    """Calculate upper percentile MCP cutoff for avoiding Andor saturation
    Andor begins saturating at ~5000 for 'signal' value.
    Set percentile cutoff of incident fluence (mcp) to the percentile where
    Andor reaches 4000, well before saturation. If that is greater than the
    99.9th percentile, set percentile cutoff to 99.9 to eliminate potential
    strong outliers in incident fluence.
    
    Note: photon energy calibrations should already be approximately correct
    before running this function.
    """
    mcp_values = []
    andor_values = []
    for run_number in run_numbers:
        current_data_path = ''.join([DATA_PATH, 'run', str(run_number), 'allevts.h5'])
        f = h5py.File(current_data_path, 'r')
        current_phot = _get_photon_energy(f, run_number)
        current_mcp = np.array(f['Acqiris2']['acq'])
        current_mcp = current_mcp[(current_phot > 781) & (current_phot < 782)]
        mcp_values.extend(current_mcp)
        current_andor = np.array(f['Andor']['signal'])
        current_andor = current_andor[(current_phot > 781) & (current_phot < 782)]
        andor_values.extend(current_andor)
    #plt.figure()
    #plt.scatter(mcp_values, andor_values)
    mcp_percentile_cutoff = min([percentileofscore(andor_values, 4000), 99.9])
    return mcp_percentile_cutoff
    
def do_mcp_nonlinearity_calibration():
    """Calculate polynomial to correct slight nonlinearity in MCP response
    Note: dark calibration should be done before this
    """
    no_sample_data_path = ''.join([DATA_PATH, 'run108allevts.h5'])
    f = h5py.File(no_sample_data_path)
    phot = _get_photon_energy(f, 108)
    mcp = np.array(f['Acqiris2']['acq'])
    andor = np.array(f['Andor']['signal']-f['Andor']['reference'])
    # Subtract dark signals:
    dark_calibration = _get_dark_calibration()
    mcp = mcp-dark_calibration['mcp']
    andor = andor-dark_calibration['andor']
    # Take data within (relatively) narrow photon energy range:
    phot_in_range = (phot > 781) & (phot < 782)
    mcp = mcp[phot_in_range]
    andor = andor[phot_in_range]
    # make sure to only take data for which andor doesn't saturate
    mcp_percentile_cutoff = min([percentileofscore(andor, 4000), 99.9])
    mcp_cutoff_value = scoreatpercentile(mcp, mcp_percentile_cutoff)
    mcp_in_range = mcp < mcp_cutoff_value
    mcp = mcp[mcp_in_range]
    andor = andor[mcp_in_range]
    correction_polynomial = np.polyfit(
        mcp, 
        andor*(np.mean(mcp)/np.mean(andor)),
        deg=3)
    plt.figure()
    plt.scatter(mcp, andor)
    plt.scatter(np.polyval(correction_polynomial, mcp), andor)
    pickle_on = open(MCP_CALIBRATION_FILE, 'wb')
    pickle.dump(correction_polynomial, pickle_on)
    pickle_on.close()
    
def _get_photon_energy(file, run_number):
    monoencoder = np.array(file['monoenc']['encoder_count'][:, 0])
    run_pars = run_cals.get_run_pars(run_number)
    phot = monoencoder*run_pars['phot_scale']
    phot = phot-run_pars['phot_shift']
    return phot