"""Return processed LCLS data from LK30 experiment which can be plotted easily
"""

import numpy as np
import os
import pickle
from scipy.optimize import curve_fit

import als_process2019
import abs_ana
import calibrate_signals
import calc_cascade_response

#### Get top point fluence cutoff percentile for the data sets
signal_calibrations = calibrate_signals.get_calibrations()
HIGH_CUTOFF = signal_calibrations['fluence_cutoff']['high_fluence']
LOW_CUTOFF = signal_calibrations['fluence_cutoff']['low_fluence']
LOWEST_CUTOFF = signal_calibrations['fluence_cutoff']['lowest_fluence']
####

#### Bottom point fluence cutoff percentile for using data with better stats
CUTOFF_BOTTOM = 30
####

#### Define runs containing different fluence data sets
HIGH_FLUENCE_RUNS = [109, 111, 117, 118]
LOW_FLUENCE_RUNS = [126, 127]
LOWEST_FLUENCE_RUNS = [124, 125]
####

#### Names of files for saving
INCIDENT_SPEC_FILE = os.path.join(os.path.dirname(__file__), '../data/processed/summ_specs.p')
CONSTANT_ABSORPTION_FILE = os.path.join(os.path.dirname(__file__), '../data/processed/constant_absorption_specs.p')
SPEC_QUANTIFICATION_FILE = os.path.join(os.path.dirname(__file__), '../data/processed/quant.p')
DAMAGE_SPEC_FILE = os.path.join(os.path.dirname(__file__), '../data/processed/damage_specs.p')
####

#### Define photon energy binning for calculating spectra
PHOT_BIN_CENTERS = np.linspace(773.5, 782, 59)
dphot_bin = PHOT_BIN_CENTERS[1]-PHOT_BIN_CENTERS[0]
PHOT_BIN_EDGES = PHOT_BIN_CENTERS-dphot_bin/2
PHOT_BIN_EDGES = np.append(PHOT_BIN_EDGES, [PHOT_BIN_EDGES[-1]+dphot_bin])
####

#### Calculate conversion factor of spectral intensity to holes
m_s = 1.96      # [\mu_B/atom], from Loic's sum rule analysis
l3_start = 764
l3_end = 790
l2_start = 790
l2_end = 810
ALS_SPEC = als_process2019.get_als_spectra()
l3_in_range = (ALS_SPEC['xmcd']['phot'] > l3_start) & (ALS_SPEC['xmcd']['phot'] < l3_end)
l2_in_range = (ALS_SPEC['xmcd']['phot'] > l2_start) & (ALS_SPEC['xmcd']['phot'] < l2_end)
l3_xmcd = np.trapz(ALS_SPEC['xmcd']['spec'][l3_in_range], x=ALS_SPEC['xmcd']['phot'][l3_in_range])
l2_xmcd = np.trapz(ALS_SPEC['xmcd']['spec'][l2_in_range], x=ALS_SPEC['xmcd']['phot'][l2_in_range])
C = (-1*l3_xmcd+2*l2_xmcd)/m_s
D = C/1.3  # Conversion factor of spectral intensity to holes for L3 edge
####

# Spectral ranges of interest
FERMI_ENERGY = 777.5
HOLES_RANGE = (FERMI_ENERGY-2, FERMI_ENERGY)
ELECS_RANGE = (FERMI_ENERGY, FERMI_ENERGY+2)
ENERGIES_RANGE = (FERMI_ENERGY-2, FERMI_ENERGY+2)
XMCD_RANGE = (777.56, 778.77)     # FWHM of XMCD
####

# conversion factor from absorbed fluence to absorbed energy density
FLUENCE2ABSORBEDENERGY = 17.3

def get_fluence_binned_data(data, fluence_bin_edges, absorbed=False):
    if absorbed is True:
        fluence_key = 'absorbed_fluence'
    else:
        fluence_key = 'fluence_in'
    binned_data = []
    for bin_edges in fluence_bin_edges:
        fluence_high_enough = data[fluence_key] > bin_edges[0]
        fluence_low_enough = data[fluence_key] < bin_edges[1]
        fluence_in_range = fluence_high_enough & fluence_low_enough
        binned_data.append(data[fluence_in_range])
    return binned_data

def get_incident_specs():
    """Return spectra binned with constant incident fluence
    """
    pickle_out = open(INCIDENT_SPEC_FILE, 'rb')
    summ_specs = pickle.load(pickle_out)
    pickle_out.close()
    return summ_specs

def save_incident_specs():
    # Fluence bins chosen to have similar statistics between data sets
    bin_edges_lowest = [(-0.399, 1)]
    bin_edges_low = [(-1, 100)]
    bin_edges_high = [(25, 65), (65, 160)]
    data_sets = []
    data_lowest_fluence = abs_ana.get_data(LOWEST_FLUENCE_RUNS, True, (CUTOFF_BOTTOM, LOWEST_CUTOFF))
    data_sets.extend(get_fluence_binned_data(data_lowest_fluence, bin_edges_lowest))
    # use higher fluences of set by having low cutoff at 50 percent rather than 30 percent:
    data_low_fluence = abs_ana.get_data(LOW_FLUENCE_RUNS, None, (CUTOFF_BOTTOM, LOW_CUTOFF))
    data_sets.extend(get_fluence_binned_data(data_low_fluence, bin_edges_low))
    data_high_fluence = abs_ana.get_data(HIGH_FLUENCE_RUNS, True, (CUTOFF_BOTTOM, HIGH_CUTOFF))
    data_sets.extend(get_fluence_binned_data(data_high_fluence, bin_edges_high))
    spectra = []
    for ind, data_set in enumerate(data_sets):
        spec = abs_ana.get_spectra(data_set, bins=PHOT_BIN_EDGES)
        spectra.append(spec)
    pickle_on = open(INCIDENT_SPEC_FILE, 'wb')
    pickle.dump(spectra, pickle_on)
    pickle_on.close()

def get_constant_absorbed_specs():
    """Return (approximately) constant absorbed fluence spectra wrt photon energy
    """
    pickle_out = open(CONSTANT_ABSORPTION_FILE, 'rb')
    const_specs = pickle.load(pickle_out)
    pickle_out.close()
    return const_specs

def save_constant_absorbed_specs():
    # Fluence bins chosen to have similar statistics between data sets
    bin_edges_lowest = [(0, 0.0329)]
    bin_edges_low = [(0, 2.743), (2.743, 11.830)]
    bin_edges_high = [(15, 25), (25, 95)]
    data_sets = []
    data_lowest_fluence = abs_ana.get_data(LOWEST_FLUENCE_RUNS, True, (CUTOFF_BOTTOM, LOWEST_CUTOFF))
    data_sets.extend(get_fluence_binned_data(data_lowest_fluence, bin_edges_lowest, absorbed=True))
    # use higher fluences of set by having low cutoff at 50 percent rather than 30 percent:
    # (when run with 30th percentile lower bound instead, produces same result
    #  [no changes from ALS] within experimental error)
    data_low_fluence = abs_ana.get_data(LOW_FLUENCE_RUNS, True, (CUTOFF_BOTTOM, LOW_CUTOFF))
    data_sets.extend(get_fluence_binned_data(data_low_fluence, bin_edges_low, absorbed=True))
    data_high_fluence = abs_ana.get_data(HIGH_FLUENCE_RUNS, True, (CUTOFF_BOTTOM, HIGH_CUTOFF))
    data_sets.extend(get_fluence_binned_data(data_high_fluence, bin_edges_high, absorbed=True))
    spectra = []
    for data_set in data_sets:
        spec = abs_ana.get_spectra(data_set, bins=PHOT_BIN_EDGES)
        spectra.append(spec)
    pickle_on = open(CONSTANT_ABSORPTION_FILE, 'wb')
    pickle.dump(spectra, pickle_on)
    pickle_on.close()
    
def save_damage_specs():
    fluence_set = [(0, 25), (25, 35), (35, 45), (45, 75), (75, 300)]
    spectra = calculate_damage_specs(fluence_set)
    pickle_on = open(DAMAGE_SPEC_FILE, 'wb')
    pickle.dump(spectra, pickle_on)
    pickle_on.close()
    
def calculate_damage_specs(fluence_set):
    max_prev_fluence_key = 'max_prev_fluence'
    data = abs_ana.get_data(HIGH_FLUENCE_RUNS, sample_filter=None, mcp_filter=(CUTOFF_BOTTOM, HIGH_CUTOFF))
    spectra = []
    for fluences in fluence_set:
        data_set = data[data[max_prev_fluence_key] >= fluences[0]]
        data_set = data_set[data_set[max_prev_fluence_key] < fluences[1]]
        spec = abs_ana.get_spectra(data_set, bins=PHOT_BIN_EDGES)
        spectra.append(spec)
    return spectra
    
def get_damage_specs():
    pickle_out = open(DAMAGE_SPEC_FILE, 'rb')
    damage_specs = pickle.load(pickle_out)
    pickle_out.close()
    return damage_specs

def save_spec_quantification():
    spec_quantification = calculate_spec_quantification()
    pickle_on = open(SPEC_QUANTIFICATION_FILE, 'wb')
    pickle.dump(spec_quantification, pickle_on)
    pickle_on.close()

def _integrate_range(x, y, integration_range, steps=500):
    """Integrate data over specificed range
    """
    x_interped = np.linspace(integration_range[0], integration_range[1], steps)
    y_interped = np.interp(x_interped, x, y)
    integral = np.trapz(y_interped, x=x_interped)
    return integral

def calculate_spec_quantification():
    als_spec = als_process2019.get_als_spectra()
    lcls_specs = get_constant_absorbed_specs()
    xmcd_norm = _integrate_range(als_spec['xmcd']['phot'], als_spec['xmcd']['spec'], XMCD_RANGE)
    energies = []
    err_energies = []
    xmcds = []
    err_xmcds = []
    fluences = []
    for spec in lcls_specs:
        holes_per_ev = spec['xas']['diff']/D
        energy_per_ev = holes_per_ev*(FERMI_ENERGY-spec['xas']['phot'])
        energy_current = _integrate_range(spec['xas']['phot'], energy_per_ev, ENERGIES_RANGE)
        xmcd_current = _integrate_range(spec['xmcd']['phot'], spec['xmcd']['spec'], XMCD_RANGE)/xmcd_norm
        variance_energy_per_ev = ((spec['plus']['err']/2)**2+(spec['minus']['err']/2)**2)*(FERMI_ENERGY-spec['xas']['phot'])**2
        variance_energy_current = _integrate_range(spec['xas']['phot'], variance_energy_per_ev, ENERGIES_RANGE)
        variance_xmcd_per_ev = (spec['plus']['err']**2+spec['minus']['err']**2)/(xmcd_norm**2)
        variance_xmcd_current = _integrate_range(spec['xmcd']['phot'], variance_xmcd_per_ev, XMCD_RANGE)
        err_energy_current = np.sqrt(variance_energy_current)
        err_xmcd_current = np.sqrt(variance_xmcd_current)
        energies.append(energy_current)
        err_energies.append(err_energy_current)
        xmcds.append(xmcd_current)
        err_xmcds.append(err_xmcd_current)
        fluences.append(spec['absorbed_fluence2'])
    fluences = np.array(fluences)
    absorbed_energies = fluences*FLUENCE2ABSORBEDENERGY
    eV_to_meV = 1000
    spec_energies = np.array(energies)*eV_to_meV
    line_zero_offset = lambda x, slope: x*slope
    fraction_within_2eV, unused = curve_fit(line_zero_offset, absorbed_energies/2, spec_energies, p0=0.8)
    cascade_duration = calc_cascade_response.calculate_cascade_time(fraction_within_2eV)
    return {'absorbed_fluences': fluences,
            'absorbed_energies': absorbed_energies,
            'pulse_averaged_absorbed_energies': absorbed_energies/2,
            'spec_energies': spec_energies,
            'err_spec_energies': np.array(err_energies)*eV_to_meV,
            'xmcds': np.array(xmcds),
            'err_xmcds': np.array(err_xmcds),
            'fraction_within_2eV': fraction_within_2eV,
            'cascade_duration': cascade_duration}

def get_spec_quantification():
    pickle_out = open(SPEC_QUANTIFICATION_FILE, 'rb')
    damage_specs = pickle.load(pickle_out)
    pickle_out.close()
    return damage_specs
