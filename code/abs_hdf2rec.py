"""Convert data from HDF5 file to numpy record array for 2016 LCLS nonlinear XAS

2016-04
Daniel Higley
"""

# Standard python modules:
import numpy as np
import h5py
import os

# Custom python modules:
import run_cals
import calibrate_signals
signal_calibrations = calibrate_signals.get_calibrations()
import als_process2019

#### CALIBRATIONS #############################################################
# The threshold of magnet voltage at which to set a filter
MAG_VOLT_THRESH = 0.13
# Path to the directory with data
DATA_PATH = os.path.join(os.path.dirname(__file__), '../data/preprocessed/')
# Data with an absorbed fluence larger than the below specified value 
# are designated as having a potentially damaged sample
FLUENCE_CUTOFF_ABS = 22.5
# Dark counts taken from run 112:
MCP_DARK = signal_calibrations['dark']['mcp']
ANDOR_DARK = signal_calibrations['dark']['andor']
# Polynomial coefficients for correcting slight mcp nonlinearity:
MCP_CORRECTION_POLYNOMIAL = signal_calibrations['mcp']
###############################################################################

def get_rec_array(runs):
    """Return data from runs as a record array
    """
    data_list = []
    for run in runs:
        # Load run specific parameters
        run_spec_pars = get_run_params(run)
        # Get data file path
        curr_data_path = ''.join([DATA_PATH, 'run', str(run),
                                  'allevts.h5'])
        run_hdf_file = h5py.File(curr_data_path, 'r')
        data_dict = {}
        run_dict = {'run_num': [run]*len(run_hdf_file['gd']['f_21_ENRC'])}
        data_dict.update(run_dict)
        data_dict.update(get_signals(run_hdf_file, run_spec_pars['burst']))
        # subtract mcp and andor dark counts taken from run 112:
        data_dict['mcp'] = data_dict['mcp']-MCP_DARK
        data_dict['andor'] = data_dict['andor']-ANDOR_DARK
        # Correct mcp nonlinearity:
        data_dict['mcp_uncorrected'] = data_dict['mcp'].copy()
        data_dict['mcp'] = np.polyval(MCP_CORRECTION_POLYNOMIAL, data_dict['mcp'])
        data_dict.update(get_mono(run_hdf_file))
        data_dict.update(get_magnet(run_hdf_file))
        data_dict.update(get_manip(run_hdf_file))
        data_dict.update(get_scaled_pars(data_dict, run_spec_pars))
        data_dict.update(get_absorbed_fluence(data_dict))
        data_dict.update(get_sample_info(run_hdf_file, run_spec_pars['burst'],
                                         data_dict['fluence_in'],
                                         data_dict['absorbed_fluence'],
                                         data_dict['magnet_dir']))
        data_dict.update({'ind': np.arange(len(data_dict['fluence_in']))})
        # Convert data from dictionary to a numpy record array:
        data = np.rec.fromarrays(data_dict.values(), names=[*data_dict.keys()])
        # Apply run-specific filter to filter out data that wasn't recorded in correct conditions
        data = run_spec_filter(data, run)
        # Append the current run's data to the list of run data
        data_list.append(data)
    data = np.hstack(data_list)
    return data

def get_run_params(run):
    """Return the parameters a run was taken with and extracted parameters
    """
    run_params = run_cals.get_run_pars(run)
    return run_params

def run_spec_filter(data, run):
    """Filter out data that was not recorded with the desired experimental conditions
    """
    if (run == 109) or (run == 111):
        # Alignment is not good for low sample x values on the
        # chip scanned in these runs -> take these out
        data = data[data['sam_x'] > 25]
    if (run == 111):
        # The very last sample in this run was blown up by X-rays before being
        #   scanned over. Take this out of analysis
        data = data[data['ind'] < 30304]
    if (run == 118):
        # DAQ lost synchronization at the end of the run -> remove
        #   unsynchronized data from the analysis
        data = data[:172000]
    if (run == 124):
        # DAQ lost synchronization at the end of the run -> remove
        #    unsynchronized data from the analysis
        data = data[:51000]
    return data

def get_signals(f, burst_mode):
    """Extract X-ray signals from the input file
    """
    # Create dictionary to hold signals:
    sigs = {}
    # Get the gas detector data after the gas attenuator
    gd21 = f['gd']['f_21_ENRC']
    gd22 = f['gd']['f_22_ENRC']
    sigs['gd2'] = (gd21+gd22)/2
    # Get the MCP values (I0)
    acq2 = np.array(f['Acqiris2']['acq'])
    sigs['mcp'] = np.array(acq2)
    # Get the Andor signals (transmitted X-rays)
    andor = np.array(f['Andor']['signal'])
    andor_ref = np.array(f['Andor']['reference'])
    sigs['andor'] = andor-andor_ref
    if burst_mode:
        # Subtract backgrounds on Andor and mcp
        beam_on = f['EVR']['150'].astype(bool)
        mcp_bg = np.mean(sigs['mcp'][~beam_on])
        andor_bg = np.mean(sigs['andor'][~beam_on])
        sigs['mcp'] = sigs['mcp']-mcp_bg
        sigs['andor'] = sigs['andor']-andor_bg
    return sigs

def get_mono(f):
    """Extract monochromator position from the input file
    """
    if 'monoenc' in f.keys():
        mono_pos = np.array(f['monoenc']['encoder_count'][:, 0])
    else:
        mono_pos = f['Epics']['SXR:MON:MMS:06.RBV']
    # Make sure the mono position is stored as floats to avoid rounding
    # errors later:
    mono_pos = np.array(mono_pos, dtype=float)
    mono = {'mono_pos': mono_pos}
    return mono

def get_magnet(f):
    """Return magnet readback voltage and applied field direction
    """
    # Get the magnet readback voltage:
    magnet_val = f['Epics']['SXR:EXP:AIN:1']
    # Set magnet dir according to the direction of the applied magnetic
    # field:
    #     -1 for negative applied field
    #     0 for field intermediate between negative and positive
    #     1 for positive applied field
    magnet_dir = np.zeros(len(magnet_val))
    magnet_dir[magnet_val > MAG_VOLT_THRESH] = 1
    magnet_dir[magnet_val < -1*MAG_VOLT_THRESH] = -1
    magnet = {'magnet_val': magnet_val,
              'magnet_dir': magnet_dir}
    return magnet

def get_manip(f):
    """Return sample manipulator parameters
    """
    manip = {}
    manip['sam_x'] = f['Epics']['SXR:EXP:MMS:01.RBV']
    manip['sam_y'] = f['Epics']['SXR:EXP:MMS:02.RBV']
    manip['sam_z'] = f['Epics']['SXR:EXP:MMS:03.RBV']
    sam_x_d = np.append([0], np.diff(manip['sam_x']))
    sam_y_d = np.append([0], np.diff(manip['sam_y']))
    sam_z_d = np.append([0], np.diff(manip['sam_z']))
    sam_change = (np.abs(sam_x_d) > 0.1) | (np.abs(sam_y_d) > 0.1)
    sam_change = sam_change | (np.abs(sam_z_d) > 0.1)
    manip['sam_change'] = sam_change
    sam_change_within_100 = np.convolve(sam_change, np.ones(100), mode='same')
    sam_change_within_100 = np.array(sam_change_within_100, dtype=bool)
    manip['sam_change_within_100'] = sam_change_within_100
    return manip

def get_sample_info(f, burst_mode, fluence, abs_fluence, magnet_dir):
    """Calculate for which shots the sample is intact and other sample info
    """
    if burst_mode:
        # If in burst mode, calculate # of shot on sample for each shot,
        # as well as whether sample is intact from the absorbed fluences it has
        # been exposed to previously:
        burst = f['EVR']['150'].astype(bool)
        beam_on = burst
        shot_num = np.zeros(len(burst))
        max_prev_fluence = np.zeros(len(burst))
        max_prev_abs_fluence = np.zeros(len(burst))
        intact_sample = np.zeros(len(burst), dtype=bool)
        # Current shot number on sample (starting from 1, 0 if out of burst mode):
        curr_shot_num = 0
        curr_max_prev_fluence = 0
        curr_max_prev_abs_fluence = 0
        for i in np.arange(len(burst)):
            if burst[i]:
                curr_shot_num += 1
            else:
                curr_shot_num = 0
            shot_num[i] = curr_shot_num
            if curr_shot_num > 1:
                curr_max_prev_fluence = np.amax([curr_max_prev_fluence, fluence[i-1]])
                curr_max_prev_abs_fluence = np.amax([curr_max_prev_abs_fluence, abs_fluence[i-1]])
            else:
                curr_max_prev_fluence = np.amin(fluence)
                curr_max_prev_abs_fluence = np.amin(abs_fluence)
            max_prev_fluence[i] = curr_max_prev_fluence
            max_prev_abs_fluence[i] = curr_max_prev_abs_fluence
        intact_sample = burst & (max_prev_abs_fluence <= FLUENCE_CUTOFF_ABS)
    else:
        # If not in burst mode, simply say the sample is always intact:
        intact_sample = np.ones(len(fluence), dtype=bool)
        shot_num = np.arange(len(fluence))
        # Set max previous fluence sample has been exposed to to 0
        max_prev_fluence = np.zeros(len(fluence))
        max_prev_abs_fluence = np.zeros(len(fluence))
        beam_on = np.ones_like(intact_sample)
    sample_info = {'beam_on': beam_on,
                   'shot_num': shot_num,
                   'max_prev_fluence': max_prev_fluence,
                   'max_prev_abs_fluence': max_prev_abs_fluence,
                   'intact_sample': intact_sample}
    return sample_info

def get_scaled_pars(data, run_spec_pars):
    """Scale pulse energy and photon energy values to real units
    """
    def scale_to_fluence(to_scale, phot, scale1, scale2, scale2_phot_zero=778):
        """Return input scaled to fluence
        """
        scaled = to_scale*(scale1+scale2*(phot-scale2_phot_zero))
        return scaled
    # Scale the monochromator position to real photon energy units:
    data['phot'] = data['mono_pos']*run_spec_pars['phot_scale']
    data['phot'] = data['phot']-run_spec_pars['phot_shift']
    # Scale MCP and Andor signals to fluences:
    data['fluence_in'] = scale_to_fluence(data['mcp'], data['phot'],
                                          run_spec_pars['mcp2fluence1'],
                                          run_spec_pars['mcp2fluence2'],
                                          778)
    data['fluence_out'] = scale_to_fluence(data['andor'], data['phot'],
                                           run_spec_pars['mcp2fluence1'],
                                           run_spec_pars['mcp2fluence2'],
                                           778)
    #   Convert XAS normalization paramters to transmission:
    norm_scale = np.exp(run_spec_pars['norm_offset'])
    norm_scale_slope = np.exp(run_spec_pars['norm_slope']*(data['phot']-778))
    #   Aplly transmission normalization parameters:
    data['fluence_out'] = data['fluence_out']*norm_scale
    data['fluence_out'] = data['fluence_out']*norm_scale_slope
    return data

def get_absorbed_fluence(data):
    als_spec = als_process2019.get_als_spectra()
    absorbed_fluences = []
    for ind in np.arange(len(data['phot'])):
        if data['magnet_dir'][ind] == 1:
            absorption = np.interp(data['phot'][ind], als_spec['plus']['phot'], als_spec['plus']['spec'])
        if data['magnet_dir'][ind] == -1:
            absorption = np.interp(data['phot'][ind], als_spec['minus']['phot'], als_spec['minus']['spec'])
        else:
            absorption = np.interp(data['phot'][ind], als_spec['xas']['phot'], als_spec['xas']['spec'])
        absorb_fraction = 1-np.exp(-1*absorption)
        absorbed_fluences.append(absorb_fraction*data['fluence_in'][ind])
    return {'absorbed_fluence': np.array(absorbed_fluences)}
