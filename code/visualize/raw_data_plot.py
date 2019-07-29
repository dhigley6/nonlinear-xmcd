"""Make a plot of the raw data for 2019 LK30 manuscript
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
import os

# LCLS psana to read raw data
# (needs to be imported if loading raw data, otherwise can be
#  commented out)
#import psana

import visualize.set_plot_params
visualize.set_plot_params.init_paper_small()

# To get processed data
import abs_hdf2rec
import abs_ana
import abs_get_processed_data

# Get ALS data
import als_process2019
ALS_SPEC = als_process2019.get_als_spectra()

# Short dash linestyle:
SD = (0, (3, 1.5))

RAW_MCP_AND_ANDOR_FILE = '../../data/raw/raw_mcp_and_andor.p'

def get_raw_mcp_and_andor_waveforms():
    """Return raw mcp and andor waveforms that were
    saved from raw data
    """
    pickle_out = open(RAW_MCP_AND_ANDOR_FILE, 'rb')
    raw_mcp_and_andor = pickle.load(pickle_out, encoding='bytes')
    pickle_out.close()
    return raw_mcp_and_andor

def get_raw_mcp_waveform():
    #psana.setConfigFile('psanacfg.cfg')
    ds = psana.DataSource('exp=sxrk3016:run=126:smd')
    det = psana.Detector('Acq02')
    evt = ds.events().next()
    evt = ds.events().next()
    cutoff = 79960      # take off the last bit of the data for plotting, since the time seems to reset to zero there
    mcp_waveform = {'x': det.wftime(evt)[0][:cutoff]*10**6,
                    'y': det.waveform(evt)[0][:cutoff]}
    return mcp_waveform

def get_andor_waveform():
    ds = psana.DataSource('exp=sxrk3016:run=126:smd')
    det = psana.Detector('andor')
    evt = ds.events().next()
    evt = ds.events().next()
    andor_waveform = {'x': np.arange(len(det.image(evt)[0])),
                      'y': det.image(evt)[0]}
    return andor_waveform

def get_no_sample_data():
    no_sample_data = abs_ana.get_data([108], mcp_filter=(0, 100))
    return no_sample_data

def make_plot():
    f, axs = plt.subplots(2, 2, figsize=(4, 4))
    waveforms = get_raw_mcp_and_andor_waveforms()
    mcp_waveform = waveforms[b'mcp']
    andor_waveform = waveforms[b'andor']
    no_sample = get_no_sample_data()
    specs = abs_get_processed_data.get_incident_specs()
    axs[0, 0].plot(mcp_waveform[b'x'], mcp_waveform[b'y'])
    axs[0, 1].plot(andor_waveform[b'x'], andor_waveform[b'y']/1E3)
    num_points = 100
    mcp_norm = np.amax(no_sample['mcp'][:num_points])
    andor_norm = np.amax(no_sample['andor'][:num_points])
    axs[1, 0].scatter(no_sample['mcp_uncorrected'][:num_points]/mcp_norm, no_sample['andor'][:num_points]/andor_norm, s=4, label='Raw')
    axs[1, 0].scatter(no_sample['mcp'][:num_points]/mcp_norm, no_sample['andor'][:num_points]/andor_norm, s=4, label='Corrected')
    lcls_low, = axs[1, 1].plot(specs[-3]['xas']['phot'], specs[-3]['xas']['spec'], color='c', label='2.7 mJ/cm$^2$')
    lcls_high, = axs[1, 1].plot(specs[-1]['xas']['phot'], specs[-1]['xas']['spec'], color='r', label='43 mJ/cm$^2$')
    als, = axs[1, 1].plot(ALS_SPEC['xas']['phot'], ALS_SPEC['xas']['spec'], color='k', linestyle=SD, label='ALS')
    axs[1, 1].legend(handles=[als, lcls_low, lcls_high], loc='lower center', fontsize=7, borderpad=0, labelspacing=0, handlelength=1.5)
    format_plot(axs, mcp_waveform, andor_waveform)
    fig_path = os.path.join(os.path.dirname(__file__), 'figure_s2.eps')
    plt.savefig(fig_path, dpi=600)
    fig_path = os.path.join(os.path.dirname(__file__), 'figure_s2.jpg')
    plt.savefig(fig_path, dpi=600)

def format_plot(axs, mcp_waveform, andor_waveform):
    axs[0, 0].set_xlabel('pass')
    axs[0, 0].axvline(mcp_waveform[b'x'][12800], color='k', linestyle='--')
    axs[0, 0].axvline(mcp_waveform[b'x'][13000], color='k', linestyle='--')
    axs[0, 0].axvline(mcp_waveform[b'x'][1000], color='k', linestyle='--')
    axs[0, 0].axvline(mcp_waveform[b'x'][9000], color='k', linestyle='--')
    axs[0, 1].axvline(andor_waveform[b'x'][400], color='k', linestyle='--')
    axs[0, 1].axvline(andor_waveform[b'x'][1100], color='k', linestyle='--')
    axs[0, 1].axvline(andor_waveform[b'x'][100], color='k', linestyle='--')
    axs[0, 1].axvline(andor_waveform[b'x'][300], color='k', linestyle='--')
    axs[0, 0].set_xlim((1, 1.8))
    axs[0, 1].set_xlim((0, 2048))
    axs[0, 0].set_xlabel('Time ($\mu s$)')
    axs[0, 0].set_ylabel('Voltage (V)')
    axs[0, 1].set_xlabel('Pixel')
    axs[0, 1].set_ylabel('Counts/10$^3$')
    axs[1, 0].set_xlabel('I$_0$ (a.u.)')
    axs[1, 0].set_ylabel('I$_1$ (a.u.)')
    axs[1, 1].set_ylabel('Intensity')
    axs[1, 1].set_xlabel('Photon Energy (eV)')
    axs[1, 1].set_xlim((745, 840))
    axs[1, 1].set_xlim((773.5, 782))
    axs[1, 1].axvline(774.5, color='k', linestyle='--')
    axs[1, 1].axvline(781, color='k', linestyle='--')
    axs[0, 0].annotate('signal', xy=(1.61, -0.25), xytext=(1.5, -0.15),
                       horizontalalignment='right', 
                       arrowprops=dict(facecolor='black', arrowstyle='->'))
    axs[0, 0].annotate('bg', xy=(1.03, -0.45), xytext=(1.3, -0.38),
                       horizontalalignment='left',
                       arrowprops=dict(facecolor='black', arrowstyle='->'))
    axs[0, 1].annotate('signal', xy=(850, 19), xytext=(1200, 23),
                       horizontalalignment='left',
                       arrowprops=dict(facecolor='black', arrowstyle='->'))
    axs[0, 1].annotate('bg', xy=(130, 10), xytext=(1400, 10),
                       horizontalalignment='left',
                       arrowprops=dict(facecolor='black', arrowstyle='->'))
    axs[1, 1].annotate('Cal.\nRegion', xy=(773.75, 0.8), xytext=(774.75, 1),
                       horizontalalignment='left',
                       arrowprops=dict(facecolor='black', arrowstyle='->'))
    axs[1, 1].annotate('Cal.\nRegion', xy=(781.75, 0.35), xytext=(780, 0.4),
                       horizontalalignment='right',
                       arrowprops=dict(facecolor='black', arrowstyle='->'))
    axs[1, 0].legend(loc='upper left', fontsize=7)
    axs[1, 1].set_ylim((-0.2, 1.25))
    axs[0, 0].text(0.9, 0.1, 'A', fontsize=10, weight='bold', horizontalalignment='center', transform=axs[0, 0].transAxes)
    axs[0, 1].text(0.9, 0.1, 'B', fontsize=10, weight='bold', horizontalalignment='center', transform=axs[0, 1].transAxes)
    axs[1, 0].text(0.9, 0.1, 'C', fontsize=10, weight='bold', horizontalalignment='center', transform=axs[1, 0].transAxes)
    axs[1, 1].text(0.95, 0.1, 'D', fontsize=10, weight='bold', horizontalalignment='center', transform=axs[1, 1].transAxes)
    plt.tight_layout()
