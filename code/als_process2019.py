"""Load Loic's Co/Pd spectra recorded at ALS
"""

import numpy as np
import matplotlib.pyplot as plt
import os

Pd_ABS = 0.315    # Calculated from CXRO

ALS_FILE = os.path.join(os.path.dirname(__file__), '../data/synchrotron/ALS-CoPd-normed.txt')
# column 1: photon energy
# column 2: Cp
# column 3: Cm
# Column 4: 3x XAS
# Column 5: XMCD

# Photon energy values which were found by hand by plotting the XMCDs:
JO_L3 = 778
JO_L2 = 793.2
ALS_L3 = 777.5
ALS_L2 = 792.8

def get_raw_als_data():
    data = np.genfromtxt(ALS_FILE)
    return {'phot': data[:, 0],
            'Cp': data[:, 1],
            'Cm': data[:, 2],
            'xas': data[:, 3]/3,
            'xmcd': data[:, 4]}

def get_als_spectra():
    """Return spectra in expected format of plotting functions.
    """
    cal_als_data = get_cal_als_data()
    plus = {'phot': cal_als_data['phot'],
            'spec': cal_als_data['Cp']}
    minus = {'phot': cal_als_data['phot'],
             'spec': cal_als_data['Cm']}
    xas = {'phot': cal_als_data['phot'],
           'spec': cal_als_data['xas']}
    xmcd = {'phot': cal_als_data['phot'],
            'spec': cal_als_data['xmcd']}
    return {'plus': plus,
            'minus': minus,
            'xas': xas,
            'xmcd': xmcd}

def get_cal_als_data():
    raw_data = get_raw_als_data()
    phot_scale = (JO_L3-JO_L2)/(ALS_L3-ALS_L2)
    new_phot = (raw_data['phot']-ALS_L3)*phot_scale+JO_L3
    # scale_factor scales spectra to get them on same scale as LCLS
    scale_factor = 0.2158    # determined empirically
    # Source polarization accounts for finite polarization of ALS
    source_polarization = 0.8596   # determined empirically
    xas = raw_data['xas']*scale_factor
    xmcd = raw_data['xmcd']*scale_factor/source_polarization
    Cp = xas+xmcd/2
    Cm = xas-xmcd/2
    return {'phot': new_phot,
            'Cp': Cp+Pd_ABS,
            'Cm': Cm+Pd_ABS,
            'xas': xas+Pd_ABS,
            'xmcd': xmcd}

def make_als_plot():
    data = get_cal_als_data()
    plt.figure()
    plt.plot(data['phot'], data['xas'], label='XAS')
    plt.plot(data['phot'], data['xmcd'], label='XMCD')
    format_als_plot()

def format_als_plot():
    plt.legend(loc='upper right')
    plt.xlabel('Photon Energy (eV)')
    plt.ylabel('Intensity (a.u.)')
    plt.axhline(linestyle='--', color='k')
    plt.xlim((746, 840))
    plt.tight_layout()
    
