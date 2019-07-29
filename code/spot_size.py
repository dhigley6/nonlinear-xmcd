"""Code for calculating spot size
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import abs_ana
import norm_data

def fit_gauss(x, y, weights=None):
    """Fit gaussian to data
    """
    def gauss(x, a, x0, sigma):
        return a*np.exp(-(x-x0)**2/(2*sigma**2))
    amp_guess = np.amax(y)
    mean_guess = np.sum(x*y)/np.sum(y)
    popt, pcov = scipy.optimize.curve_fit(gauss, x, y, p0=[amp_guess, mean_guess, 0.3])
    return popt, pcov, gauss

def plot_profile(bin_data, sigs, norms):
    """Plot beam profile with Gaussian fit
    """
    profile = norm_data.norm_disc_scan(bin_data, sigs, norms)
    popt, pcov, fit_func = fit_gauss(profile['bin_centers'], profile['norm_sig'])
    plt.figure()
    plt.plot(profile['bin_centers'], profile['norm_sig'], label='Scan')
    fwhm = 2.35482*popt[2]
    plt.plot(profile['bin_centers'],
             fit_func(profile['bin_centers'], *popt),
             label='Gaussian Fit, FWHM='+str(fwhm))
    plt.legend(loc='best')

def plot_profile_lk30(run, coord='sam_x'):
    data = abs_ana.get_data(([run], None, (0, 100)))
    plot_profile(data[coord], data['andor'], data['mcp'])
    plt.title('LK30 Spot size scan, run '+str(run))
    plt.xlabel(coord+' (mm)')
    plt.ylabel('Normalized Value')
    save_file_start = '../plots/2016_04_07_spot_size_scan_'
    save_file = ''.join([save_file_start, coord[-1], '_run', str(run), '.png'])
    plt.savefig(save_file)
