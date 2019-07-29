"""Make constant absorption plot for submission to Nature Communications
D. J. Higley (2019)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import os

import visualize.dummy_ax
import visualize.set_plot_params
visualize.set_plot_params.init_paper_small()

import als_process2019
import abs_get_processed_data

# Spectra recorded at ALS Synchrotron Light Source:
ALS_SPEC = als_process2019.get_als_spectra()
# Spectra recorded with a constant absorbed X-ray fluence:
ABSORBED_SPECS = abs_get_processed_data.get_constant_absorbed_specs()
# Short dash linestyle:
SD = (0, (3, 1.5))

FERMI_ENERGY = 777.5

def constant_plot2():
    """Plot with constant absorbed fluence spectra
    """
    f, axs = plt.subplots(2, 1, sharex=True, sharey=False, figsize=(3.37, 3.5), gridspec_kw={'height_ratios': [1, 1]})
    spec = ABSORBED_SPECS[-1]
    axs[0].plot(ALS_SPEC['xas']['phot'], (ALS_SPEC['xas']['spec']), linestyle=SD, color='b')
    axs[0].plot(spec['xas']['phot'], (spec['xas']['spec']), linestyle='-', color='b')
    Te, height = fit_fermi_dirac_change(spec['xas'])
    print(''.join(['Electronic temperature fit of ', str(Te), ' K']))
    fermi_dirac = height*fermi_dirac_change(np.arange(770, 790, 0.05), FERMI_ENERGY, Te)
    axs[1].plot(np.arange(770, 790, 0.05), 2*fermi_dirac, linestyle='-', color='b')
    yerr = np.sqrt((spec['plus']['err'])**2+(spec['minus']['err'])**2)
    axs[1].errorbar(spec['xas']['phot'], 2*spec['xas']['diff'], yerr=yerr, marker='o', linestyle='None', markeredgecolor='b', color='b', markeredgewidth=1, markersize=2, capsize=2)
    axs[0].plot(ALS_SPEC['xmcd']['phot'], ALS_SPEC['xmcd']['spec']*-1, color='r', linestyle=SD)
    axs[0].plot(spec['xmcd']['phot'], spec['xmcd']['spec']*-1, color='r', linestyle='-')
    offset = 0.25
    axs[1].plot(ALS_SPEC['xmcd']['phot'], ALS_SPEC['xmcd']['spec']*0.22-offset, color='r')
    yerr = np.sqrt((spec['plus']['err'])**2+(spec['minus']['err'])**2)
    axs[1].errorbar(spec['xmcd']['phot'], spec['xmcd']['diff']*-1-offset, yerr=yerr, marker='o', linestyle='None', markeredgecolor='r', color='r', markeredgewidth=1, markersize=2, capsize=2)
    format_constant_plot(f, axs)
    fig_path = os.path.join(os.path.dirname(__file__), 'figure3.eps')
    plt.savefig(fig_path, dpi=600)

def format_constant_plot(f, axs):
    axs[0].set_xlim((773.5, 782))
    axs[1].set_xlabel('Photon Energy (eV)')
    axs[0].set_ylabel('Intensity')
    axs[1].set_ylabel('Intensity\nChange')
    axs[0].axvline(FERMI_ENERGY, linewidth=1, color='k', linestyle='--')
    axs[1].axvline(FERMI_ENERGY, linewidth=1, color='k', linestyle='--')
    sync_line = plt.Line2D([], [], linewidth=1, color='k', linestyle=SD)
    xfel_line = plt.Line2D([], [], linewidth=1, color='k', linestyle='-')
    axs[0].legend([sync_line, xfel_line], ['ALS', 'XFEL'], loc='upper right')
    #data_line = plt.Line2D([], [], marker='o', color='k', linestyle='', markersize=2, capsize=2)
    data_line = plt.errorbar([], [], yerr=[], marker='o', linestyle='None', markeredgecolor='k', color='k', markeredgewidth=1, markersize=2, capsize=2)
    fit_line = plt.Line2D([], [], color='k', linestyle='-')
    axs[1].legend([data_line, fit_line], ['Expt.', 'Fit'], loc='upper right', ncol=2, columnspacing=1)
    axs[0].text(0.05, 0.85, 'A', fontsize=10, weight='bold', horizontalalignment='center', transform=axs[0].transAxes)
    axs[1].text(0.05, 0.85, 'B', fontsize=10, weight='bold', horizontalalignment='center', transform=axs[1].transAxes)
    axs[0].text(777, 0.9, 'XAS', horizontalalignment='right', transform=axs[0].transData, color='b')
    axs[0].text(776.8, 0.1, '-XMCD', horizontalalignment='right', transform=axs[0].transData, color='r')
    axs[1].text(774.6, 0.12, '$2x\Delta$XAS', horizontalalignment='left', transform=axs[1].transData, color='b')
    axs[1].text(774, -0.2, '-$\Delta$XMCD', horizontalalignment='left', transform=axs[1].transData, color='r')
    plt.tight_layout(h_pad=0)

def fermi_dirac_change(phot, mu, T, T0=300.0, broad=0.43):
    # Make sample vector which is at least 10*broad longer than phot on each side
    if np.size(phot) < 2:
        diff = 0.01
    else:
        diff = np.amin(np.diff(phot))
    sample_phot = np.arange(-21*broad+np.amin(phot), 21*broad+np.amax(phot)+diff/10, diff/10)
    print(len(sample_phot))
    # Co L3 lifetime broadening is 0.43 eV
    k = 8.617E-5
    fermi_dirac_start = 1/(np.exp((sample_phot-mu)/(k*T0))+1)
    fermi_dirac_end = 1/(np.exp((sample_phot-mu)/(k*T))+1)
    fermi_dirac_change = fermi_dirac_start-fermi_dirac_end
    lorentz_points = np.arange(-10*broad, 10*broad+diff, diff)
    lorentz_points = sample_phot-np.mean(sample_phot)
    lorentz_broad = lorentzian(lorentz_points, broad)
    lorentz_broad = lorentz_broad/np.sum(lorentz_broad)
    fermi_dirac_change = np.convolve(fermi_dirac_change, lorentz_broad, mode='same')
    # 260 meV broadened Gaussian to account for experimental resolving power of 3000
    gauss_broad_points = np.arange(-20*broad, 20*broad+sample_phot[1]-sample_phot[0], sample_phot[1]-sample_phot[0])
    gauss_broad = np.exp(-1*(gauss_broad_points)**2/(2*0.11)**2)
    gauss_broad = gauss_broad/np.sum(gauss_broad)
    fermi_dirac_change = np.convolve(fermi_dirac_change, gauss_broad, mode='same')
    fermi_dirac_change = np.interp(phot, sample_phot, fermi_dirac_change)
    return fermi_dirac_change

def fit_fermi_dirac_change(spec):

    def error_func(pars):
        Te = pars[0]
        height = pars[1]
        height = 0.9
        fermi_change = height*fermi_dirac_change(spec['phot'], FERMI_ENERGY, Te)
        error = spec['diff']-fermi_change
        return error

    pars, cov = leastsq(error_func, [3000, 1])
    print(pars)
    return pars

def lorentzian(lorentz_points, broad):
    return 1.0/(1.0+(lorentz_points/(broad/2))**2)
