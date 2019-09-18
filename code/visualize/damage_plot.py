"""Make plots for supplement of manuscript
"""

import numpy as np
import matplotlib.pyplot as plt
import os

import als_process2019
ALS_SPEC = als_process2019.get_als_spectra()
import abs_get_processed_data
import abs_ana

# Conversion factor from FWHM fluence to average fluence:
FLUENCE_CONVERSION = abs_ana.FWHM_FLUENCE2AVERAGE_FLUENCE

import visualize.set_plot_params
visualize.set_plot_params.init_paper_small()

def make_plot():
    fluence_set = [(0, 25), (25, 35), (35, 45), (45, 75), (75, 300)]
    spectra = abs_get_processed_data.get_damage_specs()
    f, axs = plt.subplots(1, 2, sharex=True, sharey='row', figsize=(5.5, 3))
    for ind, spec in enumerate(spectra):
        fluences = fluence_set[ind]
        label = str(round(int(FLUENCE_CONVERSION*fluences[0])))+' to '+str(round(int(FLUENCE_CONVERSION*fluences[1])))+' mJ/cm$^2$'
        axs[0].plot(spec['minus']['phot'], spec['minus']['spec'], label=label)
        axs[1].plot(spec['plus']['phot'], spec['plus']['spec'])
    plt.gcf().legend(loc=7, title='Maximum Previous\nFluence',
                     frameon=False)
    _format_damage(axs)
    fig_path = os.path.join(os.path.dirname(__file__), 'figure_s3.eps')
    plt.savefig(fig_path, dpi=600)
    fig_path = os.path.join(os.path.dirname(__file__), 'figure_s3.jpg')
    plt.savefig(fig_path, dpi=600)

def _format_damage(axs):
    axs[0].set_xlabel(' ')
    axs[1].set_xlabel(' ')
    plt.gcf().text(0.41, 0.04, 'Photon Energy (eV)', ha='center')
    axs[0].set_ylabel('Intensity (a.u.)')
    axs[0].set_xlim([773.5, 782])
    axs[0].text(0.1, 0.9, 'A', fontsize=10, weight='bold',
                horizontalalignment='center', transform=axs[0].transAxes)
    axs[1].text(0.1, 0.9, 'B', fontsize=10, weight='bold',
                horizontalalignment='center', transform=axs[1].transAxes)
    plt.tight_layout()
    plt.gcf().subplots_adjust(right=0.73)
