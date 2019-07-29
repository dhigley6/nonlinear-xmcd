"""Make overall manuscript plot for submission to Nature Communications in 2019
D. J. Higley (2019)
"""

import numpy as np
import matplotlib.pyplot as plt
import os

import visualize.dummy_ax
import visualize.set_plot_params
visualize.set_plot_params.init_paper_small()

import als_process2019
import abs_get_processed_data

# Spectra recorded at ALS Synchrotron Light Source:
ALS_SPEC = als_process2019.get_als_spectra()
# Spectra recorded with a constant incident X-ray fluence:
INCIDENT_SPECS = abs_get_processed_data.get_incident_specs()
# Spectra recorded with a constant absorbed X-ray fluence:
ABSORBED_SPECS = abs_get_processed_data.get_constant_absorbed_specs()
# Short dash linestyle:
SD = (0, (3, 1.5))

Pd_ABS = 0.315    # Calculated from CXRO

FERMI_ENERGY = 777.5

def make_plot():
    """Plot with constant incident fluence spectra
    """

    def plus_minus_plot(high_spec, ax):
        ax.plot(ALS_SPEC['minus']['phot'], ALS_SPEC['minus']['spec'], linestyle=SD, color='m', label='$\sigma_-$\nSync.')
        ax.plot(ALS_SPEC['plus']['phot'], ALS_SPEC['plus']['spec'], linestyle=SD, color='c', label='$\sigma_+$\nSync.')
        ax.plot(high_spec['minus']['phot'], high_spec['minus']['spec'], color='m', label='$\sigma_-$\nXFEL')
        ax.plot(high_spec['plus']['phot'], high_spec['plus']['spec'], color='c', label='$\sigma_+$\nXFEL')
        ax.plot([773, 783], [Pd_ABS, Pd_ABS], color='g')

    def xas_xmcd_plot(high_spec, ax):
        ax.plot(ALS_SPEC['xas']['phot'], ALS_SPEC['xas']['spec'], linestyle=SD, color='b', label='XAS\nSync.')
        ax.plot(high_spec['xas']['phot'], high_spec['xas']['spec'], label='XAS\nXFEL', color='b')
        ax.plot(ALS_SPEC['xmcd']['phot'], -1*ALS_SPEC['xmcd']['spec'], linestyle=SD, color='r', label='XMCD\nSync.')
        ax.plot(high_spec['xmcd']['phot'], -1*high_spec['xmcd']['spec'], color='r', label='XMCD\nXFEL')
        
    def plus_minus_difference_plot(specs, ax):
        offset = 0.3
        for ind, spec in enumerate(specs):
            if ind == 0:
                ax.plot(spec['minus']['phot'], spec['minus']['diff']+offset*ind, color='m', label='$\Delta\sigma_-$')
                ax.plot(spec['plus']['phot'], spec['plus']['diff']+offset*ind, color='c', label='$\Delta\sigma_+$')
            else:
                ax.plot(spec['minus']['phot'], spec['minus']['diff']+offset*ind, color='m')
                ax.plot(spec['plus']['phot'], spec['plus']['diff']+offset*ind, color='c')
            ax.axhline(offset*ind, color='k', linestyle='--')

    def xas_xmcd_difference_plot(specs, ax):
        offset = 0.3
        for ind, spec in enumerate(specs):
            if ind == 0:
                ax.plot(spec['xas']['phot'], spec['xas']['diff']*2+offset*ind, color='b', label='$2x\Delta$XAS')
                ax.plot(spec['xmcd']['phot'], (spec['xmcd']['diff'])*-1+offset*ind, color='r', label='-$\Delta$XMCD')
            else: 
                ax.plot(spec['xas']['phot'], spec['xas']['diff']*2+offset*ind, color='b')
                ax.plot(spec['xmcd']['phot'], (spec['xmcd']['diff'])*-1+offset*ind, color='r')
            ax.axhline(offset*ind, color='k', linestyle='--')

    high_spec = INCIDENT_SPECS[-1]
    f, axs = plt.subplots(2, 2, sharex=True, sharey='row', figsize=(3.37, 4), gridspec_kw={'height_ratios': [1, 1.1]})
    plus_minus_plot(high_spec, axs[0, 0])
    xas_xmcd_plot(high_spec, axs[0, 1])
    specs = [high_spec, INCIDENT_SPECS[-2], INCIDENT_SPECS[-3], INCIDENT_SPECS[-4]]
    #specs = [high_spec, INCIDENT_SPECS[-2], INCIDENT_SPECS[-3]]
    plus_minus_difference_plot(specs, axs[1, 0])
    xas_xmcd_difference_plot(specs, axs[1, 1])
    format_overall_plot(f, axs)
    fig_path = os.path.join(os.path.dirname(__file__), 'figure2.eps')
    plt.savefig(fig_path, dpi=600)

def format_overall_plot(f, axs):
    axs[0, 0].set_xlim((773.5, 782))
    axs[0, 0].set_ylim((-0.05, 1.6))
    axs[1, 0].set_ylim((-0.4, 1.1))
    axs[0, 0].set_ylabel('Intensity')
    axs[1, 0].set_ylabel('Intensity Change')
    axs[1, 0].set_xticks((775, 776, 777, 778, 779, 780, 781))
    axs[1, 0].set_xticklabels(['775', '', '777', '', '779', '', '781'])
    axs[1, 0].legend(loc='lower left', borderpad=0.3, framealpha=1, borderaxespad=0, frameon=False, fontsize=7)
    axs[1, 1].legend(loc='lower left', borderpad=0.3, framealpha=1, borderaxespad=0, frameon=False, fontsize=7)
    axs[0, 0].axvline(FERMI_ENERGY, linewidth=1, color='k', linestyle='--')
    axs[0, 1].axvline(FERMI_ENERGY, linewidth=1, color='k', linestyle='--')
    axs[1, 0].axvline(FERMI_ENERGY, linewidth=1, color='k', linestyle='--')
    axs[1, 1].axvline(FERMI_ENERGY, ymin=0.1, linewidth=1, color='k', linestyle='--')
    axs[0, 0].axhline(0, color='k', linewidth=1, linestyle='--')
    axs[0, 1].axhline(0, color='k', linewidth=1, linestyle='--')
    axs[0, 0].text(0.15, 0.9, 'A', fontsize=10, weight='bold', horizontalalignment='center', transform=axs[0, 0].transAxes)
    axs[1, 0].text(0.15, 0.91, 'C', fontsize=10, weight='bold', horizontalalignment='center', transform=axs[1, 0].transAxes)
    axs[0, 1].text(0.85, 0.9, 'B', fontsize=10, weight='bold', horizontalalignment='center', transform=axs[0, 1].transAxes)
    axs[1, 1].text(0.85, 0.91, 'D', fontsize=10, weight='bold', horizontalalignment='center', transform=axs[1, 1].transAxes)
    axs[0, 0].text(780, 0.36, 'Pd', horizontalalignment='center', transform=axs[0, 0].transData, color='g')
    axs[0, 0].text(777, 1.3, '$\sigma_-$', horizontalalignment='right', transform=axs[0, 0].transData, color='m')
    axs[0, 0].text(779, 0.6, '$\sigma_+$', horizontalalignment='right', transform=axs[0, 0].transData, color='c')
    axs[0, 1].text(775.2, 0.55, 'XAS\n' r'$\dfrac{\sigma_++\sigma_-}{2}$', horizontalalignment='center', transform=axs[0, 1].transData, color='b')
    axs[0, 1].text(775.2, 0.07, '-XMCD\n$\sigma_--\sigma_+$', horizontalalignment='center', transform=axs[0, 1].transData, color='r')
    plt.text(782.25, 0.975, '0.01 mJ/cm$^2$', horizontalalignment='center', transform=axs[1, 0].transData, color='k', clip_on=False, fontsize=7, bbox=dict(facecolor='white', edgecolor='none', pad=0))
    plt.text(782.25, 0.675, '2.7 mJ/cm$^2$', horizontalalignment='center', transform=axs[1, 0].transData, color='k', clip_on=False, fontsize=7, bbox=dict(facecolor='white', edgecolor='none', pad=0))
    #axs[1, 0].text(782, 0.33, '3.7', horizontalalignment='right', transform=axs[1, 0].transData, color='k', clip_on=False, fontsize=7)
    plt.text(782.25, 0.375, '19 mJ/cm$^2$', horizontalalignment='center', transform=axs[1, 0].transData, color='k', clip_on=False, fontsize=7, bbox=dict(facecolor='white', edgecolor='none', pad=0))
    plt.text(782.25, 0.075, '43 mJ/cm$^2$', horizontalalignment='center', transform=axs[1, 0].transData, color='k', clip_on=False, fontsize=7, bbox=dict(facecolor='white', edgecolor='none', pad=0))
    #axs[1, 1].text(773.6, 0.96, '0.01 mJ/cm$^2$', horizontalalignment='left', transform=axs[1, 1].transData, color='k', clip_on=False, fontsize=7)
    #axs[1, 1].text(773.6, 0.66, '2.5 mJ/cm$^2$', horizontalalignment='left', transform=axs[1, 1].transData, color='k', clip_on=False, fontsize=7)
    #axs[1, 1].text(773.5, 0.33, '3.7', horizontalalignment='left', transform=axs[1, 1].transData, color='k', clip_on=False, fontsize=7)
    #axs[1, 1].text(773.6, 0.34, '19\nmJ/cm$^2$', horizontalalignment='left', transform=axs[1, 1].transData, color='k', clip_on=False, fontsize=7, linespacing=0.75)
    #axs[1, 1].text(773.6, 0.03, '43\nmJ/cm$^2$', horizontalalignment='left', transform=axs[1, 1].transData, color='k', clip_on=False, fontsize=7, linespacing=0.75)

    axs[1, 0].set_xlabel('dummy')
    ax_fermi0 = axs[0, 0].twiny()
    ax_fermi1 = axs[0, 1].twiny()
    fermi_axs = [ax_fermi0, ax_fermi1]
    for fermi_ax in fermi_axs:
        fermi_ax.set_xlim((-4, 4.5))
        fermi_ax.set_xticks((-4, -2, 0, 2, 4))
        #fermi_ax.set_xlabel(r'E-E$_{core\rightarrow Fermi}$ (eV)')
        fermi_ax.set_xlabel(' ')
    plt.tight_layout(w_pad=-0.3, h_pad=0.2, pad=0.2)
    axs[1, 0].set_xlabel('')
    big_ax = visualize.dummy_ax.make_big_dummy_ax(f)
    big_ax.plot([], [], color='k', linestyle=SD, label='ALS')
    big_ax.plot([], [], color='k', linestyle='-', label='43 mJ/cm$^2$')
    big_ax.set_xlabel('Photon Energy (eV)')
    big_ax.legend(loc='upper center', frameon=True, borderpad=0.3, framealpha=1, borderaxespad=0.3)
    big_ax.set_yticks([])
    big_ax.set_xticks([0])
    big_ax.set_xticklabels([' '])
    big_ax.xaxis.set_label_coords(0.5, -0.065)
    big_fermi_ax = big_ax.twiny()
    big_fermi_ax.set_xlabel(r'Photon Energy Above E$_{core\rightarrow Fermi}$ (eV)')
    big_fermi_ax.set_xticks([0])
    big_fermi_ax.set_xticklabels([' '])
