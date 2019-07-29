"""Make quantitative plot for 2019 submission to Nat. Comm.
"""

import numpy as np
import matplotlib.pyplot as plt
import os

from visualize import dummy_ax
import abs_get_processed_data
import visualize.set_plot_params
visualize.set_plot_params.init_paper_small()

def make_figure():
    f, axs = plt.subplots(1, 2, sharex=True, figsize=(3.37, 2.4))
    quant = abs_get_processed_data.get_spec_quantification()
    dummy = np.array([0, 200])
    axs[0].errorbar(quant['pulse_averaged_absorbed_energies'], quant['spec_energies'], marker='o', label='Expt.', yerr=quant['err_spec_energies'], xerr=0.2*quant['pulse_averaged_absorbed_energies'], linestyle='None', markerfacecolor='none', markeredgecolor='k', color='k', markeredgewidth=1, markersize=3, capsize=2)
    axs[0].plot(dummy, dummy*quant['fraction_within_2eV'], label=r'Fit, $\tau_c = 13$ fs', linestyle=(0, (3, 1.5)), color='r')
    axs[1].errorbar(quant['pulse_averaged_absorbed_energies'], quant['xmcds']*100, marker='o', label='Expt.', yerr=quant['err_xmcds']*100, xerr=0.2*quant['pulse_averaged_absorbed_energies'], linestyle='None', markerfacecolor='k', markeredgecolor='k', color='k', markeredgewidth=1, markersize=3, capsize=2)
    format_figure(f, axs)
    fig_path = os.path.join(os.path.dirname(__file__), 'figure4.eps')
    plt.savefig(fig_path, dpi=600)
    
def format_figure(f, axs):
    axs[0].legend(loc='upper left', fontsize=7)
    axs[0].legend(loc='upper left', fontsize=7)
    axs[0].text(0.85, 0.15, 'A', fontsize=10, weight='bold', horizontalalignment='right', transform=axs[0].transAxes)
    axs[1].text(0.35, 0.15, 'B', fontsize=10, weight='bold', horizontalalignment='right', transform=axs[1].transAxes)
    axs[1].legend(loc='upper right', fontsize=7)
    axs[1].set_ylabel('Normalized XMCD (%)')
    axs[0].set_ylabel('Co 3d Electron Energy\n(meV/atom)')
    big_ax = visualize.dummy_ax.make_big_dummy_ax(f)
    big_ax.set_xlabel('Pulse Averaged Absorbed\nX-Ray Energy (meV/atom)')
    big_ax.xaxis.set_label_coords(0.5, -0.2)
    big_ax.set_yticks([])
    big_ax.set_xticks([0])
    big_ax.set_xticklabels([' '])
    axs[0].set_xticks([0, 50, 100, 150])
    axs[1].set_xticks([0, 50, 100, 150])
    axs[0].set_xlim([-5, 175])
    axs[0].set_ylim([-25, 175])
    plt.tight_layout(pad=0.2, w_pad=2.3, h_pad=0)
