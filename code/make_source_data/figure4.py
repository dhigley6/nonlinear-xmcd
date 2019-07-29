import numpy as np
import matplotlib.pyplot as plt
import xlsxwriter

import abs_get_processed_data

def get_and_save_source_data(workbook=None):
    source_data = get_source_data()
    save_source_data(source_data, workbook)

def get_source_data():
    quant = abs_get_processed_data.get_spec_quantification()
    source_data = { 'Pulse averaged absorbed energy (meV)': quant['pulse_averaged_absorbed_energies'],
                    'Co 3d Energy Per Atom (meV)': quant['spec_energies'],
                    'Standard error of Co 3d energy per atom (meV)': quant['err_spec_energies'],
                    'XMCD': quant['xmcds'],
                    'Standard error of XMCD': quant['err_xmcds'],
                    'Pulse averaged absorbed energy error (meV)': quant['pulse_averaged_absorbed_energies']*0.2}
    return source_data
    

def save_source_data(source_data, workbook=None):
    if workbook is None:
        workbook = xlsxwriter.Workbook('fig4.xlsx')
    fig4_sheet = workbook.add_worksheet('Fig. 4')
    for ind, (header, data) in enumerate(source_data.items()):
        fig4_sheet.write(0, ind, header)
        fig4_sheet.write_column(1, ind, data)
        
def sanity_plot():
    sd = get_source_data()
    f, axs = plt.subplots(1, 2, sharex=True)
    axs[0].errorbar(sd['Pulse averaged absorbed energy (meV)'], sd['Co 3d Energy Per Atom (meV)'],
                    xerr=sd['Pulse averaged absorbed energy error (meV)'],
                    yerr=sd['Standard error of Co 3d energy per atom (meV)'])
    axs[1].errorbar(sd['Pulse averaged absorbed energy (meV)'], sd['XMCD'],
                    xerr=sd['Pulse averaged absorbed energy error (meV)'],
                    yerr=sd['Standard error of XMCD'])
    