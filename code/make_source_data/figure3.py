import numpy as np
import matplotlib.pyplot as plt
import xlsxwriter

import als_process2019
import abs_get_processed_data

# Spectra recorded at ALS Synchrotron Light Source:
ALS_SPEC = als_process2019.get_als_spectra()
# Spectra recorded with a constant absorbed X-ray fluence:
ABSORBED_SPECS = abs_get_processed_data.get_constant_absorbed_specs()

def get_and_save_source_data(workbook=None):
    source_data = get_source_data()
    save_source_data(source_data, workbook)

def get_source_data():
    als_source_data = get_als_source_data()
    lcls_source_data = get_lcls_source_data()
    source_data = {**als_source_data, **lcls_source_data}
    return source_data

def get_als_source_data():
    als_source_data = {'ALS photon energy': ALS_SPEC['xas']['phot'],
                       'ALS XAS': ALS_SPEC['xas']['spec'],
                       'ALS XMCD': ALS_SPEC['xmcd']['spec']}
    return als_source_data

def get_lcls_source_data():
    spec = ABSORBED_SPECS[-1]
    abs_energy_density = '140 meV/atom'
    phot_str = ''.join([abs_energy_density, ' photon energy'])
    xas_str = ''.join([abs_energy_density, ' XAS'])
    xas_diff_str = ''.join([abs_energy_density, ' XAS difference from ALS'])
    xas_diff_err_str = ''.join([abs_energy_density, ' XAS difference error'])
    xmcd_str = ''.join([abs_energy_density, ' XMCD'])
    xmcd_diff_str = ''.join([abs_energy_density, ' XMCD difference from ALS'])
    xmcd_diff_err_str = ''.join([abs_energy_density, ' XMCD difference error'])
    lcls_source_data = {phot_str: spec['xas']['phot'],
                        xas_str: spec['xas']['spec'],
                        xas_diff_str: spec['xas']['diff'],
                        xas_diff_err_str: np.sqrt((spec['plus']['err'])**2+(spec['minus']['err'])**2)/2,
                        xmcd_str: spec['xmcd']['spec'],
                        xmcd_diff_str: spec['xmcd']['diff'],
                        xmcd_diff_err_str: np.sqrt((spec['plus']['err'])**2+(spec['minus']['err'])**2)}
    return lcls_source_data
    

def save_source_data(source_data, workbook=None):
    if workbook is None:
        workbook = xlsxwriter.Workbook('fig3.xlsx')
    fig3_sheet = workbook.add_worksheet('Fig. 3')
    for ind, (header, data) in enumerate(source_data.items()):
        fig3_sheet.write(0, ind, header)
        fig3_sheet.write_column(1, ind, data)
        
def sanity_plot():
    sd = get_source_data()
    f, axs = plt.subplots(2, 1, sharex=True)
    axs[0].plot(sd['ALS photon energy'], sd['ALS XAS'])
    axs[0].plot(sd['ALS photon energy'], sd['ALS XMCD'])
    axs[0].plot(sd['140 meV/atom photon energy'], sd['140 meV/atom XAS'])
    axs[0].plot(sd['140 meV/atom photon energy'], sd['140 meV/atom XMCD'])
    axs[1].errorbar(sd['140 meV/atom photon energy'], sd['140 meV/atom XAS difference from ALS'], 
                    yerr=sd['140 meV/atom XAS difference error'])
    axs[1].errorbar(sd['140 meV/atom photon energy'], sd['140 meV/atom XMCD difference from ALS'],
                    yerr=sd['140 meV/atom XMCD difference error'])
                    