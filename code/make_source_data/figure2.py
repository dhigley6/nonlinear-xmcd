import numpy as np
import matplotlib.pyplot as plt
import xlsxwriter

import als_process2019
import abs_get_processed_data

# Spectra recorded at ALS Synchrotron Light Source:
ALS_SPEC = als_process2019.get_als_spectra()
# Spectra recorded with a constant incident X-ray fluence:
INCIDENT_SPECS = abs_get_processed_data.get_incident_specs()

def get_and_save_source_data(workbook=None):
    source_data = get_source_data()
    save_source_data(source_data, workbook)

def get_als_source_data():
    als_source_data = {'ALS photon energy': ALS_SPEC['xas']['phot'],
                       'ALS plus': ALS_SPEC['plus']['spec'],
                       'ALS minus': ALS_SPEC['minus']['spec'],
                       'ALS XAS': ALS_SPEC['xas']['spec'],
                       'ALS XMCD': ALS_SPEC['xmcd']['spec']}
    return als_source_data

def get_lcls_source_data():
    spec_0p01 = INCIDENT_SPECS[-4]
    spec_2p7 = INCIDENT_SPECS[-3]
    spec_19 = INCIDENT_SPECS[-2]
    spec_43 = INCIDENT_SPECS[-1]
    specs_dict = {'0.01 mJ/cm^2': spec_0p01,
                  '2.7 mJ/cm^2': spec_2p7,
                  '19 mJ/cm^2': spec_19,
                  '43 mJ/cm^2': spec_43}
    lcls_source_data = {}
    for fluence, spec in specs_dict.items():
        phot_str = ''.join([fluence, ' photon energy'])
        plus_str = ''.join([fluence, ' plus'])
        plus_diff_str = ''.join([fluence, ' plus difference from ALS'])
        minus_str = ''.join([fluence, ' minus'])
        minus_diff_str = ''.join([fluence, ' minus difference from ALS'])
        xas_str = ''.join([fluence, ' XAS'])
        xas_diff_str = ''.join([fluence, ' XAS difference from ALS'])
        xmcd_str = ''.join([fluence, ' XMCD'])
        xmcd_diff_str = ''.join([fluence, ' XMCD difference from ALS'])
        lcls_source_data.update({phot_str: spec['xas']['phot'],
                                 plus_str: spec['plus']['spec'],
                                 plus_diff_str: spec['plus']['diff'],
                                 minus_str: spec['minus']['spec'],
                                 minus_diff_str: spec['minus']['diff'],
                                 xas_str: spec['xas']['spec'],
                                 xas_diff_str: spec['xas']['diff'],
                                 xmcd_str: spec['xmcd']['spec'],
                                 xmcd_diff_str: spec['xmcd']['diff']})
    return lcls_source_data

def get_source_data():
    als_source_data = get_als_source_data()
    lcls_source_data = get_lcls_source_data()
    source_data = {**als_source_data, **lcls_source_data}
    return source_data


def save_source_data(source_data, workbook=None):
    if workbook is None:
        workbook = xlsxwriter.Workbook('fig2.xlsx')
    fig2_sheet = workbook.add_worksheet('Fig. 2')
    for ind, (header, data) in enumerate(source_data.items()):
        fig2_sheet.write(0, ind, header)
        fig2_sheet.write_column(1, ind, data)
    

def sanity_plot():
    """Plot source data to make sure it makes sense
    """
    sd = get_source_data()
    f, axs = plt.subplots(2, 4)
    axs[0, 0].plot(sd['ALS photon energy'], sd['ALS plus'])
    axs[0, 1].plot(sd['ALS photon energy'], sd['ALS minus'])
    axs[0, 2].plot(sd['ALS photon energy'], sd['ALS XAS'])
    axs[0, 3].plot(sd['ALS photon energy'], sd['ALS XMCD'])
    fluence_strs = ['0.01 mJ/cm^2', '2.7 mJ/cm^2', '19 mJ/cm^2', '43 mJ/cm^2']
    for fluence_str in fluence_strs:
        phot_str = ''.join([fluence_str, ' photon energy'])
        plus_str = ''.join([fluence_str, ' plus'])
        plus_diff_str = ''.join([fluence_str, ' plus difference from ALS'])
        minus_str = ''.join([fluence_str, ' minus'])
        minus_diff_str = ''.join([fluence_str, ' minus difference from ALS'])
        xas_str = ''.join([fluence_str, ' XAS'])
        xas_diff_str = ''.join([fluence_str, ' XAS difference from ALS'])
        xmcd_str = ''.join([fluence_str, ' XMCD'])
        xmcd_diff_str = ''.join([fluence_str, ' XMCD difference from ALS'])
        axs[0, 0].plot(sd[phot_str], sd[plus_str])
        axs[1, 0].plot(sd[phot_str], sd[plus_diff_str])
        axs[0, 1].plot(sd[phot_str], sd[minus_str])
        axs[1, 1].plot(sd[phot_str], sd[minus_diff_str])
        axs[0, 2].plot(sd[phot_str], sd[xas_str])
        axs[1, 2].plot(sd[phot_str], sd[xas_diff_str])
        axs[0, 3].plot(sd[phot_str], sd[xmcd_str])
        axs[1, 3].plot(sd[phot_str], sd[xmcd_diff_str])