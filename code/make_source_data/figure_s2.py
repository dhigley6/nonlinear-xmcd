import numpy as np
import matplotlib.pyplot as plt
import xlsxwriter
import pickle

import abs_get_processed_data
import als_process2019
import abs_ana

RAW_MCP_AND_ANDOR_FILE = '../../data/raw/raw_mcp_and_andor.p'

def get_and_save_source_data(workbook=None):
    source_data = get_source_data()
    save_source_data(source_data, workbook)

def get_source_data():
    raw_mcp_and_andor = get_raw_mcp_and_andor_waveforms()
    no_sample_data = get_no_sample_data()
    als_spec = als_process2019.get_als_spectra()
    lcls_specs = abs_get_processed_data.get_incident_specs()
    source_data = {'Time (us)': raw_mcp_and_andor[b'mcp'][b'x'],
                   'Voltage (V)': raw_mcp_and_andor[b'mcp'][b'y'],
                   'Pixel': raw_mcp_and_andor[b'andor'][b'x'],
                   'Counts/1E3': raw_mcp_and_andor[b'andor'][b'y'],
                   'Raw I0': no_sample_data['Raw I0'],
                   'Raw I1': no_sample_data['Raw I1'],
                   'Corrected I0': no_sample_data['Corrected I0'],
                   'Corrected I1': no_sample_data['Corrected I1'],
                   'ALS Photon Energy': als_spec['xas']['phot'],
                   'ALS XAS': als_spec['xas']['spec'],
                   '2.7 mJ/cm^2 Photon Energy': lcls_specs[-3]['xas']['phot'],
                   '2.7 mJ/cm^2 XAS': lcls_specs[-3]['xas']['spec'],
                   '43 mJ/cm^2 Photon Energy': lcls_specs[-1]['xas']['phot'],
                   '43 mJ/cm^2 XAS': lcls_specs[-1]['xas']['spec']}
    return source_data

def get_raw_mcp_and_andor_waveforms():
    """Return raw mcp and andor waveforms that were
    saved from raw data
    """
    pickle_out = open(RAW_MCP_AND_ANDOR_FILE, 'rb')
    raw_mcp_and_andor = pickle.load(pickle_out, encoding='bytes')
    pickle_out.close()
    return raw_mcp_and_andor

def get_no_sample_data():
    data = abs_ana.get_data([108], mcp_filter=(0, 100))
    num_points = 100
    # Normalize data for plotting so that scaling is nice:
    mcp_norm = np.amax(data['mcp'][:num_points])
    andor_norm = np.amax(data['andor'][:num_points])
    mcp_raw = data['mcp_uncorrected'][:num_points]/mcp_norm
    mcp = data['mcp'][:num_points]/mcp_norm
    andor_raw = data['andor'][:num_points]/andor_norm
    andor = andor_raw    # Andor values are not corrected
    return {'Raw I0': mcp_raw,
            'Corrected I0': mcp,
            'Raw I1': andor_raw,
            'Corrected I1': andor}

def save_source_data(source_data, workbook):
    if workbook is None:
        workbook = xlsxwriter.Workbook('fig_s2.xlsx')
    fig_s2_sheet = workbook.add_worksheet('Supp. Fig. 2')
    for ind, (header, data) in enumerate(source_data.items()):
        fig_s2_sheet.write(0, ind, header)
        fig_s2_sheet.write_column(1, ind, data)
        
def sanity_plot():
    sd = get_source_data()
    f, axs = plt.subplots(2, 2)
    axs[0, 0].plot(sd['Time (us)'], sd['Voltage (V)'])
    axs[0, 1].plot(sd['Pixel'], sd['Counts/1E3'])
    axs[1, 0].scatter(sd['Raw I0'], sd['Raw I1'])
    axs[1, 0].scatter(sd['Corrected I0'], sd['Corrected I1'])
    axs[1, 1].plot(sd['ALS Photon Energy'], sd['ALS XAS'])
    axs[1, 1].plot(sd['2.7 mJ/cm^2 Photon Energy'], sd['2.7 mJ/cm^2 XAS'])
    axs[1, 1].plot(sd['43 mJ/cm^2 Photon Energy'], sd['43 mJ/cm^2 XAS'])