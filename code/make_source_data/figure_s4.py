"""Save source data for Fig. S4
Note that this relies on the results of notebook 'demag.ipynb'
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
import xlsxwriter
import os

DEMAG_RESULTS_FILE = os.path.join(os.path.dirname(__file__), '../../data/demag_simulation.p')

def get_and_save_source_data(workbook=None):
    source_data = get_source_data()
    save_source_data(source_data, workbook)

def get_source_data():
    pickle_out = open(DEMAG_RESULTS_FILE, 'rb')
    source_data = pickle.load(pickle_out)
    pickle_out.close()
    return source_data

def save_source_data(source_data, workbook=None):
    if workbook is None:
        workbook = xlsxwriter.Workbook('fig_s4.xlsx')
    fig_s4_sheet = workbook.add_worksheet('Supp. Fig. 4')
    for ind, (header, data) in enumerate(source_data.items()):
        fig_s4_sheet.write(0, ind, header)
        fig_s4_sheet.write_column(1, ind, data)

def sanity_plot():
    sd = get_source_data()
    f, axs = plt.subplots(2, 1, sharex=True)
    axs[0].plot(sd['time (fs)'], sd['pulse'])
    axs[0].plot(sd['time (fs)'], sd['Co M'])
    axs[0].plot(sd['time (fs)'], sd['Co/Pt M'])
    axs[1].plot(sd['time (fs)'], sd['Co Te (K)'])
    axs[1].plot(sd['time (fs)'], sd['Co/Pt Te (K)'])