import numpy as np
import matplotlib.pyplot as plt
import xlsxwriter

from visualize import loop_plotter

def get_and_save_source_data(workbook=None):
    source_data = get_source_data()
    save_source_data(source_data, workbook)
        
def get_source_data():
    data = loop_plotter.get_data()
    data = loop_plotter.process_data(data)
    source_data = {'Applied Magnetic Field (T) (Increasing Applied Field)': data['increasing']['x'],
                    'XMCD (Increasing Applied Field)': data['increasing']['y'],
                    'Applied Magnetic Field (T) (Decreasing Applied Field)': data['decreasing']['x'],
                    'XMCD (Decreasing Applied Field)': data['decreasing']['y']}
    return source_data

def save_source_data(source_data, workbook=None):
    if workbook is None:
        workbook = xlsxwriter.Workbook('fig_s3.xlsx')
    fig_s1_sheet = workbook.add_worksheet('Supp. Fig. 3')
    for ind, (header, data) in enumerate(source_data.items()):
        fig_s1_sheet.write(0, ind, header)
        fig_s1_sheet.write_column(1, ind, data)
        
def sanity_plot():
    sd = get_source_data()
    f = plt.figure()
    plt.plot(sd['Applied Magnetic Field (T) (Increasing Applied Field)'],
             sd['XMCD (Increasing Applied Field)'])
    plt.plot(sd['Applied Magnetic Field (T) (Decreasing Applied Field)'],
             sd['XMCD (Decreasing Applied Field)'])