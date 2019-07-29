import numpy as np
import matplotlib.pyplot as plt
import xlsxwriter

import abs_get_processed_data

def get_and_save_source_data(workbook=None):
    source_data = get_source_data()
    save_source_data(source_data, workbook)

def get_source_data():
    specs = abs_get_processed_data.get_damage_specs()
    source_data = {'0 to 11 mJ/cm^2 Photon Energy': specs[0]['minus']['phot'],
                   '0 to 11 mJ/cm^2 Plus': specs[0]['plus']['spec'],
                   '0 to 11 mJ/cm^2 Minus': specs[0]['minus']['spec'],
                   '11 to 15 mJ/cm^2 Photon Energy': specs[1]['minus']['phot'],
                   '11 to 15 mJ/cm^2 Plus': specs[1]['plus']['spec'],
                   '11 to 15 mJ/cm^2 Minus': specs[1]['minus']['spec'],
                   '15 to 19 mJ/cm^2 Photon Energy': specs[2]['minus']['phot'],
                   '15 to 19 mJ/cm^2 Plus': specs[2]['plus']['spec'],
                   '15 to 19 mJ/cm^2 Minus': specs[2]['minus']['spec'],
                   '19 to 33 mJ/cm^2 Photon Energy': specs[3]['minus']['phot'],
                   '19 to 33 mJ/cm^2 Plus': specs[3]['plus']['spec'],
                   '19 to 33 mJ/cm^2 Minus': specs[3]['minus']['spec'],
                   '33 to 132 mJ/cm^2 Photon Energy': specs[4]['minus']['phot'],
                   '33 to 132 mJ/cm^2 Plus': specs[4]['plus']['spec'],
                   '33 to 132 mJ/cm^2 Minus': specs[4]['minus']['spec']}
    return source_data
    
def save_source_data(source_data, workbook=None):
    if workbook is None:
        workbook = xlsxwriter.Workbook('fig_s3.xlsx')
    fig_s3_sheet = workbook.add_worksheet('Supp. Fig. 3')
    for ind, (header, data) in enumerate(source_data.items()):
        fig_s3_sheet.write(0, ind, header)
        fig_s3_sheet.write_column(1, ind, data)
    
def sanity_plot():
    specs = abs_get_processed_data.get_damage_specs()
    f, axs = plt.subplots(1, 2, sharex=True)
    for spec in specs:
        axs[0].plot(spec['minus']['phot'], spec['minus']['spec'])
        axs[1].plot(spec['minus']['phot'], spec['plus']['spec'])