"""Getting and modifying run parameters
"""

import numpy as np
import os
CAL_PATH = os.path.join(os.path.dirname(__file__), '../cals/abs_run_params3.csv')

def get_run_pars(run):
    pars = get_pars()
    row = get_run_row_num(run, pars)
    return pars[row]

def change_run_pars(run, par_dict):
    pars = get_pars()
    row = get_run_row_num(run, pars)
    for key in par_dict.keys():
        pars[row][key] = par_dict[key]
    save_pars(pars)

def get_pars():
    pars = np.genfromtxt(CAL_PATH, names=True, delimiter=',')
    return pars

def get_run_row_num(run, pars):
    run_start_col = pars['run_start']
    run_end_col = pars['run_end']
    rows_greater_or_equal = run >= run_start_col
    rows_less_or_equal = run <= run_end_col
    row_nums = np.arange(len(pars))
    row_num_want = row_nums[rows_greater_or_equal & rows_less_or_equal]
    row_num_want = int(row_num_want)
    return row_num_want
    
def save_pars(pars):
    header = make_header(pars)
    np.savetxt(CAL_PATH, pars, delimiter=',', header=header)

def make_header(pars):
    names = pars.dtype.names
    header = ''
    for name in names:
        header = header+name+','
    # Take out final comma:
    header = header[:-1]
    return header
