"""Submit multiple LCLS nonlinear absorption runs to the batch farm for pre-processing

Usage:
    abs_submit_multiple_runs.py <start_run> <end_run> [<maxevt>]
"""


import argparse
import os

# Configuration parameters for how to run the jobs to extract data
# Number of cores to run each job on
NUMBER_OF_CORES = 6
# Number of processers on each host for each job
PROCESSERS_PER_HOST = 2
# Queue to use
QUEUE_NAME = 'psanaq'

# Read in command line arguments
parser = argparse.ArgumentParser(
    description='Submit multiple runs to be pre-processed to LCLS batch farm')
parser.add_argument('start_run', type=int,
                    help='The first run to process')
parser.add_argument('end_run', type=int,
                    help='The last run to process')
parser.add_argument('-maxevt', type=int, default=None,
                    help='The maximum number of events to process for each run')

args = parser.parse_args()
start_run = args.start_run
end_run = args.end_run
maxevt = args.maxevt

# Submit the jobs
log_file_path_start = '../data/log/'
command_start = ''.join(['bsub -q ',
                         QUEUE_NAME,' ',
                         '-R "span[ptile='+str(PROCESSERS_PER_HOST)+']" ',
                         '-x -a mympi ',
                         '-n '+str(NUMBER_OF_CORES)+' '])
if maxevt is None:
    # No maximum number of events specified, process all events
    for curr_run in range(start_run, end_run+1):
        log_file_path = ''.join([log_file_path_start,
                                 str(curr_run),
                                 '_allevts.log'])
        command = ''.join([command_start,
                           '-o '+str(log_file_path)+' ',
                           '"python abs_XTC2HDF.py ',
                           str(curr_run)+'"'])
        os.system(command)
else:
    # Maximum event specifed, only process each run up to that event number
    for curr_run in range(start_run, end_run+1):
        log_file_path = ''.join([log_file_path_start,
                                 str(curr_run),
                                 '_'+str(maxevt)+'evts.log'])
        command = ''.join([command_start,
                           '-o '+str(log_file_path)+' ',
                           '"python abs_XTC2HDF.py ',
                           str(curr_run)+' ',
                           '-maxevt ',str(maxevt)+'"'])
        os.system(command)

print "All selected jobs are submitted to the batch farm"
