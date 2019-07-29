"""A pre-processor of data from 2016 LK30 nonlinear absorption exp.

Usage:
    abs_XTC2HDF.py <run> [<maxevt>] [-h]

Options:
    <run>        Run Number
    <maxevt>     Maximum number of events to process
    -h|--help    See help information

## Important Note:
## Sometimes a particular dataset will be non-existant in a set of runs. In
## that case, it needs to be disabled when pre-processing the XTC data for
## that run set. To do this, change the det_ext_list variable below. More
## information is provided where this variable is defined.
"""

print "IMPORTING PYTHON MODULES..."
import numpy as np
import h5py
import argparse
import time
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
num_ranks = comm.Get_size()
print "DONE"

# LCLS psana to read data
print "IMPORTING psana..."
import psana
print "DONE"

# psExtMod to handle data from the LCLS detectors
print "IMPORTING PSExtMod..."
import psExtMod
print "DONE"

#### COMMAND LINE ARGUMENTS ###################################################
exp_desc = 'Pre-process data from 2016 LCLS nonlinear absorption experiment'
parser = argparse.ArgumentParser(description=exp_desc)
parser.add_argument("run", type=int,
                    help='The run number to process')
parser.add_argument('-maxevt', type=int, default=None,
                    help=''.join(['The maximum number of events to process. ',
                                  'If not provided, processes all events']))
args = parser.parse_args()
run_num = args.run
maxevt = args.maxevt
#### END COMMAND LINE ARGUMENTS ###############################################

#### GLOBAL CALIBRATIONS ######################################################
PSANA_DS = 'exp=SXR/sxrk3016:run='
DATA_SAVE_START = '../data/pre_proc/'
#### END GLOBAL CALIBRATIONS ##################################################

#### RUN-SPECIFIC CALIBRATIONS ################################################
# path to csv file to load run specific calibrations from
RUN_SPECIFIC_CALS_FILE = '../cals/abs_run_reference.csv'
# List of run specific calibration coefficients to make sure are present:
RUN_SPECIFIC_CALS_NAMES = ['run_start', 'run_end',
                           'mcp_int_start', 'mcp_int_end',
                           'mcp_dead_start', 'mcp_dead_end',
                           'andor_inty_start', 'andor_inty_end',
                           'andor_intx_start', 'andor_intx_end',
                           'ref_andor_inty_start', 'ref_andor_inty_end',
                           'ref_andor_intx_start', 'ref_andor_intx_end']

# Load run-specific calibrations
print 'Loading run-dependent constants for analysis...'
run_specific_cals = np.genfromtxt(RUN_SPECIFIC_CALS_FILE,
                                  names=True, delimiter=',')
# check that all the necessary calibration coefficient names are present
cal_names_found = run_specific_cals.dtype.names
cal_names_present = [cal_name_req in cal_names_found
                     for cal_name_req in RUN_SPECIFIC_CALS_NAMES]
if not all(cal_names_present):
    raise ValueError("Not all required calibration coefficents found")
run_start_col = run_specific_cals['run_start']
run_end_col = run_specific_cals['run_end']
# Find the specified run number calibrations as the row which has a starting
# number less than or equal to the run number of interest and a ending
# number greater than or equal to the run number of interest (assumes that a
# given run will not fulfill [and should not fulfill!] these conditions for
# more than one row).
rows_greater_or_equal = run_num >= run_start_col
rows_less_or_equal = run_num <= run_end_col
row_num_want = rows_greater_or_equal & rows_less_or_equal
run_cals = run_specific_cals[row_num_want]
print "Found Run Specific Calibrations of"
print run_cals.dtype.names
print run_cals
print "DONE"
#### END RUN-SPECIFIC CALIBRATIONS ############################################

#### MAIN EVENT LOOP PREPARATIONS #############################################
print ''.join(['<<< Start Analysis of run ',
               str(run_num),
               ' of nonlinear absorption experiment, sxrk3016 >>>'])
if (maxevt is None):
    print 'Processing all events'
else:
    print 'Processing the first ' + str(maxevt) + ' events'

if rank == 0:
    start_time = time.time()
# PSANA Config File
psana.setConfigFile('psanacfg.cfg')
# PSANA Data Source
dsStr = PSANA_DS + str(run_num)
ds = psana.DataSource(dsStr+':idx')
run = ds.runs().next()
epics = ds.env().epicsStore()
times = run.times()
# select a subset of events by their times,
# so each cpu-core ("rank") works on a separate set of events
evts_per_rank = len(times)/num_ranks
mytimes = times[rank*evts_per_rank:(rank+1)*evts_per_rank]
if maxevt is not None:
    # If a maximum event number is specified in command line, calculate
    # the number of events each rank should process to reach this maximum
    # number of events (if there are enough events, otherwise default to
    # just process all events).
    maxevt_per_rank = maxevt/num_ranks
    evts_per_rank = min(maxevt_per_rank, evts_per_rank)
print 'Calculated times of events for each rank to process'
#### SETUP DETECTORS ##########################################################
# This section sets up the detectors (data sets) to extract data from the XTC
# file for.
det_ext_list = []
ext_evtid = psExtMod.PSGetAllData(dataset_name='evtid',
                                  det_type=psana.EventId)
ext_ebeam = psExtMod.PSGetAllData(dataset_name='ebeam',
                                  det_type=psana.Bld.BldDataEBeamV7,
                                  det_source='BldInfo(EBeam)')
gd_source = 'BldInfo(FEEGasDetEnergy)'
ext_gd = psExtMod.PSGetAllData(dataset_name='gd',
                               det_type=psana.Bld.BldDataFEEGasDetEnergyV1,
                               det_source=gd_source)
ext_gmd = psExtMod.PSGetAllData(dataset_name='gmd',
                                det_type=psana.Bld.BldDataGMDV2,
                                det_source='BldInfo(GMD)')
monoenc_src = 'DetInfo(SxrEndstation.0:USDUSB.0)'
ext_monoenc = psExtMod.PSGetAllData(dataset_name='monoenc',
                                    det_type=psana.UsdUsb.DataV1,
                                    det_source=monoenc_src)
ext_evr = psExtMod.PSGetEVR()
# List of EPICS PVs to extract
EPICS_PVS = ['SXR:MON:MMS:06.RBV', 'SXR:GMD:SRG:01:Calib:Pressure:Calc',
             'SXR:EXP:AIN:1',
             'SXR:EXP:MMS:01.RBV', 'SXR:EXP:MMS:02.RBV',
             'SXR:EXP:MMS:03.RBV',
             'SXR:MNT:MMS:01.RBV', 'SXR:MNT:MMS:02.RBV', 
             'SXR:FLX:MMS:01.RBV', 'SXR:EXP:MMS:12.RBV']
#EPICS_PVS = ['SXR:MON:MMS:06.RBV']
# Description of EPICS PVs extracted:
# SXR:MON:MMS:06.RBV: Position of the mono, corresponds to photon energy (but 
# use encoded values where possible, since this PV updates slowly)
# SXR:GMD:SRG:01:Calib:Pressure:Calc: GMD presure read out by the spinning
# rotor gauge
# SXR:EXP:AIN:1      : Magnet voltage readback
# SXR:EXP:MMS:01.RBV: x of RSXS manipulator
# SXR:EXP:MMS:02.RBV: y of RSXS manipulator
# SXR:EXP:MMS:03.RBV: z of RSXS manipulator
# SXR:MNT:MMS:01.RBV: y of upstream filter translator
# SXR:MNT:MMS:02.RBV: y of downstream filter translator
# SXR:FLX:MMS:01.RBV: y of FLX filter manipulator
# SXR:EXP:MMS:12.RBV: y of filter manipulator in front of RCI
ext_epics = psExtMod.PSGetEPICS(EPICS_PVS, epics)
ext_mcp = psExtMod.PSAcqInt(acq_num=2,
                            acq_chan=0,
                            acq_dead_start=run_cals['mcp_dead_start'],
                            acq_dead_end=run_cals['mcp_dead_end'],
                            acq_int_start=run_cals['mcp_int_start'],
                            acq_int_end=run_cals['mcp_int_end'])
sig_roi = [run_cals['andor_inty_start'], run_cals['andor_inty_end'],
           run_cals['andor_intx_start'], run_cals['andor_intx_end']]
ref_roi = [run_cals['ref_andor_inty_start'], run_cals['ref_andor_inty_end'],
           run_cals['ref_andor_intx_start'], run_cals['ref_andor_intx_end']]
andor_roi_list = [sig_roi, ref_roi]
andor_roi_names = ['signal', 'reference']
ext_andor = psExtMod.PSANDORROIs(alias='andor',
                                 roi_list=andor_roi_list,
                                 roi_names=andor_roi_names,
                                 psana_env=ds.env())

# Select detectors to extract data from (if, for some reason, one does not
# want to extract the data from a particular data set or that data is
# non-existant for a particular run, simply temporarily comment out the
# appropriate line below to disable it [Be sure to uncomment it again when
# processing other runs])
det_ext_list = []
det_ext_list.append(ext_evtid)
det_ext_list.append(ext_ebeam)
det_ext_list.append(ext_gd)
#det_ext_list.append(ext_gmd)
det_ext_list.append(ext_monoenc)
det_ext_list.append(ext_evr)
det_ext_list.append(ext_epics)
det_ext_list.append(ext_mcp)
det_ext_list.append(ext_andor)
#### END SETUP DETECTORS ######################################################
#### END MAIN EVENT LOOP PREPARATIONS #########################################

#### MAIN EVENT LOOP ##########################################################
ngood = 0      # Counter for number of good events processed
print '#######################################################################'
print 'Start Main Event Loop on rank ', rank
for eventCounter, curr_time in enumerate(mytimes):
    # grab the current event
    evt = run.event(curr_time)
    # Every 500 events, print out the number of good/bad events
    if ((eventCounter)%500 == 0):
        print
        print 'Rank ', rank, ' Event Counter:'
        print 'Events: ', eventCounter, ' Good: ', ngood
        print
    # Check the event
    evt_checks = [det_data.check_evt(evt) for det_data in det_ext_list]
    if np.all(evt_checks) == False:
        # This event is bad, go to the next event
        continue
    # Bad event checks have been performed, this is a good event
    ngood += 1
    # Process the event's data
    for det_ext in det_ext_list:
        det_ext.process_evt(evt)
    # Break out of the loop if we have gone through the maximum specifed
    # number of events per core
    if maxevt is not None and eventCounter >= (evts_per_rank-1):
        break
##### END MAIN EVENT LOOP #####################################################

##### POSTPROCESS/SAVE DATA ###################################################
# Gather all the data from the different ranks:
for det_ext in det_ext_list:
    det_ext.gather_data(comm)
# Save data from each detector:
if rank == 0:
    print "Saving data..."
    if maxevt is None:
        save_file_path = DATA_SAVE_START+'run'+str(run_num)+'allevts.h5'
    else:
        save_file_path = DATA_SAVE_START+'run'+str(run_num)+'first'+str(maxevt)+'evts.h5'
    save_file = h5py.File(save_file_path, 'w')
    for det_ext in det_ext_list:
        save_file.create_dataset(det_ext.dataset_name,
                                 data=det_ext.get_data())
    save_file.close()
    end_time = time.time()
    print 'Done with everything'
    print 'Elapsed time of ', end_time-start_time, ' seconds'
#### END POSTPROCESS/SAVE DATA ################################################
