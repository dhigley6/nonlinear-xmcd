"""Classes to help with data extraction from LCLS XTC files

Daniel Higley - 09/2015 -- 01/2016
"""

import numpy as np
import psana

class BasePS(object):
    """A base class to inherit from
    """

    def __init__(self):
        self.reset_data()
        self.data_names = ['']

    def reset_data(self):
        """Reset data and accumulation statistics
        """
        self.data = []
        self.ngood = 0
        self.nbad = 0
        self.nrecord = 0

    def _get_evt_det(self, evt):
        """Return psana detector object for the event
        """
        raise NotImplementedError

    def check_evt(self, evt):
        """Check whether the event is adequate
        """
        evt_det = self._get_evt_det(evt)
        if evt_det is None:
            self.nbad += 1
            if self.nbad%20 == 0:
                print ''.join([str(self.nbad),
                               ' no ', self.dataset_name,
                               ' events so far on this rank'])
            return False
        else:
            self.ngood += 1
            return True

    def process_evt(self):
        """To be executed on each new, good event
        """
        raise NotImplementedError

    def gather_data(self, comm):
        """Gather data between different cores
        """
        self.data = comm.gather(self.data)
        if self.data is not None:
            self.data = [item for sublist in self.data for item in sublist]

    def get_data(self):
        """Return the accumulated data as a numpy array
        """
        data_array = np.array(self.data, self.data_dtype)
        return data_array

class PSGetEPICS(BasePS):
    """Extract specified EPICS PVs
    """

    def __init__(self, pv_list, epics):
        BasePS.__init__(self)
        self.pv_list = pv_list
        self.epics = epics
        self.dataset_name = 'Epics'
        self.data_dtype = None

    def check_evt(self, evt):
        epics_data = [self.epics.value(pv) for pv in self.pv_list]
        if None in epics_data:
            self.nbad += 1
            if self.nbad%20 == 0:
                print ''.join([str(self.nbad),
                               ' no EPICS events so far on this rank'])
            return False
        else:
            self.ngood += 1
            return True

    def process_evt(self, evt):
        evt_data = [self.epics.value(pv) for pv in self.pv_list]
        if self.data_dtype is None:
            # If the data type has not been set yet, set it
            rec_dtype = np.dtype({'names': self.pv_list,
                                  'formats': [(np.array(datum).dtype,
                                               np.size(np.array(datum)))
                                              for datum in evt_data]})
            self.data_dtype = rec_dtype
        self.data.append(tuple(evt_data))
        self.nrecord += 1
            

class PSGetAllData(BasePS):
    """Extract all data from a specified psana object (detector)
    IMPORTANT: intended to only be used for scalar data
    """

    def __init__(self, dataset_name, det_type, det_source=None):
        BasePS.__init__(self)
        self.dataset_name = dataset_name
        self.det_type = det_type
        self.det_source = det_source
        self.data_dtype = None

    def _get_evt_det(self, evt):
        if self.det_source is None:
            det_data = evt.get(self.det_type)
        else:
            det_data = evt.get(self.det_type,
                               psana.Source(self.det_source))
        return det_data

    def process_evt(self, evt):
        # Value is taken out from the methods, because it makes no sense
        # it is an attribute of the delay encoder
        det = self._get_evt_det(evt)  
        det_attrs = filter(lambda a: ((not a.startswith('__'))
                                      and (a != 'value')
                                      and (a != 'DamageMask')
                                      and (a != 'idxtime')), dir(det))
        det_methods = filter(lambda a: callable(det.__getattribute__(a)),
                             det_attrs)
        det_params = filter(lambda a:
                            not callable(det.__getattribute__(a)),
                            det_attrs)
        det_attr_names = np.append(det_params, det_methods)
        evt_data = [det.__getattribute__(det_param)
                    for det_param in det_params]
        evt_data.extend([det.__getattribute__(det_method)()
                         for det_method in det_methods])
        if self.data_dtype is None:
            # If the data type has not been set yet, set it
            rec_dtype = np.dtype({'names': det_attr_names,
                                  'formats': [(np.array(datum).dtype,
                                               np.size(np.array(datum)))
                                              for datum in evt_data]})
            self.data_dtype = rec_dtype
        self.data.append(tuple(evt_data))
        self.nrecord += 1

class PSGetEVR(BasePS):
    """Extract event codes
    """

    def __init__(self):
        self.ngood = 0
        self.nbad = 0
        self.evr_src = psana.Source('DetInfo(NoDetector.0:Evr.0)')
        self.data = []
        self.dataset_name = 'EVR'

    def _get_evt_det(self, evt):
        return evt.get(psana.EvrData.DataV4, self.evr_src)

    def process_evt(self, evt):
        evr_data = self._get_evt_det(evt)
        evt_event_codes = [event_code.eventCode()
                           for event_code in evr_data.fifoEvents()]
        self.data.append(evt_event_codes)

    def get_data(self):
        # Flatten the list of lists of event codes to a single list:
        #event_codes_all = [evt_codes
        #                   for rank_evt_codes in self.data
        #                   for evt_codes in rank_evt_codes]
        event_codes_all = self.data
        # Get all the unique event codes present in the data
        event_codes_sets = [set(evt_event_codes)
                            for evt_event_codes in event_codes_all]
        event_codes_present = set.union(*event_codes_sets)
        # Make the array
        event_code_names = [str(item) for item in event_codes_present]
        arr_dtype = np.dtype({'names': event_code_names,
                              'formats': [bool]*len(event_code_names)})
        arr = np.recarray(len(event_codes_all), arr_dtype)
        for code_num in event_codes_present:
            arr[str(code_num)] = [code_num in evt_event_codes
                                  for evt_event_codes in event_codes_all]
        return arr

class PSAcqInt(BasePS):
    """Extract integral of an acqiris trace
    """
    # TODO: Upgrade this to do multiple acqiris integrals

    def __init__(self, acq_num, acq_chan,
                 acq_dead_start, acq_dead_end,
                 acq_int_start, acq_int_end,
                 invert=True):
        BasePS.__init__(self)
        self.data_dtype = np.dtype({'names':['acq'], 'formats': [float]})
        self.acq_num = acq_num
        self.acq_chan = acq_chan
        self.dataset_name = ''.join(['Acqiris',str(acq_num)])
        self.dead_start = int(acq_dead_start)
        self.dead_end = int(acq_dead_end)
        self.int_start = int(acq_int_start)
        self.int_end = int(acq_int_end)
        self.invert = invert

    def _get_evt_det(self, evt):
        if self.acq_num == 1:
            acq1_src = psana.Source('DetInfo(SxrEndstation.0:Acqiris.1)')
            acq_det = evt.get(psana.Acqiris.DataDescV1,
                              acq1_src)
        elif self.acq_num == 2:
            acq2_src = psana.Source('DetInfo(SxrEndstation.0:Acqiris.2)')
            acq_det = evt.get(psana.Acqiris.DataDescV1,
                              acq2_src)
        return acq_det

    def process_evt(self, evt):
        acq_det = self._get_evt_det(evt)
        acq_waveform = acq_det.data(self.acq_chan).waveforms()[0]
        # Integrate signal and subtract background
        acq_dead = acq_waveform[self.dead_start:self.dead_end]
        acq_sig = acq_waveform[self.int_start:self.int_end]
        background_per_point = np.median(acq_dead)
        signal = np.sum(acq_sig-background_per_point)
        if self.invert:
            signal = -1*signal
        self.data.append((signal,))
        self.nrecord += 1

class PSMethod(BasePS):
    """Extract data from a method of a psana detector
    """

    def __init__(self, dataset_name,
                 det_type, det_source,
                 det_method):
        BasePS.__init__(self)
        self.dataset_name = dataset_name
        self.det_type = det_type
        self.det_source = det_source
        self.det_method = det_method

    def _get_evt_det(self, evt):
        det = evt.get(self.det_type, psana.Source(self.det_source))
        return det

    def get_data(self):
        data_dtype = np.dtype({'names':[self.det_method], 'formats': [type(self.data[0][0])]})
        data_array = np.array(self.data, data_dtype)
        return data_array

    def process_evt(self, evt):
        det = self._get_evt_det(evt)
        evt_data = det.__getattribute__(self.det_method)()
        self.data.append((evt_data,))
        self.nrecord += 1

class PSCCDROIs(BasePS):
    """Extract ROIs from calibrated CCD images
    """

    def __init__(self, src, key, roi_list, roi_names,
                 dataset_name='CCD'):
        BasePS.__init__(self)
        self.data_dtype = np.dtype({'names':roi_names, 'formats': [float]*len(roi_names)})
        self.src = src
        self.key = key
        self.roi_list = roi_list
        self.dataset_name = dataset_name

    def _get_evt_data(self, evt):
        ccd_src = psana.Source(self.src)
        ccd_data = np.copy(evt.get(psana.ndarray_float64_2,
                                   ccd_src, self.key))
        return ccd_data

    def check_evt(self, evt):
        """Check whether the event is adequate
        """
        evt_data = self._get_evt_data(evt)
        if (evt_data is None) or (np.size(evt_data) == 1):
            self.nbad += 1
            if self.nbad%20 == 0:
                print ''.join([str(self.nbad),
                               ' no ', self.dataset_name,
                               ' events so far on this rank'])
            return False
        else:
            self.ngood += 1
            return True

    def process_evt(self, evt):
        data = self._get_evt_data(evt)
        roi_sums = [np.sum(data[int(roi[0]):int(roi[1]),
                                int(roi[2]):int(roi[3])])
                         for roi in self.roi_list]
        self.data.append(tuple(roi_sums))
        self.nrecord += 1

class PSANDORROIs(BasePS):
    """Extract ROIs from RAW ANDOR images
    """

    def __init__(self, alias, roi_list, roi_names,
                 psana_env, dataset_name='Andor'):
        BasePS.__init__(self)
        self.data_dtype = np.dtype({'names':roi_names, 'formats': [float]*len(roi_names)})
        self.alias = alias
        self.roi_list = roi_list
        self.dataset_name = dataset_name
        self.detector = None
        self.andor = psana.Detector(self.alias, psana_env)

    def _get_evt_data(self, evt): 
        ccd_data = np.copy(self.andor.raw(evt))
        return ccd_data

    def check_evt(self, evt):
        """Check whether the event is adequate
        """
        evt_data = self._get_evt_data(evt)
        if (evt_data is None) or (np.size(evt_data) == 1):
            self.nbad += 1
            if self.nbad%20 == 0:
                print ''.join([str(self.nbad),
                               ' no ', self.dataset_name,
                               ' events so far on this rank'])
            return False
        else:
            self.ngood += 1
            return True

    def process_evt(self, evt):
        data = self._get_evt_data(evt)
        if np.size(np.shape(data)) == 2:
            roi_sums = [np.mean(data[int(roi[0]):int(roi[1]),
                                    int(roi[2]):int(roi[3])])
                        for roi in self.roi_list]
        else:
            roi_sums = [np.mean(data[int(roi[2]):int(roi[3])])
                        for roi in self.roi_list]
        self.data.append(tuple(roi_sums))
        self.nrecord += 1
