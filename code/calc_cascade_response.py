import numpy as np
import matplotlib.pyplot as plt

def calculate_cascade_time(percent_within_2eV):
    cascade_durations = np.arange(0.0001, 40, 0.1)
    integrated_responses = list(map(
        calc_cascade_response, 
        cascade_durations))
    response_at_or_lower_than_calculation = integrated_responses < percent_within_2eV
    cascade_duration = cascade_durations[response_at_or_lower_than_calculation][0]
    print('Caclulated fraction of energy within 2 eV of Fermi level of ')
    print(percent_within_2eV)
    print('Calculated cascasde duration of')
    print(''.join([str(round(cascade_duration)), ' fs']))
    return cascade_duration

def calc_cascade_response(cascade_duration=20):
    """Return expected energy within 2 eV of Fermi level
    Uses model described in manuscript
    """
    time = np.arange(-1000, 2000, 0.1)
    pulse = (time > 0) & (time <= 40)
    cascade_time = np.arange(-400, 400, 0.1)
    cascade = np.ones_like(cascade_time)
    cascade[cascade_time < 0] = 0
    cascade[cascade_time > cascade_duration] = 0
    cascade = cascade/np.sum(cascade)
    response = np.convolve(pulse, cascade, mode='same')
    integrated_response = np.sum(response*pulse)/np.sum(pulse)
    return integrated_response

def test_plot():
    cascade_durations = np.arange(0.0001, 100, 10)
    integrated_responses = list(map(calc_cascade_response, cascade_durations))
    plt.figure()
    plt.plot(cascade_durations, integrated_responses)
    plt.xlabel('Cascade Duration (fs)')
    plt.ylabel('Integrated Response')
