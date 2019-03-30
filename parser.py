import pandas as pd
import numpy as np
import math, sys

base = 'https://api.astrocats.space/'
b_conditional = '/photometry/time+magnitude+e_magnitude+band?format=csv&band=B'
data = pd.read_csv('https://api.astrocats.space/catalog/claimedtype?claimedtype=Ia-(.*)&format=csv', error_bad_lines=False)

# constants
interpolation_constant = 5
continuity_constant = 20
phillips_delta = 15
b_band_error_margin = 0.3

band_magnitude = "b"
time_accurate_sne = []


def calc_phillips_expected_mag(band_mag_15):
    """
    Calculate expected b-band peak magnitude based on the band's magnitude after 15 days
    using the Phillips relationship

    band_mag_15 -- the 
    """


def print_sn_info(sn_name, b_mag_15):
    band_info = band_magnitude + " band ---- "
    snid = "Supernova ID: " + sn_name
    mag15 = band_magnitude + "_bang after 15D: " + b_mag_15
    print(band_info + snid + ", " + mag15)
    

def run():
    for sn in data["event"]:
        #Get data - convert columns to numerical values
        sn_data = pd.read_csv(base + sn + b_conditional)
        sn_data['magnitude'] = pd.to_numeric(sn_data['magnitude'])
        time_criteria = False
        delta_time_tracker = 0
        sne_obj = {}
        #check for criteria number #1
        if not sn_data["magnitude"].empty:

            # Find time at highest magnitude
            peak_index = sn_data["magnitude"].idxmin()
            peak_time = np.floor(sn_data.at[peak_index, 'time'])
            curr_time = peak_time

            # Check for 1 reading every 2 days and for 20 days of readings after peak
            for reading in sn_data["time"]:
                if reading < peak_time:
                    continue
                elif reading - interpolation_constant > curr_time:
                    break
                elif reading > peak_time + continuity_constant:
                    time_criteria = True
                    break
                curr_time = reading

                #check for closest to phillips_delta value
                # if reading - phillips_delta == 0

        #All evaluated supernovae satisfy criteria #1
        if time_criteria:
            time_accurate_sne.append(sn)

    print("Total SNe that satisfy light curve requirements: " + str(len(time_accurate_sne)))
    print("List of SNe that satisfy light curves")
    for sn in time_accurate_sne:
        print(sn)


if __name__ == '__main__':
    run()


    