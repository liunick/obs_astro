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



def calc_phillips_expected_mag(band_mag_15):
    """
    Calculate expected b-band peak magnitude based on the band's magnitude after 15 days
    using the Phillips relationship

    band_mag_15 -- the band's magnitude 15 days after the peak 
    """

def interpolate(date1, mag1, date2, mag2, date15):
    """
    Calculate interpolated magnitude for 15 days after peak

    date1   -- date of first mag reading
    mag1    -- value of first mag reading
    date2   -- date of second mag reading
    mag2    -- value of second mag reading
    date15  -- exact date of 15 days after peak
    """
    slope = (mag2 - mag1) / (date2 - date1)
    days_to_15 = date15 - date1
    mag_delta = slope * days_to_15
    return mag1 + mag_delta

def print_sn_info(sn_name, b_mag_15):
    band_info = band_magnitude + " band ---- "
    snid = "Supernova ID: " + sn_name
    mag15 = band_magnitude + "_bang after 15D: " + b_mag_15
    print(band_info + snid + ", " + mag15)
    

def filter_time_accurate():
    """
    Find list of supernova that fit light curve requirements
    """
    time_accurate_sne = []
    for sn in data["event"]:
        #Get data - convert columns to numerical values
        sn_data = pd.read_csv(base + sn + b_conditional)
        sn_data['magnitude'] = pd.to_numeric(sn_data['magnitude'])
        time_criteria = False
        
        sne_obj = {}
        
        print("Currently looking at: " + str(sn))
        sys.stdout.flush()

        #check for criteria number #1
        if not sn_data["magnitude"].empty:
            neg_delta_time_tracker = [-float("inf"), None]
            pos_delta_time_tracker = [float("inf"), None]
            # Find time at highest magnitude
            peak_index = sn_data["magnitude"].idxmin()
            peak_time = np.floor(sn_data.at[peak_index, 'time'])
            curr_time = peak_time
            # Check for 1 reading every 2 days and for 20 days of readings after peak
            for reading in sn_data[["time", "magnitude"]].values:
                if reading[0] < peak_time:
                    continue
                elif reading[0] - interpolation_constant > curr_time:
                    break
                elif reading[0] > peak_time + continuity_constant:
                    time_criteria = True
                    break
                curr_time = reading[0]
                
                #check for closest to phillips_delta value
                curr_delta = peak_time - (reading[0] - phillips_delta)
                neg_delta = (peak_time - (neg_delta_time_tracker[0] - phillips_delta))
                pos_delta = (peak_time - (pos_delta_time_tracker[0] - phillips_delta))
                if curr_delta > 0 and curr_delta < neg_delta: 
                    neg_delta_time_tracker = [reading[0], reading[1]]
                elif curr_delta < 0 and curr_delta > pos_delta:
                    pos_delta_time_tracker = [reading[0], reading[1]]

            # print("\tCurr neg: [" + str(neg_delta_time_tracker[0]) + ", " + str(neg_delta_time_tracker[1]) + "]")
            # print("\tCurr pos: [" + str(pos_delta_time_tracker[0]) + ", " + str(pos_delta_time_tracker[1]) + "]")
            # sys.stdout.flush()

            #linearly interpolate between 2 points to reach mag15
            if (
                neg_delta_time_tracker[0] is not None and
                neg_delta_time_tracker[1] is not None and
                pos_delta_time_tracker[0] is not None and
                pos_delta_time_tracker[1] is not None
            ):

                mag15 = interpolate(neg_delta_time_tracker[0], neg_delta_time_tracker[1], pos_delta_time_tracker[0], pos_delta_time_tracker[1], peak_time + phillips_delta)
                print("\tInterpolated mag: " + str(mag15))
                sys.stdout.flush()


        #All evaluated supernovae satisfy criteria #1
        if time_criteria:
            time_accurate_sne.append(sn)
            
    print("Total SNe that satisfy light curve requirements: " + str(len(time_accurate_sne)))
    print("List of SNe that satisfy light curves")
    for sn in time_accurate_sne:
        print(sn)    
    return time_accurate_sne

def run():
    time_accurate_sne = filter_time_accurate()

if __name__ == '__main__':
    run()


    