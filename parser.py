import pandas as pd
import numpy as np
import math, sys

base = 'https://api.astrocats.space/'
b_conditional = '/photometry/time+magnitude+e_magnitude+band?format=csv&band=B'
data = pd.read_csv('https://api.astrocats.space/catalog/lumdist+claimedtype?claimedtype=Ia-(.*)&format=csv', error_bad_lines=False)

# constants
interpolation_constant = 5
continuity_constant = 20
phillips_delta = 15
b_band_error_margin = 0.8

band_magnitude = "b"

def calc_absolute_magnitude(lumdist, app_mag):
    return app_mag - (5 * (math.log10(lumdist*1000000) - 1))

def calc_apparent_magnitude(lumdist, abs_mag):
    return abs_mag + (5 * (math.log10(lumdist*1000000) - 1))

def calc_phillips_expected_mag(lumdist, band_mag_15):
    """
    Calculate expected b-band peak absolute magnitude based on the band's magnitude after 15 days
    using the Phillips relationship

    band_mag_15 -- the band's magnitude 15 days after the peak 
    """
    a = -21.726
    b = 2.698
    expected_mag = a + b * band_mag_15
    return calc_apparent_magnitude(lumdist, expected_mag)

def within_margin(margin, actual, expected):
    return actual < expected + margin and actual > expected - margin

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
    for sn in data[["event", "lumdist"]].values:

        #Get data - convert columns to numerical values
        sn_data = pd.read_csv(base + sn[0] + b_conditional)
        sn_data['magnitude'] = pd.to_numeric(sn_data['magnitude'])
        time_criteria = False
        margin_criteria = False

        #check for criteria number #1
        if not sn_data["magnitude"].empty:
            neg_delta_time_tracker = [-float("inf"), None]
            pos_delta_time_tracker = [float("inf"), None]
            # Find time at highest magnitude
            peak_mag = sn_data["magnitude"].min()
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

            #linearly interpolate between 2 points to reach mag15
            if (
                neg_delta_time_tracker[0] is not None and
                neg_delta_time_tracker[1] is not None and
                pos_delta_time_tracker[0] is not None and
                pos_delta_time_tracker[1] is not None
            ):
                
                mag15 = interpolate(neg_delta_time_tracker[0], neg_delta_time_tracker[1], pos_delta_time_tracker[0], pos_delta_time_tracker[1], peak_time + phillips_delta)
                expected_mag = calc_phillips_expected_mag(sn[1], mag15 - peak_mag)
                margin_criteria = within_margin(b_band_error_margin, peak_mag, expected_mag)
                print("Currently observing: " + str(sn[0]))
                print("\tInterpolated mag: " + str(mag15))
                print("\tPeak mag: " + str(peak_mag))
                print("\tLumdist: " + str(sn[1]))
                print("\tExpected mag: " + str(expected_mag))
                sys.stdout.flush()


        #All evaluated supernovae satisfy criteria #1
        if time_criteria:
            time_accurate_sne.append([sn[0], margin_criteria])
            
    within_margin_ratio = 0
    print("List of SNe that satisfy light curves")
    for sn in time_accurate_sne:
        if sn[1]:
            within_margin_ratio += 1
        print(sn[0])
    total_number_of_sne = len(time_accurate_sne)
    print("--------------------------------------------------------")    
    print("Total SNe that satisfy light curve requirements: " + str(total_number_of_sne))
    print("Percentage of SNe that satisfy B-band error on expected peak mag: " + str(within_margin_ratio/total_number_of_sne))
    return time_accurate_sne

def run():
    time_accurate_sne = filter_time_accurate()

if __name__ == '__main__':
    run()


    