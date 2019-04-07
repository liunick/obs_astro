import pandas as pd
import numpy as np
import math, sys


# CONSTANTS #
BASE = 'https://api.astrocats.space/'
B_CONDITIONAL = '/photometry/time+magnitude+e_magnitude+band?format=csv&band=B'
V_CONDITIONAL = '/photometry/time+magnitude+e_magnitude+band?format=csv&band=V'
I_CONDITIONAL = '/photometry/time+magnitude+e_magnitude+band?format=csv&band=I'
BAND_CONDITIONALS = [B_CONDITIONAL, V_CONDITIONAL, I_CONDITIONAL]
DATA = pd.read_csv('https://api.astrocats.space/catalog/lumdist+claimedtype?claimedtype=Ia-(.*)&format=csv', error_bad_lines=False)


B_BAND = 0
V_BAND = 1
I_BAND = 2
BANDS = [B_BAND, V_BAND, I_BAND]
BAND_LABELS = ["B", "V", "I"]

B_BAND_ERROR_MARGIN = 0.8
V_BAND_ERROR_MARGIN = 0.6
I_BAND_ERROR_MARGIN = 0.5
BAND_ERROR_MARGINS = [B_BAND_ERROR_MARGIN, V_BAND_ERROR_MARGIN, I_BAND_ERROR_MARGIN]

INTERPOLATION_CONSTANT = 5
CONTINUITY_CONSTANT = 20
PHILLIPS_DELTA = 15

def calc_absolute_magnitude(lumdist, app_mag):
    return app_mag - (5 * (math.log10(lumdist*1000000) - 1))

def calc_apparent_magnitude(lumdist, abs_mag):
    return abs_mag + (5 * (math.log10(lumdist*1000000) - 1))

def calc_phillips_expected_mag(lumdist, band_mag_15):
    """
    Calculate expected b-band peak absolute magnitude BASEd on the band's magnitude after 15 days
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
    
def check_criteria(sn, sn_data, band_type):
    criteria = [False, False]
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
            elif reading[0] - INTERPOLATION_CONSTANT > curr_time:
                break
            elif reading[0] > peak_time + CONTINUITY_CONSTANT:
                criteria[0] = True
                break
            curr_time = reading[0]
            
            #check for closest to PHILLIPS_DELTA value
            curr_delta = peak_time - (reading[0] - PHILLIPS_DELTA)
            neg_delta = (peak_time - (neg_delta_time_tracker[0] - PHILLIPS_DELTA))
            pos_delta = (peak_time - (pos_delta_time_tracker[0] - PHILLIPS_DELTA))
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
            mag15 = interpolate(neg_delta_time_tracker[0], neg_delta_time_tracker[1], pos_delta_time_tracker[0], pos_delta_time_tracker[1], peak_time + PHILLIPS_DELTA)
            expected_mag = calc_phillips_expected_mag(sn[1], mag15 - peak_mag)
            criteria[1] = within_margin(BAND_ERROR_MARGINS[band_type], peak_mag, expected_mag)
            print("Currently observing: " + str(sn[0]))
            print("\tInterpolated mag: " + str(mag15))
            print("\tPeak mag: " + str(peak_mag))
            print("\tLumdist: " + str(sn[1]))
            print("\tExpected mag: " + str(expected_mag))
            sys.stdout.flush()
    return criteria

def run():
    """
    Find list of supernova that fit light curve requirements
    """
    b_criteria_sne = []
    v_criteria_sne = []
    i_criteria_sne = []
    
    for sn in DATA[["event", "lumdist"]].values:
        for band in BANDS:
            #Get data - convert columns to numerical values
            sn_data = pd.read_csv(BASE + sn[0] + BAND_CONDITIONALS[band])
            sn_data['magnitude'] = pd.to_numeric(sn_data['magnitude'])
            time_criteria, margin_criteria = check_criteria(sn, sn_data, band)

            #All evaluated supernovae satisfy criteria #1
            if time_criteria:
                if band == B_BAND:
                    b_criteria_sne.append([sn[0], margin_criteria])
                elif band == V_BAND:
                    v_criteria_sne.append([sn[0], margin_criteria])
                elif band == I_BAND:
                    i_criteria_sne.append([sn[0], margin_criteria])
    
    criteria_sne = [b_criteria_sne, v_criteria_sne, i_criteria_sne]
    for band in BANDS:
        within_margin_ratio = 0
        # print("List of SNe that satisfy light curves")
        for sn in criteria_sne[band]:
            if sn[1]:
                within_margin_ratio += 1
        total_number_of_sne = len(criteria_sne[band])
        print("--------------------------------------------------------")    
        print("Total SNe that satisfy light curve requirements for " + BAND_LABELS[band] + "-band: " + str(total_number_of_sne))
        print("Percentage of SNe that satisfy" + BAND_LABELS[band] + "-band error on expected peak mag: " + str(within_margin_ratio/total_number_of_sne))

if __name__ == '__main__':
    run()


    