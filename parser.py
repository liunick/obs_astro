import pandas as pd
import numpy as np
import math, sys

base = 'https://api.astrocats.space/'
b_conditional = '/photometry/time+magnitude+e_magnitude+band?format=csv&band=B'
data = pd.read_csv('https://api.astrocats.space/catalog/claimedtype?claimedtype=Ia-(.*)&format=csv', error_bad_lines=False)

interpolation_constant = 5
continuity_constant = 20
number_of_readings = 0
time_accurate_sne = []

for sn in data["event"]:
    #Get data - convert columns to numerical values
    sn_data = pd.read_csv(base + sn + b_conditional)
    sn_data['magnitude'] = pd.to_numeric(sn_data['magnitude'])
    
    time_criteria = False
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
            number_of_readings += 1

    #All evaluated supernovae satisfy criteria #1
    if time_criteria:
        time_accurate_sne.append(sn)

print("Total SNe that satisfy light curve requirements: " + str(len(time_accurate_sne)))
print("List of SNe that satisfy light curves")
for sn in time_accurate_sne:
    print(sn)

    