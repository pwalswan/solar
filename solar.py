from datetime import datetime, timedelta
from math import sin, cos
import pandas as pd
import matplotlib.pyplot as plt

def processPlan(planName, planStartTime, planEndTime):
    print('processPlan ' + planName + ' ...')
    # start 2009-11-28 03:17:30.850
    # end   2009-11-28 07:45:27.870

    # At the Downscan start time, the grating position in instrumentTelemetry goes from 0 to non-zero: 107189.77591053583
    # instrumentTelementry line 29159
    # At the Downscan end time, the grating position goes back to zero
    # instrumentTelementry line 31696

    # UpScan starts on
    # instrumentTelemetry line 38938
    # UpScan end on 
    # instrumentTelementry eof

    # Loop through instrumentTelemetry since it has the first factor we need
    instrumentTelemetry = pd.read_csv('./data/instrumentTelemetry.txt', skipinitialspace=True)
    detectorTemps = pd.read_csv('./data/detectorTemp.txt', skipinitialspace=True)
    referenceSpectrum = pd.read_csv('./data/referenceSpectrum.txt', skipinitialspace=True)

    # New list and dataframe for results
    resultsList = []

    # Get the start and end index for the plan matched on microseconds
    startIndex = instrumentTelemetry[instrumentTelemetry['microsecondsSinceGpsEpoch'] == planStartTime].index[0]
    try:
        endIndex = instrumentTelemetry[instrumentTelemetry['microsecondsSinceGpsEpoch'] == planEndTime].index[0]
    except:
        endIndex = instrumentTelemetry.index[-1] # Get the last index if the planEndTime exceeds the last intrumentTelemetry reading

    i = 1

    for telemetryRowIndex in range(startIndex, endIndex):
        instrumentTelementryRow = instrumentTelemetry.iloc[telemetryRowIndex]
        gratingPosition = instrumentTelementryRow['gratPos']

        # Wavelength (the grating equation)
        offset = 239532.38
        stepSize = 2.4237772022101214E-6 # [rad/step]
        d = 277.77777777777777 # [nm]
        phiGInRads = 0.08503244115716374 # [rad]
        ang1 = (offset - gratingPosition) * stepSize # [rad]
        wavelength = 2 * d * sin(ang1) * cos(phiGInRads / 2.0) # [nm]

        count_rate = instrumentTelementryRow['counts']

        # Get the detectorTemp. All the data provided is time-tagged. Using the telemetry index for the detectorTemps file
        detectorTempsRow = detectorTemps.iloc[telemetryRowIndex]
        detectorTemp = detectorTempsRow['temp (C)']

        # Photon Count Rate (counts/s/m2/nm)
        tempCorrFactor = 0.0061628 # [counts/degC]
        count_rate_corr = count_rate * (1.0 + tempCorrFactor * (20.0 - detectorTemp))
        # This is looping through an UpScan or a DownScan, so don't expect to see a Dark count here
        # Using the first intrumentTelemetry count for Dark plan in plans.txt.
        # Hard coding until resolved
        dark_counts = 504.61167548376125
        dark_integrationTime = 9.434192873600002E14
        dark_count_rate = dark_counts / dark_integrationTime
        dark_count_rate_corr = dark_count_rate * (1.0 + tempCorrFactor * (20.0 - detectorTemp))
        # Seems dark_count_rate_corr should be a collection?
        # median_dark_count_rate = median(dark_count_rate_corr)
        median_dark_count_rate = dark_count_rate_corr # Just set it to the dark_count_rate_corr
        apertureArea = .01 / (1E2 * 1E2) # [m^2] (aperature area from cm^2 to m^2)
        photonsPerSecondPerM2 = (count_rate_corr - median_dark_count_rate) / apertureArea # [photons/sec/m^2/nm]

        h = 6.62606957E-34 # [J*s]
        c = 299792458.0 # [m/s]
        wavelengthInMeters = wavelength / 10**9
        energyPerPhoton = h * c / wavelengthInMeters # [J]
        wattsPerM2 = photonsPerSecondPerM2 * energyPerPhoton # [watts/m^2/nm]

        # Look up wavelength by irradiance
        # My exact wattsPerM2 don't match to the 19 decimals, so try to get a match using 6.
        truncatedWattsPerM2 = int(wattsPerM2 * 10**6) / 10**6 # Truncate to six decimals
        truncatedWattsPerM2_str = '{:.6f}'.format(truncatedWattsPerM2) # convert to strings for trunc comparison (not a round)
        referenceSpectrum['irradiance_str'] = referenceSpectrum['irradiance (watts/m^2/nm)'].apply(lambda x: '{:.6f}'.format(x)) # trunc 6
        try: 
            referenceIndex = referenceSpectrum[referenceSpectrum['irradiance_str'] == truncatedWattsPerM2_str].index[0] # create string copy of irradiance
        except:
            # print('No matching irradiance in referenceSpectrum. Likely due to the dark counts not being included.')
            continue

        # We have the index, there may be more than one, or they may not be a match due to truncation
        referenceSpectrumRow = referenceSpectrum.iloc[referenceIndex]
        referenceSpectrumWavelength = referenceSpectrumRow['wavelength(nm)']
        
        resultsList.append([
            planName,
            telemetryRowIndex, 
            instrumentTelementryRow['microsecondsSinceGpsEpoch'],
            wattsPerM2,
            referenceSpectrumWavelength
            ]
        )

        i += 1

        if i == 30:
            break 

    return resultsList

def plotResults():
    
    print('plotResults ...')
    
    planNames = []
    wavelengths = []
    wattsPerM2s = []

    with open('results.txt') as file:

        next(file) # skip headers

        for result in file:
            planName, telemetryRowIndex, microsecondsSinceGpsEpoch, wattsPerM2, wavelength = result.split(', ')
            wavelengths.append(float(wavelength))
            wattsPerM2s.append(float(wattsPerM2.strip('\n')))

    plt.figure(figsize=(12, 8))
    plt.plot(wavelengths, wattsPerM2s, marker='o', color='b', linestyle='-', linewidth=2, markersize=8)
    plt.xlabel('Wavelength')
    plt.ylabel('WattsPerM2')
    plt.title('WattsPerM2 as a Function of Wavelength (Dark counts not calibrated)')
    plt.grid(True)
    plt.show()

    # file.close()

#Start
print("Executing ...")

plans = pd.read_csv('./data/plans.txt', skipinitialspace=True)

plans['startTimeSeconds'] = plans['startTime'].apply(lambda startTime: startTime / 1000000)
plans['endTimeSeconds'] = plans['endTime'].apply(lambda endTime: endTime / 1000000)
plans['startDateTime'] = plans['startTime'].apply(lambda startTime: datetime(1980, 1, 6, 0, 0, 0, 0) + timedelta(microseconds=startTime))
plans['endDateTime'] = plans['endTime'].apply(lambda endTime: datetime(1980, 1, 6, 0, 0, 0, 0) + timedelta(microseconds=endTime))

headers = [
    'planName', 
    'telemetryRowIndex', 
    'microsecondsSinceGpsEpoch',
    'wattsPerM2',
    'wavelength'
]

with open('results.txt', 'w') as resultsFile:
    resultsFile.write(', '.join(headers) + '\n')

    for index, row in plans.iterrows():

        if row['planName'] == 'Dark':
            continue

        if row['planName'] == 'UpScan': # Format one plan plot nicely before adding UpScan
            continue

        returnedList = processPlan(row['planName'], row['startTime'], row['endTime'])

        for result_row in returnedList:
            formatted_row = ', '.join(map(str, result_row))
            resultsFile.write(formatted_row + '\n')


print('Complete  ...')

plotResults()
