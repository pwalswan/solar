from datetime import datetime, timedelta
from math import sin, cos
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def getDarkCorrection(planName, planStartTime, planEndTime):
    print('getDarkCorrection ...')

    instrumentTelemetry = pd.read_csv('./data/instrumentTelemetry.txt', skipinitialspace=True)
    detectorTemps = pd.read_csv('./data/detectorTemp.txt', skipinitialspace=True)

    startIndex = instrumentTelemetry[instrumentTelemetry['microsecondsSinceGpsEpoch'] == planStartTime].index[0]
    endIndex = instrumentTelemetry[instrumentTelemetry['microsecondsSinceGpsEpoch'] == planEndTime].index[0]
    
    tempCorrFactor = 0.0061628 # [counts/degC]

    darkCountRates = []

    i = 0 # i is the sample size counter for dark rate
    
    for telemetryRowIndex in range(startIndex, endIndex):
        instrumentTelementryRow = instrumentTelemetry.iloc[telemetryRowIndex]

        # Get the detectorTemp. All the data provided is time-tagged. Using the telemetry index for the detectorTemps file
        detectorTempsRow = detectorTemps.iloc[telemetryRowIndex]
        detectorTemp = detectorTempsRow['temp (C)']

        dcr = instrumentTelementryRow['counts']
        # Convert to seconds
        dark_integration_time = instrumentTelemetry['microsecondsSinceGpsEpoch'] / 1_000_000

        # Write dark_count_rate
        dark_count_rate = dcr / dark_integration_time
        dark_count_rate_corr = dark_count_rate * (1.0 + tempCorrFactor * (20.0 - detectorTemp))

        darkCountRates.append(dark_count_rate_corr)

        i += 1
        
        if i == 10:
            break # Just need a few dark samples to get a typical baseline value

    median_dark_count_rate = np.median(darkCountRates)

    return median_dark_count_rate

def processPlan(planName, planStartTime, planEndTime, median_dark_count_rate_arg):
    print('processPlan ' + planName + ' ...')
    
    # start 2009-11-28 03:17:30.850
    # end   2009-11-28 07:45:27.870

    # At the Downscan start time, the grating position in instrumentTelemetry goes from 0 to non-zero: 107189.77591053583
    # instrumentTelementry line 29159
    # At the Downscan end time, the grating position goes back to zero
    # instrumentTelementry line 31696 (no 31687)
    # The difference between the indexes in the DownScan period is 31969 - 29159 = 2528
    # The number of rows in the results.txt for DownScan is a match: 2528

    # UpScan starts on
    # instrumentTelemetry line 38938
    # UpScan end on 
    # instrumentTelementry eof
    # integrationTime.txt line 23515 and currently set integration time changes from 1000 to 1750
    # There are exact matches for the times when the plans start, and the integration times change.

    # integrationTime interval changes from 1000 to 1750 at the start of the Downscan 
    # The time value is 9.43413450850000 at row 11841
    # Notice on the rows before 11841, the interval is set to 1000, and each time is different
    # by only one number, seconds. You can see the only difference is counding by tens.
    # The currently set integration time changes here.
    # So at this point and after, until it changes, the counts are not necessarily per second,
    # but per integration time.
    # So you need to convert the counts by the currently set integration time
    # The next change in integrationTime is, at the end of Downscan
    # 16263 9.434178728500002E14, 1750.0
    # There is an exact match for the timestamps when a new integration time starts in both integrationTime
    # and instrumentTelemetry.
    # Notice the counts units in intrumentTelemetry.txt are not counts/s
    # The DownScan and UpScan plots overall seem the same (similar). If reduced to 30 rows, the plotting
    # appearance shows, they are plotting the different result text files.

    # Loop through instrumentTelemetry since it has the first factor we need
    instrumentTelemetry = pd.read_csv('./data/instrumentTelemetry.txt', skipinitialspace=True)
    detectorTemps = pd.read_csv('./data/detectorTemp.txt', skipinitialspace=True)
    referenceSpectrum = pd.read_csv('./data/referenceSpectrum.txt', skipinitialspace=True)
    integrationTime = pd.read_csv('./data/integrationTime.txt', skipinitialspace=True)

    # New list and dataframe for results
    resultsList = []

    # Get the start and end index for the plan matched on microseconds
    startIndex = instrumentTelemetry[instrumentTelemetry['microsecondsSinceGpsEpoch'] == planStartTime].index[0]
    try:
        endIndex = instrumentTelemetry[instrumentTelemetry['microsecondsSinceGpsEpoch'] == planEndTime].index[0]
    except:
        endIndex = instrumentTelemetry.index[-1] # Get the last index if the planEndTime exceeds the last intrumentTelemetry reading

    for telemetryRowIndex in range(startIndex, endIndex):
        instrumentTelementryRow = instrumentTelemetry.iloc[telemetryRowIndex]
        gratingPosition = instrumentTelementryRow['gratPos']
        microsecondsSinceGpsEpoch = instrumentTelementryRow['microsecondsSinceGpsEpoch']

        # Wavelength (the grating equation)
        offset = 239532.38
        stepSize = 2.4237772022101214E-6 # [rad/step]
        d = 277.77777777777777 # [nm]
        phiGInRads = 0.08503244115716374 # [rad]
        ang1 = (offset - gratingPosition) * stepSize # [rad]
        wavelength = 2 * d * sin(ang1) * cos(phiGInRads / 2.0) # [nm]  ~178nm

        count_rate = instrumentTelementryRow['counts']

        intTime = 1750  # should be derived from integrationTime.txt intTime (milli-seconds)
        count_rate = count_rate * 1000/intTime # should be derived from integrationTime

        # Get the detectorTemp. All the data provided is time-tagged. Using the telemetry index for the detectorTemps file
        detectorTempsRow = detectorTemps.iloc[telemetryRowIndex]
        detectorTemp = detectorTempsRow['temp (C)']

        # Photon Count Rate (counts/s/m2/nm)
        tempCorrFactor = 0.0061628 # [counts/degC]
        count_rate_corr = count_rate * (1.0 + tempCorrFactor * (20.0 - detectorTemp))

        median_dark_count_rate = median_dark_count_rate_arg
        apertureArea = .01 / (1E2 * 1E2) # [m^2] (aperature area from cm^2 to m^2)
        photonsPerSecondPerM2 = (count_rate_corr - median_dark_count_rate) / apertureArea # [photons/sec/m^2/nm]

        h = 6.62606957E-34 # [J*s]
        c = 299792458.0 # [m/s]
        wavelengthInMeters = wavelength / 10**9
        energyPerPhoton = h * c / wavelengthInMeters # [J]
        wattsPerM2 = photonsPerSecondPerM2 * energyPerPhoton # [watts/m^2/nm]

        referenceSpectrum['irradiance (watts/m^2/nm)'] = pd.to_numeric(referenceSpectrum['irradiance (watts/m^2/nm)'])
        closestIndex = (referenceSpectrum['irradiance (watts/m^2/nm)'] - wattsPerM2).abs().idxmin()

        referenceSpectrumRow = referenceSpectrum.iloc[closestIndex]
        referenceSpectrumWavelength = referenceSpectrumRow['wavelength(nm)']
        
        resultsList.append([
            planName,
            telemetryRowIndex, 
            instrumentTelementryRow['microsecondsSinceGpsEpoch'],
            wattsPerM2,
            referenceSpectrumWavelength
            ]
        )

    resultsList = sorted(resultsList, key=lambda x: x[4]) # sort by wavelength

    return resultsList

def plotIrradianceResult(planName):
    
    print('plotIrradianceResults ...')
    
    planNames = []
    wavelengths = []
    wattsPerM2s = []
    irradianceRatios = []

    with open('results' + planName + '.txt') as file:

        next(file) # skip headers

        for result in file:
            planName, telemetryRowIndex, microsecondsSinceGpsEpoch, wattsPerM2, wavelength = result.split(', ')
            if float(wavelength.rstrip('\n')) > 180 and float(wavelength.rstrip('\n')) < 183:
                wavelengths.append(float(wavelength))
                wattsPerM2s.append(float(wattsPerM2.strip('\n')))
                irradianceRatios.append(float(wattsPerM2.strip('\n')) / float(wavelength))

    plt.figure(figsize=(12, 8))
    plt.plot(wavelengths, wattsPerM2s, marker='o', color='b', linestyle='-', linewidth=2, markersize=8)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Irradiance (WattsPerM2)')
    plt.title('Irradiance as a Function of Wavelength')
    plt.grid(True)
    plt.savefig('irradiance_plot_' + planName + '.png')
    # plt.show()

    plt.plot(irradianceRatios, wavelengths)
    plt.xlabel('Irradiance/Wavelength (WattsPerM2/Wavelengh (nm))')
    plt.ylabel('Wavelength (nm)')
    plt.title('Ratio of Irradiance per Wavelength by Wavelength')
    plt.savefig('irradiance_ratio_' + planName + '.png')

# Start
print("Executing ...")

plans = pd.read_csv('./data/plans.txt', skipinitialspace=True)

# Add additional datetime formats to plans
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

median_dark_count_rate = 0

# Get the dark count correction
for index, row in plans.iterrows():

    if row['planName'] == 'UpScan' or row['planName'] == 'DownScan':
        continue

    if row['planName'] == 'Dark':
        median_dark_count_rate = getDarkCorrection(row['planName'], row['startTime'], row['endTime'])

for index, row in plans.iterrows():

    if row['planName'] == 'Dark':
        continue

    with open('results' + row['planName'] + '.txt', 'w') as resultsFile:
        resultsFile.write(', '.join(headers) + '\n')

        returnedList = processPlan(row['planName'], row['startTime'], row['endTime'], median_dark_count_rate)

        for result_row in returnedList:
            formatted_row = ', '.join(map(str, result_row))
            resultsFile.write(formatted_row + '\n')

    print('Complete ' + row['planName'] + ' ...')

    plotIrradianceResult(row['planName'])
