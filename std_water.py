#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from peak_fitter_so import peak_fitter_so

# ========================================================================
#  
#                    Standard samples processing
#  
#                                                       Yoon Research Group
#  -------------------------------------------------------------------------
#      Dilution factors (sample : Square lake water)
#         20x: 0.2 ml : 4 ml  -> 21 x
#         10x: 0.4 ml : 4 ml  -> 11 x
#          5x: 0.8 ml : 4 ml  ->  6 x
#          2x:   2 ml : 2 ml  ->  2 x
# =========================================================================

# =========================================================================
# Configuration- directory data
# =========================================================================
#campaign = 'SinkRise2024'  
#campaign='MLAC2023'
campaign='BearSpring2023'


archiveDir = '/home/public/dyeTracingData'
archiveDir = os.path.normpath(archiveDir)
resultsDir = os.path.join(os.getcwd(),"std", campaign)  # Save data to '/{campaign}/'
resultsDir_btc = os.path.join(os.getcwd(),"btc", campaign)  # Save data to .... you


#### Create directories if they don't exist
os.makedirs(resultsDir, exist_ok=True)
figDir = os.path.join(resultsDir, 'fig/')
os.makedirs(figDir, exist_ok=True)


# ========================File name recognition================================
dyeTypes = ['SrB', 'uranine']  # 'SrB', 'uranine', 'RWT', 'eosine'  ### This are the names of the dyes to be reading. 

for dyeType in dyeTypes:
    dataDir = os.path.join(archiveDir, campaign, f'STD_{dyeType}/')
    
    # Check if directory exists
    if not os.path.isdir(dataDir):
        print(f"[WARNING] Directory not found: {dataDir}")
        continue
    
    flNameList = [f for f in os.listdir(dataDir) if f.endswith('.txt')]
    # Skip if no .txt files found
    if not flNameList:
        print(f"[WARNING] No .txt files found in: {dataDir}")
        continue   
    nSample = len(flNameList)

    # =========================================================================
    # Standard solutions
    # =========================================================================
    cnc = np.full(nSample, np.nan)  # Concentration
    dilutionFactor = np.ones(nSample)  # Dilution factor

    for iSample in range(nSample):
        flName = flNameList[iSample][:-4]  # Remove '.txt' extension
        note = flName.split('_')

        # Extract concentration
        cnc[iSample] = float(f"{note[1]}.{note[2]}")

        # Extract dilution factor
        if len(note) >= 4:
            if 'x' in note[3]:
                dilutionFactor[iSample] = float(note[3].split('x')[0])
            elif 'X' in note[3]:
                dilutionFactor[iSample] = float(note[3].split('X')[0])

    # Adjust dilution factors
    dilutionFactor[dilutionFactor == 20] = 21
    dilutionFactor[dilutionFactor == 10] = 11
    dilutionFactor[dilutionFactor == 5] = 6

    # =========================================================================
    # Peak Fitting
    # =========================================================================
    area = []
    for iSample in range(nSample):
        nIter = 10 ##10   for standars 1
        if cnc[iSample] > 50: ##50
            nIter = 10 ##50  #for std =1
        elif cnc[iSample] > 1:
            nIter = 50 ##50

        # Call peakFitter function
        area.append(peak_fitter_so(dataDir, figDir, flNameList[iSample], nIter))
        print(f'Sample {iSample+1}/{nSample}') ### samples tracking label

    # Convert area list to a numpy array
    area = np.array(area)

    # Select the appropriate column based on dye type used 
    if dyeType == 'uranine':
        area = area[:, 0]   ###column 0 for uranine
    elif dyeType == 'eosin':
         area = area[:, 1]  ###columm 1 for RWT
    elif dyeType == 'RWT':
        area = area[:, 2]  ###columm 1 for RWT
    elif dyeType == 'SrB':
        area = area[:, 3] ### column 2 for SrB


    # Adjust area by dilution factor
    area = dilutionFactor * area

    # =========================================================================
    # Standard data analysis
    # =========================================================================
    # Define cost function for optimization
    def cost_fn_op(x,area,cnc):
        return np.sqrt(np.sum((area * x - cnc) ** 2))
    
# Optimize to find the slope
    res = minimize(cost_fn_op, 1e-3, args=(area,cnc),options={'maxiter': 1000})  ## optmization using Nelder-Mead method
    slope = res.x[0]
    intercept=res.x[0]
    equation = f'y = {slope:.3e}x + {intercept:.3e}'

    # Calculate SSR and R^2
    SSR = np.sum((area * slope - cnc) ** 2)
    R2 = 1 - SSR / np.sum((cnc - np.mean(cnc)) ** 2)

    # Save results, format npy on STD folder
    flNameSTD = os.path.join(resultsDir, f'STD_{dyeType}.npy')
    np.save(flNameSTD, {'area': area, 'cnc': cnc, 'slope': slope})

    # Save results, format npu on BTC folder, for a direct access of BTC script
    flNameSTD_btc = os.path.join(resultsDir_btc, f'STD_{dyeType}.npy')
    # Ensure the directory exists
    os.makedirs(resultsDir_btc, exist_ok=True)
    np.save(flNameSTD_btc, {'area': area, 'cnc': cnc, 'slope': slope})


    # =========================================================================
    # Plotting
    # =========================================================================
    plt.figure(figsize=(8.5, 11))

    # Linear scale plot (left subplot)
    plt.subplot(121)  # 1 row, 2 columns, 1st subplot
    plt.plot(area, slope * area, '-', label='Fit')  # Plot the fitted line
    plt.scatter(area, cnc, label='Data')  # Plot the data points
    plt.xlabel('Spectral peak area')
    plt.ylabel('Concentration [ppb]')
    plt.title(f'{dyeType} (R^2 = {R2:.4f})')
    plt.legend()
    plt.text(0.05, 0.85, equation, transform=plt.gca().transAxes, fontsize=10, verticalalignment='top')


    # Log scale plot (right subplot)
    plt.subplot(122)  # 1 row, 2 columns, 2nd subplot
    plt.loglog(area, slope * area, '-', label='Fit')  # Plot the fitted line on a log scale
    plt.scatter(area, cnc, label='Data')  # Plot the data points on a log scale
    plt.xlabel('Spectral peak area')
    plt.ylabel('Concentration [ppb]')
    plt.title('Log Scale')
    plt.legend()
    
    # Display the figure
    plt.tight_layout()  # Adjust layout to prevent overlap
   
     # Save figure
    flNameFig = os.path.join(figDir, f'{dyeType}_calibration.png')
    plt.savefig(flNameFig, dpi=300)
    plt.show()
    plt.close()
####
