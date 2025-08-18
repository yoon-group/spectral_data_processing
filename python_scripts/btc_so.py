#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# ========================================================================
#         Dye tracing water samples
#
#                                                      Yoon Research Group
# -------------------------------------------------------------------------
# The data files naming rule: "trip#_sampler#_bottle#_CorrectionData.txt"
# Identify the format campaign name files
# Format file name: "1_1_1_CorrectionData.txt",  ntrip_nsampler_nbottle_ 
# campaign = 'SinkRise2024'  "for example"
#
# Format structure  "MLAC.HZ.220914.1000_CorrectionData.txt" or "MLAC.HZ.220914.2000.10x_CorrectionData.txt"  name.place.date.time.dilutionfactor_  format date (YY,MM,DAY) (HHMM)
# "MLAC.HZ.220918.1000.10x_CorrectionData.txt" or  "MLAC.HZ.220918.1000_CorrectionData.txt" for the same sample, but one is diluted, the diluted will be analized.
# campaign = 'MLAC2022'  "for example"
#
# Format file name:  "MLAC.HZ.230911.1030_CorrectionData.txt"  or  "MLAC.HZ.230912.1930.10x_CorrectionData.txt"  name.place.date.time.dilutionfactor_ 
# campaign = 'MLAC2023' "for example"
#
# Forman file name: "Bear.230502.2330.C12.50x_CorrectionData.txt" or   "Bear.230506.1115.3D15_CorrectionData.txt"  name.date.time.sample.sampler.dilutionfactor_
# campaign = "BearSpring2023"
# =========================================================================
# ============================ Define campaign  =============================

###LIbrarires required
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from matplotlib.ticker import LogLocator, FuncFormatter, LogFormatterSciNotation
from peak_fitter_so import peak_fitter_so
import re
from collections import defaultdict


##Uncomment to run the campaign 
#campaing= "XXXX"

#campaign = "SinkRise2024"
#campaign = "BearSpring2023"
campaign = 'MLAC2023'

# ============================ set up directory =============================
archiveDir = '/home/public/dyeTracingData'
resultsDir = os.path.join(os.getcwd(),"btc", campaign)  # Save data to .... your location 

##Create directories if they don't exist  
os.makedirs(resultsDir, exist_ok=True)
figDir = os.path.join(resultsDir, 'fig')
os.makedirs(figDir, exist_ok=True)
data_dir = os.path.join(archiveDir, campaign, 'sampleData')
sample_type = 'water'  # either {'water', 'eluent'}


# ==================== Defining systematic structure ====================
## Define the next parameters
t_SrBInj = datetime(2024, 11, 18, 14, 10)  # datetime for SrB injection (yyyy,mm,dd,HH,MM)
t_uranineInj = datetime(2024, 11, 18, 14, 32)  # datetime for uranine injection
nTrip = 4  ### num of campaigns
nSampler = 4  ## num of autosamplers
nBottle = 24  ## num of bottles per sampler

#### for the case 1_1_1_2x_
dates = np.empty((nTrip, nSampler, nBottle), dtype=object)  
for iTrip in range(nTrip):
    if iTrip == 0:
        dt_sampler = timedelta(minutes=30)  ##time resolution of every sample
        dt_bottle = timedelta(hours=2)      #delay of every sampler to start
        t0 = datetime(2024, 11, 18, 18, 0)  ##start day of the campaign 1
    elif iTrip == 1:
        dt_sampler = timedelta(minutes=30) ##time resolution of every sample
        dt_bottle = timedelta(hours=2)
        t0 = datetime(2024, 11, 20, 16, 0) ##start day of the campaign 2
    elif iTrip == 2:
        dt_sampler = timedelta(minutes=45)
        dt_bottle = timedelta(hours=3)
        t0 = datetime(2024, 11, 22, 13, 0) ##start day of the campaign 3
    else:#iTrip == 3        ### add more bloques in the case of more samplers and  change "else:" by "elif" 
        dt_sampler = timedelta(minutes=105)
        dt_bottle = timedelta(hours=7)
        t0 = datetime(2024, 11, 25, 13, 0) ##start day of the campaign 4
    # else:#iTrip == 4        ### 
    #     dt_sampler = timedelta(minutes=105)
    #     dt_bottle = timedelta(hours=7)
    #     t0 = datetime(2024, 11, 25, 13, 0)   ### for more campaigns  
    
    for iSampler in range(nSampler):
        for iBottle in range(nBottle):
            dates[iTrip, iSampler, iBottle] = t0 + (iSampler * dt_sampler) + (iBottle * dt_bottle)

# Initialize
sampleDate = []
fl_name_list = []
results = []
file_dict = defaultdict(list)

# Regex patterns ()
pattern_trip = re.compile(r'^(\d+)_(\d+)_(\d+)_CorrectionData\.txt$')  # Case: "1_1_1" Sink_Rise for examaple
pattern_named = re.compile(r'^(.*?\.\d{6}\.\d{4}(?:\.[^.]+)*?)(?:\.(\d+)x[a-zA-Z]*)?_CorrectionData\.txt$')  #Other "BS, mlac"


# ==================== File Loading & Matching with the systematic structure ====================
# ==================== File Loading & Matching with the systematic structure ====================

fl_name_struct = os.listdir(data_dir) ### list of all file names from the data_dir
print (fl_name_struct[:2])  ##chech the structure name of your data

# Step 1: Process files
for fl_name in fl_name_struct:
    match_trip = pattern_trip.match(fl_name)
    match_named = pattern_named.match(fl_name)

    if match_trip:
        try:
            iTrip, iSampler, iBottle = map(int, match_trip.groups())
            sample_date = dates[iTrip - 1, iSampler - 1, iBottle - 1]  ##sample_date, variable just for 1_1_1 format
            fl_name_list.append(fl_name)
            sampleDate.append(sample_date)
        except Exception as e:
            print(f"Error processing trip-format file {fl_name}: {e}")
    elif match_named:
        base, mult = match_named.group(1), match_named.group(2)
        file_dict[base].append((fl_name, int(mult) if mult else None))

# Step 2: Select highest dilution per sample (if there area multiple dilutions from the same sample)
for base, files in file_dict.items():
    if any(mult is not None for _, mult in files):
        fl_name_list.append(max(files, key=lambda x: x[1] or 0)[0])          # Keep highest multiplier
    else:
        fl_name_list.append(files[0][0])

# Step 3: Parse date/time from filenames (for non-trip files)
for fl_name in fl_name_list[len(sampleDate):]:  # Only files from named pattern
    try:
        parts = fl_name.replace('_CorrectionData.txt', '').split('.')
        date_part = next((p for p in parts if re.fullmatch(r'\d{6}', p)), None)
        time_part = next((p for p in parts if re.fullmatch(r'\d{4}', p)), None)
        if not date_part or not time_part:
            raise ValueError("Missing date/time part in filename")
        sample_date = datetime.strptime(f"{date_part} {time_part}", "%y%m%d %H%M")
        sampleDate.append(sample_date)
    except Exception as e:
        print(f"Error parsing date from filename {fl_name}: {e}")

# Step 4: Sort chronologically
sampleDate, fl_name_list = zip(*sorted(zip(sampleDate, fl_name_list)))
sampleDate, fl_name_list = list(sampleDate), list(fl_name_list)

# Step 5: Compute elapsed time in hours
time = pd.to_datetime(sampleDate)
time = ((time - time[0]).total_seconds()) / 3600

# Step 6: Initialize area and results
indNow = list(range(len(sampleDate)))
area = [None] * len(sampleDate)


# ==================== Load Existing Results (if available) ====================
# ==================== Load Existing Results (if available) ====================
pre_results = os.path.join(resultsDir, f'field_samples_{campaign}.npz')
if os.path.exists(pre_results):
    previous_data = np.load(pre_results, allow_pickle=True)     # Loading the previously saved data
    sampleDate_pre = pd.to_datetime(previous_data['date'])  # Ensure datetime type  "date
    area_pre = previous_data['area']
    n_iter_pre = previous_data.get('n_iter', 0) ###assign zero to run for the first time.

if 'area_pre' in locals():
    for i, t in enumerate(sampleDate):
        match_idx = next((j for j, tp in enumerate(sampleDate_pre) if tp == t), None)
        if match_idx is not None and match_idx < len(area_pre):
            area[i] = area_pre[match_idx]

# Extract dilution factors
def extract_dilution(name):
    match = re.search(r"(\d+)x[a-zA-Z]*_", name)
    if match:
        return int(match.group(1))
    elif re.search(r"\bx[a-zA-Z]*_", name):
        return 1
    return 1

#df = [extract_dilution(name) for name in fl_name_list]
df = np.array([extract_dilution(name) for name in fl_name_list])


## ##CORRECTION BY VOLUME,some samples requires dilution by water. Here you can consider those
#df[df == 3] = 4
#df[df == 4] = 5   
df[df == 5] = 6
df[df == 10] = 11
df[df == 20] = 21
# df[df == 25] = 26
# df[df == 50] = 51
# df[df == 100] = 101    


#Assemble results from the last campaign 
n_iter_pre = previous_data.get('n_iter', {}) if 'previous_data' in locals() else {}
results = [{
    "area": area[i],
    "sample": fl_name_list[i],
    "date": sampleDate[i],
    "time": time[i],
    "n_iter": n_iter_pre[i] if i < len(n_iter_pre) else 0,
    "df": df[i]
} for i in range(len(sampleDate))]



# ==================== Process New Samples ====================
# If samples have already been processed in previous runs, they will not be reprocessed
# unless the iteration factor (n_iter) is increased.
# For example, if a campaign was processed and fitted with n_iter=2, those samples
# will only be reprocessed again if n_iter > 2.
# By default, n_iter is set to 15.

if 'results' not in locals():
    results = []

if 'area' not in locals() or len(area) != len(sampleDate):
    area = [None] * len(sampleDate)

# Extract processed sample names from results that have area values
processed_samples = {
    fl_name_list[i]
    for i in range(len(area))
    if area[i] is not None}

# Convert results to a dict for fast lookup by sample name (if not already)
if not isinstance(results, dict):
    results_dict = {res["sample"]: res for res in results if res.get("area") is not None}
else:
    results_dict = results

# Extract processed sample names
processed_samples = {
    name for name, res in results_dict.items()
    if res.get("area") is not None}

to_process = []

# Loop over current indices and process only new samples or ones needing update
#for iSample in [indNow[133]]:  ## uncomment for discrete samples
for iSample in indNow:  ## uncomment to all the samples present in the folder 
    current_sample = fl_name_list[iSample]
    n_iter =15  ### number of fitting times, by default =15
    prev_result = results_dict.get(current_sample)
    
    n_iter_mismatch = (
    prev_result is not None and
    "n_iter" in prev_result and
    prev_result["n_iter"] != n_iter)
    
    if current_sample in processed_samples and not n_iter_mismatch:
        print(f"Skipping, already processed sample: {current_sample}")
        continue
    to_process.append(iSample)  # Only add unprocessed samples
    
# Second loop: Process only samples that need it, with a progress display
for count, iSample in enumerate(to_process, start=1):
    current_sample = fl_name_list[iSample]
    print(f'{count}/{len(to_process)} samples - Processing {current_sample}')
    
    fitted_area = peak_fitter_so(data_dir, figDir, current_sample,n_iter)
    area[iSample] = fitted_area

    # Save updated or new result
    results_dict[current_sample] = {
        "area": fitted_area, ##area for the Uranine,RWT,SrB, currently just Ura, SrB
        "sample": current_sample, ##sample name
        "date": sampleDate[iSample], ##date
        "time": time[iSample], ## hours
        "n_iter": n_iter,  ## by default 15
        "df":df[iSample]  ## dilution factor (e.g., 3x, 10x....4xa... etc) base on the name
    }

# ==================== Save Processed area results, with the meta data of the sample ====================
# Convert dict back to list if needed later
results = list(results_dict.values())
# Safely trim and convert area array
if to_process:
    area = area[:max(to_process)+1]
else:
    area = area[:0]  # or leave area unchanged if already initialized
area = np.array(area)

# Extract columns from results list
area_column = [res["area"] for res in results]
sample_column = [res["sample"] for res in results]
date_column = [res["date"] for res in results]
time_column = [res["time"] for res in results]
n_iter_column = [res.get("n_iter", None) for res in results]  # Handle None for n_iter
df_column=[res["df"] for res in results]

# Save all columns into a single npz file
np.savez(os.path.join(resultsDir, f'field_samples_{campaign}.npz'), 
         area=area_column, 
         sample=sample_column, 
         date=date_column, 
         time=time_column, 
         n_iter=n_iter_column,
         df=df_column)

# ==================== Convert Peak Area to Concentration and ploting ====================
# Here is require the standards for the each dye, (script 2), as results, the'STD_{dyeType}.npy' file is required
# Check the location of that file, it should be in the main directory to be properly called

dyeTypes = ['uranine',"SrB"]  ### Uranine, SrB, RWT
# Loop through dyeTypes
for dyeType in dyeTypes:
    
    ###Call the STD file created with the std_water.py script, if npy file is not found,
    ## copy directly from the std folder and place into the btc folder.
    stdFlName = os.path.join(resultsDir, f'STD_{dyeType}.npy') 
    ###load of slope from standards
    if os.path.exists(stdFlName):
        std_data = np.load(stdFlName, allow_pickle=True).item()
        slope = std_data['slope']
    else:
        print(f"File {stdFlName} not found!")
        continue  # Skip this iteration if the file does not exist

        ##all areas from the peak fitter 
    if dyeType == 'uranine':
        areas = np.array([r['area'][0] for r in results])
    elif dyeType == 'eosin':
            areas = np.array([r['area'][1] for r in results])    
    elif dyeType == 'RWT':
        areas = np.array([r['area'][2] for r in results])
    elif dyeType == 'SrB':
        areas = np.array([r['area'][3] for r in results])
    
    cnc = areas * slope
    
    # Add concentration to results
    for i, r in enumerate(results):
        r[f'cnc_{dyeType}'] = cnc[i]*df[i]  ## concentration applied with the dilution factor
            
    # Extract values from results
    time =np.array([r["time"] for r in results])
    cnc =np.array( [r[f'cnc_{dyeType}'] for r in results])  # Assuming 'area' is the CNC value
    df = np.array([r["df"] for r in results])

    # Plot concentration vs time, normal scale 
    plt.figure(figsize=(14, 6))
    plt.subplot(121)
    plt.plot(time, cnc, 'k-', linewidth=2)
    plt.xlabel('Hours')
    plt.ylabel('Concentration [ppb]')
    plt.xticks(rotation=0)
    plt.grid(True)
    plt.title(f'{dyeType} ({campaign})')

    # Plot concentration vs log(time) scale
    plt.subplot(122)
    plt.plot(time, cnc, 'k-', linewidth=2)
    plt.xscale('log')
    plt.yscale('log')

    # Major ticks at powers of 10
    plt.gca().xaxis.set_major_locator(LogLocator(base=10.0, numticks=10))
    plt.gca().yaxis.set_major_locator(LogLocator(base=10.0, numticks=10))

    # Minor ticks (optional: keep or remove)
    plt.gca().xaxis.set_minor_locator(LogLocator(base=10.0, subs='auto', numticks=10))
    plt.gca().xaxis.set_minor_formatter(plt.NullFormatter())  # Hide minor labels
    
    plt.gca().yaxis.set_minor_locator(LogLocator(base=10.0, subs='auto', numticks=10))
    plt.gca().yaxis.set_minor_formatter(plt.NullFormatter())
    
    # Format major tick labels to avoid scientific notation clutter
    def format_func(value, tick_number):
        if value == 0:
            return "0"
        elif value >= 1:
            return f'{value:.0f}'
        else:
            return f'{value:.3f}'
    
    plt.gca().xaxis.set_major_formatter(FuncFormatter(format_func))
    plt.gca().yaxis.set_major_formatter(LogFormatterSciNotation(base=10.0))  ##scient notation for axis y
    
    # Labels and grid
    plt.xlabel('Travel time [hours]')
    plt.ylabel('Concentration [ppb]')
    plt.grid(True, which='both', linestyle='--', alpha=0.3)
    plt.title(f'{dyeType} log scale')
    plt.xticks(rotation=0)
    plt.tick_params(axis='y', which='major', length=6)

    ##Save the figure in the resultsDir directory
    save_path = os.path.join(resultsDir, f'{dyeType}_concentration_plot.png')
    plt.savefig(save_path, dpi=450)  # Save as PNG with high resolution
    plt.savefig(f'{dyeType}_concentration_plot.png', dpi=450)  # Save as PNG (you can change file type if needed)

    
    #show the plot
    plt.show()
    plt.close()

# ==================== Export Final Results with concentrations calibrated ====================
def flatten_results(results):
    flat_data = []

    for entry in results:
        # Copy and flatten 'area' values
        item = entry.copy()
        area_values = item.pop("area")

        # Dynamically name area columns based on how many values are present
        for i, val in enumerate(area_values, start=1):
            item[f"area{i}"] = val

        flat_data.append(item)

    return pd.DataFrame(flat_data)

processed = flatten_results(results)
save_path = os.path.join(resultsDir, f'results_{campaign}.csv')
processed.to_csv(save_path, index=True)
