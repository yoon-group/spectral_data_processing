## Deconvolution approach
This code has been prepared to use in the following order:
1. "Peak_fitter_so.py" function should be in the same directory as the std_water and btc_so. This function will be called on the next routines.
2. The script "std_water.py" processes all dye standards prepared in water. It computes a linear correlation between dye concentrations and peak dye areas, which is then used to calibrate water samples for each campaign.
3. "btc_so.py" will provide the BTC curves, for all the water samples on every campaign. All the results will be storage in a folder, created automatically, and named as "btc". Inside will be another folder with the name of the campaign.

## Calling water samples campaings
All the dye water samples including dye standards can be download by the hydroShare website, data labeled as "Spectrofluorophotometer dye water samples"
1. Ensure that the directory path corresponds to the final location of the files on your computer prior to executing btc_so.py.
2. Difine also the final location of the results, by default, the script will create a folder called "results" and storage in the work directory defined on previous section.
