##Decovolution approach
This code has been prepared to use in the following order:
1. "Peak_fitter.py" is function should be in the same directory as the std_water and btc_so. This function will be called on the next routines.
2. "std_water.py" routine will proccess all the standards prepared in water. This script produce the equation to calibrate water samples for every campaign.
3. "btc_so.py" will provide the BTC curves, for all the water samples on every campaign. All the results will be storage in a folder, created automatically, and named as "btc". Inside will be another folder with the name of the campaign.
