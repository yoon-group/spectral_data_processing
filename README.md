# Automated Spectral Deconvolution for Dye Sample Processing
[![DOI](https://zenodo.org/badge/860567411.svg)](https://doi.org/10.5281/zenodo.17371869)
  
This repository provides a set of Python scripts for processing and analyzing dye-tracer data using a deconvolution-based approach.

---

## File Overview and Execution Order

1. **`Peak_fitter_so.py`**  
   This function file must be located in the same directory as `std_water.py` and `btc_so.py`.  
   It is called internally by these scripts to perform peak fitting during data processing.

2. **`std_water.py`**  
   This script processes all dye standard samples prepared in water.  
   It computes a linear correlation between dye concentration and peak area, which is then used to calibrate water samples for each campaign.

3. **`btc_so.py`**  
   This script generates breakthrough curves (BTCs) for all water samples in every campaign.  
   Results are automatically saved in a folder named **`btc/`**, with subfolders created for each campaign.

---

## Input Data

All dye water sample data (including dye standards) can be downloaded from **HydroShare** at the following link:

🔗 [Spectrofluorophotometer Dye Water Samples](https://www.hydroshare.org/resource/25df1ed10eee4c9da2595a663b87c67b/)

Before running the scripts, ensure that the directory paths in your code correspond to the final locations of these files on your computer.

---

## Output Location

By default, `btc_so.py` will create a folder named **`results/`** in the working directory defined in the script.  
You may modify this directory path to specify a different output location if desired.

---

## Workflow Summary

1. Place all three scripts (`Peak_fitter_so.py`, `std_water.py`, `btc_so.py`) in the same directory.  
2. Run `std_water.py` to calibrate dye concentration relationships.  
3. Run `btc_so.py` to generate BTC curves for all campaigns.  
4. Review the results saved in the automatically generated **`results/`** and **`btc/`** folders.

---

## Requirements

- Python 3.8+  
- Common scientific libraries:  
  `numpy`, `pandas`, `scipy`, `matplotlib`

---

## License

This project is licensed under the [MIT License](LICENSE).  
© 2025 Yoon Research Group, University of Florida


