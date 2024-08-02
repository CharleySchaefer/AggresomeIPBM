RAW DATA:
.xlsx - Raw data provided by Yingying
.csv  - Same data as in the .xlsx; but exported to .csv format

DATA PROCESSING
AnalyseSceenTape.m - data processing script - 
  PURPOSE:
    - Import raw intensity data
    - Use control data to define a baseline
    - Filter intensity data using the baseline
    - Calculate ratio between intensities + error propagation
    - Fit model to the intensity ratio
    - export processed_data.txt
  USAGE:
    Can be run in octave or matlab.
    To run in matlab, line 16 should be modified.
    Expected output is shown in the screenshots.

Processed_data.txt
  - Generated from .csv using AnalyseSceenTape.m.
  - columns: RNA length, intensity ratio, error bars
  
  
  
  cutoff 27
  
  R2=0.992474
Prefac   =0.238185 +/- 0.062176
DeltaH/kT=0.005962 +/- 0.000649
R2=0.992474
  
  
