# Oyster_salinity_var
Modeling oyster predator-prey dynamics under increasing salinity variability

Matlab code used in the Commander et al. manuscript describing oyster predator-prey dynamics under increasing salinity variability in Apalachicola Bay

Files:

# Setting up the different salinity scenarios:
TimeSeries_Resampling.m Analyzes salinity data from the NERR database to create a climatology and associated data needed for the resampling analysis
Salinity_runs.m is a wrapper function that loops over scenario options to generate the simulated timeseries
Create_Mock_Salinity.m uses the Denny et al. method to generate new simulated timeseries based on the Apalachicola climatology

# Running the oyster model
Oyster_run.m is a wrapper function that calls all of the needed functions for a given simulation scenario
Oyster_IPM_Disturbance.m performs one simulation of the model, given the inputs called by Oyster_run; it calls a number of functions listed below in #model files
Oyster_IPM_Disturbance_Analysis.m creates the figures in the manuscript

# model files
Oyster_Params.m creates a structure file with all parameters used in the model
mkkern.m creates IPM kernels of different types by calling kernmatSimp.m
kernmatSimp.m creates IPM kernels, given options
makeSimpVec.m creates a vector that will perform Simpson's integration 
predator_salinity_penalty.m calculates the appropriate predator attack rate, given a salinity and temperature


