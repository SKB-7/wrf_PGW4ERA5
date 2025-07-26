# PGW-Simulation for CMOR input data

This is an attempt at a very practical explaination of how to set up a PGW simulation using global climate model data in the CMOR-Format as input (for example CMIP5 or CMIP6 data).

**What data to get?**

You will need data for the following variables: hur, ta, ua, va, zg, pa, hurs, tas, ts, tos

**What time resolution should one choose?**

Monthly mean data is the easiest. This is called e.g. Amon in CMOR. tos is part of Omon in CMOR.

**How to preprocess the data?**

For all variables we need to know how they will change under climate change. This needs to be expressed as a mean annual cycle of changes.
In practice we can get a time slice of the "historical" period (HIST) and from a future period under a certain emission scenario (SCEN) such as "rcp85". A typical example: For the historical period, get data from 1971-2000. Then construct the mean annual cycle for 1971-2000, for example using the [cdo-command "ymonmean"](https://code.mpimet.mpg.de/projects/cdo/embedded/index.html#x1-5370002.8.33). Repeat for 2070-2099 and the rcp85 data. 
Lastly, subtract the historical monthly-mean annual cycle from the future monthly-mean annual cycle. Save the result from the subtraction, or the difference between the two periods, as a single netcdf-file per variable (e.g. delta_ta.nc, delta_hurs.nc, ....).
These netcdf files are needed as input for [setp_02_preproc_deltas.py](/setp_02_preproc_deltas.py), and the naming convention can be specified in [settings.py](/settings.py) (look for the dictionary "file_name_bases").
