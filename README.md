# ssp-chart
Sound Speed Profile charts, The Labrador sea, 2019

Data origin:
“These data were collected and made freely available by the International Argo Program and the national programs that contribute to it.  (http://www.argo.ucsd.edu, http://argo.jcommops.org).  The Argo Program is part of the Global Ocean Observing System. Argo (2000). Argo float data and metadata from Global Data Assembly Centre (Argo GDAC). SEANOE. http://doi.org/10.17882/42182”

nc_filter.py - script for filtering profiles; 
settings.txt - filter conditions, coordinate range; 
script nc_filter.py should be put into the folder with input files.


Input: netCDF files. Source files for drifter measurements in The Labrador sea, 2019: ftp://ftp.ifremer.fr/ifremer/argo/geo/atlantic_ocean/2019/; argo_source_files.txt - list of files.

Output: sub folder, named in agreement with search conditions from settings.txt, which contains filtered netCDF files with profiles satisfying to search conditions

data_adj.py - script for data procession and converting to csv. 

Input: Script reads netCDF files from subfolder named nc (just rename output folder obtained by filtering with nc_filter.py); 
filtered_20190101_prof.nc - example of input file. 
Sound speed is calculated using UNESCO equation. Horizons (depths, meter) are calculated from pressure taking into account the latitude/ Linear interpolation on the standard horizons was performed for sound speed, salinity and temperature.

Output: labrador_sea_full.csv, the folliwing data are provided:
sound speed, m/sec; 
temperature, C;
salinity, PSU.

indices (column 0): UTC datatime;
The column 1 is 'name' to indicate kind of data ('SVEL' for sound speed, m/sec, 'TEMP' for temperature, C, 'SALT' for salinity, PSU);
other colunns are the standard horizons (depth) from 5 m to 2000 m;

labrador_sea_full.csv - result of data processing.




