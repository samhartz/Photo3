# PHOTO-3
PHOTO-3 model describing C3, C4, and CAM photosynthesis

# Model Structure
The bulk of the photo-3 model is contained in the file model.py. The file engine.py takes user inputs, calls the file model.py, and plots the data. 
To run the model, simply run the file engine.py. This file can be edited to change certain variables (such as the duration of the simulation), change the source of the weather inputs, and plot or output specific results of interest. 
The sample_data folder contains sample weather inputs for locations in Waco, TX; Zacatecas, Mexico; Karnal, India; Sicily, Italy, Addis Ababa, Ethiopia, and Pernambuco, Brazil.

# Model requirements
This model was developed in Python 2.7 with the following packages: SciPy, NumPy, Pandas, Tkinter, Matplotlib. We suggest intalling a Python distribution such as [Anaconda][An] to meet these requirements.

[An]: https://www.continuum.io/downloads

# Instructions for obtaining and formatting weather data for the model input
Hourly solar radiation, temperature, and specific humidity is preferred to run the Photo-3 model. This hourly data is then interpolated to the model timestep of 10 minutes using the Matlab code Data_Fill_Interpolate.m, which will produce a new excel file which can be read into the model through the engine.py script. In order to use Data_Fill_Interpolate.m, the data must first be formatted with the following columns: year, month, day, hours, minutes, seconds, relative hum., temp. Data cannot contain gaps larger than 23 hours.

The final weather data file given to the model should be in the format:
Year, Month, Day, Hour, Minute, Temp, Rel Humidity, Global Solar Rad,
and should match the model timestep of 10 minutes (see example files in the sample_data folder).

Recommended sources of data:

Solar radiation, temperature, relative humidity
Solar data, along with temperature and relative humidity data, may be retrieved from the [National Renewable Energy Laboratory database (NSRDB)][nsrdb] (MTS1-3) where it is available for a large number of sites across North America.
Download GHI, temperature, and relative humidity
Solar data across the globe may be obtained from the [Helicolim-3 database by the Solar Radiation Data Service (SoDa)][soda]. 

[nsrdb]: https://maps.nrel.gov/nsrdb-viewer//?aL=UdPEX9
[soda]: http://soda-pro.com/web-services/radiation/helioclim-3-for-free

Temperature, relative humidity, rainfall
Temperature and relative humidity data may be downloaded from the [Iowa Environmental Mesonet (IEM) network][iem], which also includes a link to a Python script for scraping data (**perhaps modify the script to get the data automatically in the right form?). Windspeed and precipitation is also available from this source.**)

[iem]: http://mesonet.agron.iastate.edu/request/download.phtml


