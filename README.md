# Photo3
Photo3 model describing C3, C4, and CAM photosynthesis

# Model Structure
The bulk of the Photo3 model is contained in the file model.py. The file engine.py takes user inputs, calls the file model.py, and plots the data. 
To run the model, simply run the file engine.py. A graphical user interface will display allowing the user to choose a plant species, soil type, soil moisture, duration of the simulation, and a data file containing weather inputs (solar radiation, temperature, and humidity). The engine.py file itself may also be edited to change these variables or output specific results of interest. 
The sample_data folder contains sample weather inputs for a location in Temple, TX, and the sample_output folder contains results from a few simulations in this area.

# Model requirements
This model was developed in Python 2.7 with the following packages: SciPy, NumPy, Pandas, Tkinter, Matplotlib. We suggest intalling a Python distribution such as [Anaconda][An] to meet these requirements.

[An]: https://www.continuum.io/downloads

# Instructions for obtaining and formatting weather data for the model input
Hourly solar radiation, temperature, and specific humidity is preferred to run the Photo3 model. This hourly data is then interpolated to the model timestep of 10 minutes using the Matlab code Data_Fill_Interpolate.m, which will produce a new excel file which can be read into the model through the engine.py script. In order to use Data_Fill_Interpolate.m, the data must first be formatted with the following columns: year, month, day, hours, minutes, seconds, relative hum., temp. Data cannot contain gaps larger than 23 hours.

The final weather data file given to the model should be in the format:
Year, Month, Day, Hour, Minute, Temp, Rel Humidity, Global Solar Rad,
and should match the model timestep of 10 minutes (see example files in the sample_data folder).

Recommended sources of data:

Solar radiation, temperature, relative humidity

Solar data, along with temperature and relative humidity data, may be retrieved from the [National Renewable Energy Laboratory database (NSRDB)][nsrdb] (MTS1-3) where it is available for a large number of sites across North America.


Solar data across the globe may be obtained from the [Helicolim-3 database by the Solar Radiation Data Service (SoDa)][soda]. 

[nsrdb]: https://maps.nrel.gov/nsrdb-viewer//?aL=UdPEX9
[soda]: http://soda-pro.com/web-services/radiation/helioclim-3-for-free

Temperature, relative humidity, rainfall

Temperature and relative humidity data may be downloaded from the [Iowa Environmental Mesonet (IEM) network][iem], which also includes a link to a Python script for scraping data. Windspeed and precipitation is also available from this source. 

[iem]: http://mesonet.agron.iastate.edu/request/download.phtml


