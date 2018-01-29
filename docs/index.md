# Photo3
The Photo3 model describes C3, C4, and CAM photosynthesis in a consistent manner. To do this, it uses a model which is built on the Farquhar et al. model for carbon assimilation and incorporates a carbon concentrating mechanism in the case of C4 photosynthesis and a circadian rhythm oscillator and malic acid storage mechanism in the case of CAM photosynthesis. It then calculates stomatal conductance and transpiration based on stomatal optimization theory. The model incorporates soil and atmospheric conditions through a representation of the soil-plant-atmosphere continuum with the option of including plant water storage. Given soil moisture, air temperature, humidity and solar radiation, the model calculates net carbon assimilation and transpiration, among other variables of interest. Photo3 is currently parameterized for three representative species, one from each photosynthetic type: winter wheat (*Triticum aestivum* L.), C3; sorghum (*Sorghum bicolor*), C4; and prickly pear (*Opuntia ficus-indica*), CAM.

# Model Structure

To run the model, simply run the file engine.py. A graphical user interface will display allowing the user to choose a plant species, soil type, soil moisture, duration of the simulation, and a data file containing weather inputs (solar radiation, temperature, and humidity). The engine.py file then calls the file model.py to generate the results, which it then exports and/or plots for the user. 

The sample_data folder contains sample weather inputs for a location in Temple, TX, and the sample_output folder contains results from a few simulations in this area. An academic article describing the model details is currently under review.

# Model requirements
Photo3 was developed in Python 2.7 with the following packages: SciPy, NumPy, Pandas, Tkinter, Matplotlib. We suggest intalling a Python distribution such as [Anaconda][An] to meet these requirements.

[An]: https://www.continuum.io/downloads

# Instructions for formatting weather data for the model input
Half hourly data for solar radiation, temperature, and specific humidity is needed to run the Photo3 model. Hourly data may be obtained from a variety of sources and then interpolated to the model timestep of 30 minutes using the script cleanData.py, which will produce a new excel file which can be read into the model. In order to use cleanData.py, the data must first be formatted with the following columns: Year, Month, Day, Hour, Minute, GHI, Temperature, Relative Humidity. Data cannot contain gaps larger than 23 hours.

The final weather data file supplied to the model should have the headings: Temperature, Relative Humidity, GHI, and should match the model timestep of 30 minutes (see example files in the sample_data folder).

# Recommended sources of data

Solar data, along with temperature and relative humidity data, may be retrieved from the [National Renewable Energy Laboratory database (NSRDB)][nsrdb] (MTS1-3) where it is available for a large number of sites across North America.

Solar data across the globe may be obtained from the [Helicolim-3 database by the Solar Radiation Data Service (SoDa)][soda]. 

[nsrdb]: https://maps.nrel.gov/nsrdb-viewer//?aL=UdPEX9
[soda]: http://soda-pro.com/web-services/radiation/helioclim-3-for-free

Temperature and relative humidity data may be downloaded from the [Iowa Environmental Mesonet (IEM) network][iem], which also includes a link to a Python script for scraping data. Windspeed and precipitation is also available from this source. 

[iem]: http://mesonet.agron.iastate.edu/request/download.phtml

# Authorship

The Photo3 model was created by Samantha Hartzell, Mark Bartlett, and Amilcare Porporato in the Princeton University Department of Civil and Environmental Engineering and the Princeton Environmental Institute.
