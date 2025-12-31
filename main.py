from scipy.optimize import fsolve
from sympy import *
import numpy as np
import pandas as pd
from math import exp, pi, sqrt, log
import matplotlib.pyplot as plt
from dics import *
from functions import *
from soil import *
from photosynthesis import *
from hydraulics import *
from species_traits import *
from defs import *
import importlib as importlib

# Model settings for duration, input and output file locations, and model timestep
duration = 10 # simulation duration in days
weatherFile = "sample_data\TempleApril2015Interp30.xlsx" # location of weather data with air temperature, air relative humidity, and solar radiation
resultsFile = 'sample_output/test' # default value for the location where results are saved
timestepM = 30 # Model change in time at each step (min)
timestepD = 30 # timestep of input data (min)
dt = timestepM*60. # no. of seconds in timestep, used to advance differential equations (-)

# Import weather data and perform unit conversion
df = pd.read_excel(weatherFile)
tempC = df['Temperature'] # extracts temperature column (C)
taInp = tempC + 273. # convert temperature from C to K
rh = df['Relative Humidity'] # extracts rh column (%)
psat = A_SAT*np.exp((B_SAT*(tempC))/(C_SAT + tempC)) # saturated vapor pressure in Pa
qaInp = 0.622*rh/100.*psat/P_ATM # calculate specific humidity in kg/kg
qaInp = list(qaInp.values)
taInp = list(taInp.values)
phiInp = list(df['GHI'].values)  # extracts global solar radiation column from Excel Worksheet in W/m^2

# Enter simulation setup and initial conditions
sinit = 0.5 # initial volumetric soil moisture content (unitless)
vwi = 0.9 # initial plant volumetric water storage fraction (unitless)
species = Oficu() # plant species (see species_traits.py for species options)
atmosphere = Atmosphere(phiInp[0], taInp[0], qaInp[0]) # atmospheric conditions (solar radiation, temperature, and specific humidity)
soil = Soil(Loam(), DrydownSoil(), species.ZR, sinit) # soil type and moisture content (see soil.py for options)
photo = CAM(species, atmosphere) # plant photosynthesis (see photosynthesis.py for photosynthesis options)
hydro = HydroCap(species, atmosphere, soil, photo, vwi) # plant hydraulics (see hydraulics.py for hydraulics options)
plant = Simulation(species, atmosphere, soil, photo, hydro) # combine all components into a single simulation object

# Perform simulation
for i in range(steps(duration, int(timestepM))):

	plant.update(dt, phiInp[i], taInp[i], qaInp[i])

results = plant.output()

# Export data
data = pd.DataFrame.from_dict(results)
data.to_pickle(resultsFile) # Save data as pandas dataframe
data.to_csv(resultsFile) # Save data as csv file

# Plot results
startDay = 2 # Start time of graph (days)
endDay = duration # End time of graph (days)
dispDuration = endDay-startDay # Duration of graph (days)
daySteps = 60//timestepM*24 # Model timesteps in one day (-)
timevec = np.linspace(0,duration,duration*daySteps) # Time vector for graph (days)
timevecHr = np.linspace(0,duration*24,duration*daySteps) # Time vector for graph (hours)

#Plot soil moisture
anp = plt.figure()
plt.title("Soil moisture")
plt.xlabel("time (d)")
plt.ylabel("s (-)")
plt.plot(timevec[0:daySteps*dispDuration], results['s'][daySteps*startDay:daySteps*endDay])
for i in range(0, dispDuration):
	plt.axvspan(i-0.25, i+0.25, facecolor = 'k', alpha = 0.2)
plt.xlim(0, dispDuration)
#plt.xticks([0.,6.,12.,18.,24.])
#plt.legend()
anp.show()

#Plot leaf water potential
anp = plt.figure()
plt.title("Leaf water potential")
plt.xlabel("time (d)")
plt.ylabel("psi_l (MPa)")
plt.plot(timevec[0:daySteps*dispDuration], results['psi_l'][daySteps*startDay:daySteps*endDay])
for i in range(0, dispDuration):
	plt.axvspan(i-0.25, i+0.25, facecolor = 'k', alpha = 0.2)
plt.xlim(0, dispDuration)
#plt.xticks([0.,6.,12.,18.,24.])
#plt.legend()
anp.show()

#Plot carbon assimilation
anp = plt.figure()
plt.title("Carbon assimilation")
plt.xlabel("time (d)")
plt.ylabel("An (umol/m2/s)")
plt.plot(timevec[0:daySteps*dispDuration], results['a'][daySteps*startDay:daySteps*endDay])
#plt.axvspan(0., 0.25, facecolor = 'k', alpha = 0.3)
for i in range(0, dispDuration):
	plt.axvspan(i-0.25, i+0.25, facecolor = 'k', alpha = 0.2)
plt.xlim(0, dispDuration)
#plt.xticks([0.,6.,12.,18.,24.])
#plt.legend()
anp.show()

# Plot transpiration rate
anp = plt.figure()
plt.title("Transpiration")
plt.xlabel("time (d)")
plt.ylabel("Ev (um/s)")
plt.plot(timevec[0:daySteps*dispDuration], results['ev'][daySteps*startDay:daySteps*endDay])
#plt.axvspan(0., 0.25, facecolor = 'k', alpha = 0.3)
for i in range(0, dispDuration):
	plt.axvspan(i-0.25, i+0.25, facecolor = 'k', alpha = 0.2)
plt.xlim(0, dispDuration)
#plt.xticks([0.,6.,12.,18.,24.])
#plt.legend()
anp.show()


