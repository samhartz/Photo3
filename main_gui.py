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
import gui

timestepM = 30 # Model change in time at each step (min)
timestepD = 30 # timestep of input data 
dt = timestepM*60. # no. of seconds in timestep, used to advance differential equations

duration = gui.duration
weatherFile = gui.weatherLoc
resultsFile = gui.resultsLoc

df = pd.read_excel(weatherFile)
tempC = df['Temperature']
taInp = tempC + 273. # convert to K
rh = df['Relative Humidity'] # extracts rh column (%)
psat = A_SAT*np.exp((B_SAT*(tempC))/(C_SAT + tempC)) # saturated vapor pressure in Pa
qaInp = 0.622*rh/100.*psat/P_ATM # needs to be in kg/kg
qaInp = list(qaInp.values)
taInp = list(taInp.values)
phiInp = list(df['GHI'].values)  # extracts global solar radiation column from Excel Worksheet in W/m^2

#use with GUI
sinit = gui.s0
vwi  = 0.9 # initial water content, as a %
species = gui.species()
atmosphere = Atmosphere(phiInp[0], taInp[0], qaInp[0])
#soil = Soil(gui.sType(), DrydownSoil(), species.ZR, sinit)
soil = Soil(gui.sType(), gui.soilDynamics, species.ZR, sinit)
photo = species.PTYPE(species, atmosphere)
hydro = gui.hydrology(species, atmosphere, soil, photo, vwi)

plant = Simulation(species, atmosphere, soil, photo, hydro)


for i in range(steps(duration, int(timestepM))):

	plant.update(dt, phiInp[i], taInp[i], qaInp[i])

results = plant.output()

data = pd.DataFrame.from_dict(results)
data.to_pickle(resultsFile)
# Save data as csv file
data.to_csv(resultsFile)

#Plot results
startDay = 2
endDay = duration
dispDuration = endDay-startDay
daySteps = 60//timestepM*24
timevec = np.linspace(0,duration,duration*daySteps)
timevecHr = np.linspace(0,duration*24,duration*daySteps)

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
