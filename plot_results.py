from scipy.optimize import fsolve
from sympy import *
import numpy as np
import pandas as pd
from math import exp, pi, sqrt
import matplotlib as matplotlib
import matplotlib.pyplot as plt


C3 = pd.read_pickle('sample_output/Chile_combined2012')
#C4 = pd.read_pickle('sample_output/C4_Temple_drydown_s05_OptStom')
#CAM = pd.read_pickle('sample_output/CAM_Temple_drydown_s05_OptStom')

LAI = {'C3': 5., 'CAM': 3.,'C4': 5.} #Leaf area index (m2/m2)


duration = 5; # this is the length of the simulation, in days
timestepM = 30; # this is the model timestep, in minutes
daySteps = 60//timestepM*24

timevec = np.linspace(0, duration, duration*daySteps)
timevecHr = np.linspace(0, duration*24, duration*daySteps)

startDay = 0 # start day of display
endDay = duration # end day of display

#plot results

# the vector 'a' is originally in umol/m2/s, this converts it to mol/m2/d (per unit leaf area)

anp = plt.figure()
plt.title("Carbon assimilation")
plt.xlabel("t (h)")
plt.ylabel("An (mol/m2/d)")

plt.plot(timevecHr[daySteps*startDay:daySteps*endDay], 3.6*24/1000*C3['a'][daySteps*startDay:daySteps*endDay], 'k', label = 'C3')
# plt.hold(True)
# plt.plot(timevecHr[0:145], C4['a'][0:145], 'k--', label = 'C4')
# plt.hold(True)
# plt.plot(timevecHr[0:145], CAM['a'][0:145], 'k-.', label = 'CAM')
# plt.legend()
anp.show()

# the vector 'ev' is originally in um/sec, this converts it to mm/d (per unit ground area)

evp = plt.figure()
plt.title("Transpiration")
plt.xlabel("t (h)")
plt.ylabel("E (mm/d)")

plt.plot(timevecHr[daySteps*startDay:daySteps*endDay], 3.6*24*C3['ev'][daySteps*startDay:daySteps*endDay], 'k', label = 'C3')
# plt.hold(True)
# plt.plot(timevecHr[0:145], C4['ev'][0:145], 'k--', label = 'C4')
# plt.hold(True)
# plt.plot(timevecHr[0:145], CAM['ev'][0:145], 'k-.', label = 'CAM')
# plt.legend()
evp.show()

# gsv is already in mol/m2/s. this is the stomatal conductance to water.

gsp = plt.figure()
plt.title("Stomatal conductance")
plt.xlabel("t (h)")
plt.ylabel("gs (mol/m2/s)")

plt.plot(timevecHr[daySteps*startDay:daySteps*endDay], C3['gsv'][daySteps*startDay:daySteps*endDay], 'k', label = 'C3')
# plt.hold(True)
# plt.plot(timevecHr[0:145], C4['gsv'][0:145], 'k--', label = 'C4')
# plt.hold(True)
# plt.plot(timevecHr[0:145], CAM['gsv'][0:145], 'k-.', label = 'CAM')
# plt.legend()
gsp.show()

