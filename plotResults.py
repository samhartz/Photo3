from scipy.optimize import fsolve
from sympy import *
import Tkinter as tk
import numpy as np
import pandas as pd
from math import exp, pi, sqrt
import matplotlib.pyplot as plt

C3 = pd.read_pickle('sample_output/C3_Temple_drydown_s05_OptStom')
C4 = pd.read_pickle('sample_output/C4_Temple_drydown_s05_OptStom')
CAM = pd.read_pickle('sample_output/CAM_Temple_drydown_s05_OptStom')

LAI = {'C3': 5., 'CAM': 3.,'C4': 5.} #Leaf area index (m2/m2)


Duration = 30;
timestepM=10;
timevec = np.linspace(0,Duration,len(C3['s']));
timevecHr = np.linspace(0,Duration*24,len(C3['s']));

def Steps(Duration, TimeStep):
    """Change Duration of Simulation to to number of timesteps according to timestep value"""
    return(Duration*24*60)/TimeStep
sDay = 24.*60./timestepM


def AdailyTot(pT):
	"""mol/m2, per unit leaf area"""
	return [sum([x/sDay for x in pT['An']][t*Steps(1, timestepM):(t+1)*Steps(1, timestepM)]) for t in np.arange(Duration)]

def EdailyTot(pT):
	"""mm/m2, per unit ground area"""
	return [sum([x/sDay for x in pT['Ev']][t*Steps(1, timestepM):(t+1)*Steps(1, timestepM)]) for t in np.arange(Duration)]

def WUEdailyTot(pT, lai):
	return  np.divide(AdailyTot(pT), [x/18.02/lai for x in EdailyTot(pT)])

timevecDay =np.arange(Duration)

# #Plot long-term trends
# pW = plt.figure()
# plt.title("Daily water use efficiency")
# plt.xlabel("time (d)")
# plt.ylabel("WUE (mmol/mol)")
# plt.plot(timevecDay,WUEdailyTot(C3, 5.), 'k', label = 'C3')
# plt.hold(True)
# plt.plot(timevecDay,WUEdailyTot(C4, 5.), 'k--', label = 'C4')
# plt.hold(True)
# plt.plot(timevecDay,WUEdailyTot(CAM, 3.), 'k-.', label = 'CAM')
# plt.legend()
# pW.show()

# #plot cumulative totals
# antot = plt.figure()
# plt.title("Cumulative carbon assimilation")
# plt.xlabel("time (d)")
# plt.ylabel("A (mol/m2)")
# plt.plot(timevec, [sum([x/sDay for x in C3['An']][0:t]) for t in np.arange(Steps(Duration, timestepM))], 'k', label = 'C3')
# plt.hold(True)
# plt.plot(timevec, [sum([x/sDay for x in C4['An']][0:t]) for t in np.arange(Steps(Duration, timestepM))], 'k--', label = 'C4')
# plt.hold(True)
# plt.plot(timevec, [sum([x/sDay for x in CAM['An']][0:t]) for t in np.arange(Steps(Duration, timestepM))], 'k-.', label = 'CAM')
# plt.legend()
# antot.show()

# evtot = plt.figure()
# plt.title("Cumulative water use")
# plt.xlabel("time (d)")
# plt.ylabel("Ev (mm/m2)")
# plt.plot(timevec, [sum([x/sDay for x in C3['Ev']][0:t]) for t in np.arange(Steps(Duration, timestepM))], 'k', label = 'C3')
# plt.hold(True)
# plt.plot(timevec, [sum([x/sDay for x in C4['Ev']][0:t]) for t in np.arange(Steps(Duration, timestepM))], 'k--',label = 'C4')
# plt.hold(True)
# plt.plot(timevec, [sum([x/sDay for x in CAM['Ev']][0:t]) for t in np.arange(Steps(Duration, timestepM))], 'k-.', label = 'CAM')
# plt.legend()
# evtot.show()

#plot 1 representative day
anp = plt.figure()
plt.title("Carbon assimilation")
plt.xlabel("t (h)")
plt.ylabel("An (mol/m2/d)")
plt.xticks([0.,6.,12.,18.,24.])
plt.plot(timevecHr[0:145], C3['An'][0:145], 'k', label = 'C3')
plt.hold(True)
plt.plot(timevecHr[0:145], C4['An'][0:145], 'k--', label = 'C4')
plt.hold(True)
plt.plot(timevecHr[0:145], CAM['An'][0:145], 'k-.', label = 'CAM')
plt.legend()
anp.show()

evp = plt.figure()
plt.title("Transpiration")
plt.xlabel("t (h)")
plt.ylabel("E (mm/d)")
plt.xticks([0,6,12,18,24])
plt.plot(timevecHr[0:145], C3['Ev'][0:145], 'k', label = 'C3')
plt.hold(True)
plt.plot(timevecHr[0:145], C4['Ev'][0:145], 'k--', label = 'C4')
plt.hold(True)
plt.plot(timevecHr[0:145], CAM['Ev'][0:145], 'k-.', label = 'CAM')
plt.legend()
evp.show()

gsp = plt.figure()
plt.title("Stomatal conductance")
plt.xlabel("t (h)")
plt.ylabel("gs (mm/s)")
plt.xticks([0,6,12,18,24])
plt.plot(timevecHr[0:145], C3['gsv'][0:145], 'k', label = 'C3')
plt.hold(True)
plt.plot(timevecHr[0:145], C4['gsv'][0:145], 'k--', label = 'C4')
plt.hold(True)
plt.plot(timevecHr[0:145], CAM['gsv'][0:145], 'k-.', label = 'CAM')
plt.legend()
gsp.show()
