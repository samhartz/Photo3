# Photo3 Plot Compare by Aidan Matthews, Princeton University


import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Reading data
# List files is the list of output files desired to be compared (as .csv files)

list_file = ["sample_output/CAM","sample_output/C4","sample_output/C3"]
# list_file = ["sample_output/C3_highTemp","sample_output/C3_normTemp"]

n = len(list_file)
dict_dfs = {}
for i in range(n):
    dict_dfs[list_file[i]] = pd.read_csv(list_file[i])


# Time steps
duration = 10; # this is the length of the simulation, in days
timestepM = 30; # this is the model timestep, in minutes
daySteps = 60//timestepM*24

timevec = np.linspace(0, duration, duration*daySteps)
timevecHr = np.linspace(0, duration*24, duration*daySteps)

startDay = 0 # start day of display
endDay = duration # end day of display


##% Plots the variable of choice
# 'a' - assimilation (umol/m2/s)
# 'ev' - transpiration (um/sec)
# 'gsv' - stomatal conductance to water (mol/m2/s)
# 'psi_l' - leaf water potential (MPa)
# 'tl' - leaf temperature (K)
# 'gp' - plant conductance (um/s/MPa)

### for model with plant water storage (HydroCap() option)
# 'vw' - stored water volume (m3/m2)
# 'qw' - flux from stored water to connection node x (um/sec)
# 'qs' - flux from soil to connection node x (positive upward) (um/sec)


# Plot 1 variable on different graphs
fig, ax = plt.subplots(n)

xticks = np.linspace(0,duration,6)
for i in range(n):
    ax[i].plot(timevecHr[daySteps*startDay:daySteps*endDay],dict_dfs[list_file[i]]['a'][daySteps*startDay:daySteps*endDay])
    ax[i].set_title(list_file[i][14:])
    #ax[i].set_xlabel("t (h)")
    ax[i].set_ylabel("An (mol/m2/d)")
    ax[i].set_xticks(xticks*24,xticks)
fig.tight_layout()

ax[n-1].set_xlabel("t (days)")



# Plot 1 variable on same graph, 

fig, ax = plt.subplots()

for i in range(n):
    ax.plot(timevecHr[daySteps*startDay:daySteps*endDay],dict_dfs[list_file[i]]['s'][daySteps*startDay:daySteps*endDay],label=list_file[i][14:])
    ax.set_ylabel("s")
fig.tight_layout()
ax.set_xlabel("t (days)")
ax.set_xticks(xticks*24,xticks)
ax.legend()
plt.show()