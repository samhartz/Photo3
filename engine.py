from scipy.optimize import fsolve
from sympy import *
import Tkinter as tk
import numpy as np
import pandas as pd
from math import exp, pi, sqrt
import matplotlib.pyplot as plt
import modelCAMint as model
reload(model)

root = tk.Tk()  #This creates a window, but it won't show up

#Default settings
weatherOption = "ACTUAL"  #options: ACTUAL, AVG, BLM, STEP, CONST 
lightOption = "ACTUAL" #options: ACTUAL, STEP, PAR
rainOption = "DRYDWN" #options: ACTUAL, STOCH, DRYDWN, CONST
duration = tk.IntVar(value = 30) #days (must be an integer)
sInit = tk.DoubleVar(value = 0.5)
soilVariable = tk.StringVar(value="Sandy loam") # default value
capVariable = tk.StringVar(value = "No")
spVariable = tk.StringVar(value = "Triticum aestivum") # default value
locVariable = tk.StringVar(value = "Zacatecas, Mexico") # default value
rainVariable = tk.StringVar(value = "Drydown")

capOptions = {"No": False, "Yes": True}
speciesOptions = {"Triticum aestivum": "T. aest", "Sorghum bicolor": "S. bico", "Opuntia ficus-indica": "O. ficu"}
rainOptions = {"Drydown":"DRYDWN", "Constant": "CONST", "Input":"ACTUAL", "Stochastic": "STOCH"}
pType = {'T. aest': 'C3', 'P. menz': 'C3', 'O. ficu': 'CAM', 'A. tequ': 'CAM', 'S. bico': 'C4'} #Photosynthetic type

#Take inputs from GUI
def makemenu(root, caption, rowno, default, *options):
    lbl = tk.Label(root, text = caption)
    lbl.grid(row = rowno, column = 0)
    menu = tk.OptionMenu(root, default, *options)
    menu.grid(row = rowno, column = 1)

def makeentry(root, caption, rowno, default):
    lbl = tk.Label(root, text = caption)
    lbl.grid(row = rowno, column = 0)
    entry = tk.Entry(root, textvariable = default)
    entry.grid(row = rowno, column = 1)
    #entry.insert(0, default)

def ok():
    global capOn, species, sType, locs, Duration, s0
    Duration = duration.get()
    s0 = sInit.get()
    capOn = capOptions[capVariable.get()]
    species = speciesOptions[spVariable.get()]
    sType = soilVariable.get()
    rainOption = rainOptions[rainVariable.get()]
    root.destroy()
    

speciesmenu =  makemenu(root, 'Plant species:', 0, spVariable,  *speciesOptions.keys())
capmenu = makemenu(root, 'Plant water storage?', 1, capVariable,  *capOptions.keys())
soilmenu = makemenu(root, 'Choose soil type:', 2, soilVariable,  "Sandy loam", "Loamy sand","Loam", "Clay")
#locationsmenu = makeentry(root, 'Location:', 3, locVariable,  "Temple, TX", "Zacatecas, Mexico","Sicily, Italy",
#                     "Pernambuco, Brazil", "Addis Ababa, Ethiopia", "Karnal, India")
rainmenu = makemenu(root, 'Rain:', 3, rainVariable,  *rainOptions.keys())
durationinput = makeentry(root, 'Duration (days):', 4, duration)
moistinput = makeentry(root, 'Init. soil moisture:', 5, sInit)

button = tk.Button(root, text="Run", command=ok)
button.grid(row=6, column = 1)
root.mainloop()                 #This command will tell the window to come out

#Load weather data from excel files
# if locVariable.get() == "Temple, TX":
#     df = pd.read_csv("sample_data\TempleAprilInterp.csv"); 
# elif locVariable.get() == "Zacatecas, Mexico":
#     df = pd.read_csv("sample_data\MexicoDataApril1.csv");
# elif locVariable.get() == "Addis Ababa, Ethiopia":
#     df = pd.read_csv("sample_data\EthiopiaDataApril1.csv");
# elif locVariable.get() == "Karnal, India":
#     df = pd.read_csv("sample_data\IndiaDataApril1.csv"); 
# elif locVariable.get() == "Pernambuco, Brazil":
#     df = pd.read_csv("sample_data\BrazilDataApril1.csv");
# else:
#     df = pd.read_csv("sample_data\ItalyDataApril1.csv");
# df = xl.parse('Average')

df = pd.read_csv("sample_data\TempleAprilInterp.csv");

if weatherOption == 'ACTUAL':
    tempC = df['Temperature']; #extracts temperature column from Excel Worksheet in Celcius
    tempInp = tempC + 273.; #convert to K
    tempInp = list(tempInp.values);
    rh = df['Relative Humidity']; #extracts rh column from Excel Worksheet (%)
    po = 101.325*10**3.; #Atmospheric pressure (Pa)
    psat = 610.78*np.exp(tempC / (tempC + 238.3)*17.2694);  #saturated vapor pressure in Pa
    qaInp = 0.622*rh/100.*psat/po; #needs to be in kg/kg
    qaInp = list(qaInp.values);
else:
    tempInp = 0;
    qaInp = 0;

if lightOption == 'ACTUAL':
    globalsr = df['GHI']; #extracts global solar radiation column from Excel Worksheet in W/m^2
    srInp = list(globalsr.values);
else:
    srInp = 0;

if rainOption == 'ACTUAL':
    rain = df['Rain']; 
    rInp = list(rain.values)
else:
    rInp = 0;

#Run model
data = model.run(species, sType, capOn, weatherOption, lightOption, rainOption, Duration, tempInp, qaInp, srInp, rInp, s0);

#Save data
data.to_pickle('sample_output/CAM_Mmaxfull_Temple_drydown_s05_OptStom')


#Display output graphically
startDay = 0
endDay = 30
displayDays = endDay-startDay;
timevec = np.linspace(startDay,endDay,144*displayDays);

# if pType[species] == "CAM":
#     Mp = plt.figure()
#     plt.xlabel("time (d)")
#     plt.ylabel("M")
#     plt.title("malic acid concentration")
#     plt.plot(timevec, data['M'][144*startDay:144*endDay])
#     Mp.show()
#     zp = plt.figure()
#     plt.xlabel("time (d)")
#     plt.ylabel("z")
#     plt.title("Tonoplast order")
#     plt.plot(timevec, data['z'][144*startDay:144*endDay])
#     zp.show()
# else:
#     pass

# Atot = sum([x*24*3.6/1000 for x in Ana][288:575])/288; #gives total assimilation for one day in mol/m2/day
# Etot = sum([x*24*3.6 for x in Ev][288:575])/288; #gives total evaporation for one day in mm/day
# WUE = Atot/(Etot/18.02) # units of mmol CO2/mol H2O

# #plt.figure()
# # plt.hold(True)
# # plt.plot(z)
# # plt.show()
# # plt.plot(np.array(M)/Mmax[species])
# # plt.show()
# sp = plt.figure()
# plt.title("soil moisture")
# plt.xlabel("time (d)")
# plt.ylabel("s (%)")
# plt.plot(timevec, data['s'])
# sp.show()

# sp = plt.figure()
# plt.title("soil potential")
# plt.xlabel("time (d)")
# plt.ylabel("psi_s (MPa)")
# plt.plot(timevec, data[1])
# sp.show()

# lp = plt.figure()
# plt.title("leaf water potential")
# plt.xlabel("time (d)")
# plt.ylabel("psi_l (MPa)")
# plt.plot(timevec, data[2])
# plt.ylim(-10,1)
# lp.show()


# ap = plt.figure()
# plt.title("Temperature")
# plt.xlabel("time (d)")
# plt.ylabel("T(K)")
# plt.plot(timevec, data['Tl'], label = "Tl")
# plt.hold(True)
# del tempInp[-1]
# del tempInp[-2]
# plt.plot(timevec, tempInp, label = "Ta")
# plt.legend()
# ap.show()

# ap = plt.figure()
# plt.title("atmospheric potential")
# plt.xlabel("time (d)")
# plt.ylabel("psi(MPa)")
# plt.plot(timevec, data['psi_atma'], label = "psi_atm")
# plt.hold(True)
# plt.plot(timevec, data['psi_ll'], label = "psi_l")
# plt.hold(True)
# plt.plot(timevec, data['psi_s'], label = "psi_s")
# plt.legend()
# ap.show()

# psinp = plt.figure()
# plt.title("phi")
# plt.xlabel("time (d)")
# plt.ylabel("phi")
# plt.plot(timevec, srInp[0:144*displayDays])
# psinp.show()

# psinp = plt.figure()
# plt.title("Ta")
# plt.xlabel("time (d)")
# plt.ylabel("Ta")
# plt.plot(timevec, tempInp[0:144*displayDays], label = "Ta")
# plt.hold(True)
# plt.plot(timevec, data['Tl'], label = "Tl")
# plt.legend()
# psinp.show()

# qp = plt.figure()
# plt.title("q")
# plt.xlabel("time (d)")
# plt.ylabel("q(kg/kg)")
# plt.plot(timevec, data['qa'], label = "qa")
# #plt.hold(True)
# #plt.plot(timevec, data['qi'], label ="qi")
# #plt.legend()
# qp.show()

# cp = plt.figure()
# plt.title("[C02]")
# plt.xlabel("time (d)")
# plt.ylabel("C02 (ppm)")
# plt.plot(timevec, data['Cbs'], label = "Cbs")
# plt.hold(True)
# plt.plot(timevec, data['Cm'], label ="Cm")
# # if species =="O. ficu":
# #     plt.hold(True)
# #     plt.plot(timevec, data['Cc'], label ="Cc")
# plt.legend()
# cp.show()

# ep = plt.figure()
# plt.title("Transpiration")
# plt.xlabel("time (d)")
# plt.ylabel("Ev (mm/d)")
# plt.plot(timevec, data['Ev'])
# ep.show()

# # if capOn==True:
# #     wp = plt.figure()
# #     plt.title("stored water")
# #     plt.xlabel("time (d)")
# #     plt.ylabel("w (%)")
# #     plt.plot(timevec, data['w'])
# #     wp.show()
# # else:
# #     pass

# stp = plt.figure()
# plt.title("Stomatal conductance to water vapor")
# plt.xlabel("time (d)")
# plt.ylabel("gsv (mm/s)")
# plt.plot(timevec, data['gsv'])
# stp.show()

# stp = plt.figure()
# plt.title("VPD")
# plt.xlabel("time (d)")
# plt.ylabel("VPD (Pa)")
# plt.plot(timevec, data['VPD'])
# stp.show()

# # pcp = plt.figure()
# # plt.title("plant conductance")
# # plt.xlabel("time (d)")
# # plt.ylabel("gp (mm/s)")
# # plt.plot(timevec, data['gp'])
# # pcp.show()

# # tep = plt.figure()
# # plt.title("leaf temp")
# # plt.xlabel("time (d)")
# # plt.ylabel("Tl (K)")
# # plt.plot(timevec, data['Tl'])
# # tep.show()

anp = plt.figure()
plt.title("Carbon assimilation")
plt.xlabel("time (d)")
plt.ylabel("An (mol/m2/d)")
#plt.plot(timevec, data['Vp'][144*startDay:144*endDay], label = 'Vp')
plt.plot(timevec, data['An'][144*startDay:144*endDay], label = 'An')
plt.legend()
anp.show()

# flp = plt.figure()
# plt.title("Carbon fluxes")
# plt.xlabel("time (d)")
# plt.ylabel("F (mol/m2/d)")
# plt.plot(timevec, data['Asv'][144*startDay:144*endDay], label = 'Asv')
# plt.hold(True)
# plt.plot(timevec, data['Avc'][144*startDay:144*endDay], label = 'Avc')
# plt.hold(True)
# plt.plot(timevec, data['Asc'][144*startDay:144*endDay], label = 'Asc')
# plt.hold(True)
# plt.plot(timevec, data['Rdc'][144*startDay:144*endDay], label = 'Rdc')
# plt.hold(True)
# plt.plot(timevec, data['Rdv'][144*startDay:144*endDay], label = 'Rdv')
# plt.legend()
# flp.show()

# qp = plt.figure()
# plt.title("qa")
# plt.xlabel("time (d)")
# plt.ylabel("(kg/kg)")
# plt.plot(timevec, qaInp)
# qp.show()

# if capOn==True:
#     plt.figure()
#     plt.title("plant water content")
#     plt.xlabel("time (d)")
#     plt.ylabel("w (%)")
#     plt.plot(ta, np.array(Vwa)/Vwt[species])
#     plt.show()
# else:   
#     pass



# plt.figure()
# plt.title("leaf water potential")
# plt.xlabel("time (d)")
# plt.ylabel("psil (MPa)")
# plt.plot(ta[(Duration-displayDays)*24*60/timestepM:Duration*24*60/timestepM], psi_ll[(Duration-displayDays)*24*60/timestepM:Duration*24*60/timestepM])
# plt.show()

# plt.figure()
# plt.title("Total carbon assimilation")
# plt.xlabel("time (d)")
# plt.ylabel("A (mol/m2)")
# xd = [x/1000000*dt for x in Ana]
# xsum = [sum(xd[0:t]) for t in arange(Steps(Duration, timestepM))]


# # plt.plot(ta[(Duration-displayDays)*24*60/timestepM:Duration*24*60/timestepM], xsum[(Duration-displayDays)*24*60/timestepM:Duration*24*60/timestepM])
# # plt.show()
