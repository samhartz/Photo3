from scipy.optimize import fsolve
from sympy import *
import Tkinter as tk
import tkFileDialog
import numpy as np
import pandas as pd
from math import exp, pi, sqrt, log
import matplotlib.pyplot as plt
from dictionaries import *
from functions import *
#import phenology as phe

root = tk.Tk()  #This creates a window, but it won't show up

heatCap = False
phenologyOn = False

#Default settings
weatherOption = "ACTUAL"  #options: ACTUAL, AVG, BLM, STEP, CONST, VPD
lightOption = "ACTUAL" #options: ACTUAL, STEP, PAR
durationInp = tk.IntVar(value = 90) #days (must be an integer)
sInit = tk.DoubleVar(value = 0.6)
soilVariable = tk.StringVar(value="Loam") # default value
capVariable = tk.StringVar(value = "No")
spVariable = tk.StringVar(value = "Triticum aestivum") # default value
rainVariable = tk.StringVar(value = "Constant")
weatherFile = tk.StringVar(value="sample_data\TempleAprilMay2015.xlsx") #default value for the location of the weather data
repeatWeather = 1 #option to repeat the data in weatherFile x number of times...default value is 1 for no repeat
resultsFile = tk.StringVar(value = 'sample_output\phenologyCAM') #default value for the location where results are saved

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

def ok():
    global capOn, species, sType, locs, duration, s0, rainOption
    duration = durationInp.get()
    s0 = sInit.get()
    capOn = capOptions[capVariable.get()]
    species = speciesOptions[spVariable.get()]
    sType = soilVariable.get()
    rainOption = rainOptions[rainVariable.get()]
    root.destroy()

def chooseinputfile():
    global weatherFile
    root.filename = tkFileDialog.askopenfilename(title = "Select file",filetypes = ( ("Xlsx files","*.xlsx"), ("All files","*")))
    weatherFile.set(root.filename)
    print(weatherFile.get())

def chooseoutputfile():
    global resultsFile
    root.filename = tkFileDialog.asksaveasfilename(title = "Select file",filetypes = (  ("Xlsx files","*.xlsx"), ("All files","*")))
    resultsFile.set(root.filename)
    print(resultsFile.get())

speciesmenu =  makemenu(root, 'Plant species:', 0, spVariable,  *speciesOptions.keys())
capmenu = makemenu(root, 'Plant water storage:', 1, capVariable,  *capOptions.keys())
wLbl = tk.Label(root, text = "Weather data file")
wLbl.grid(row = 2, column = 0)
weatherEntry = tk.Entry(root, textvariable = weatherFile)
weatherEntry.grid(row = 2, column = 1)
filebutton = tk.Button(root, text="Select file", command=chooseinputfile)
filebutton.grid(row=2, column = 2)
soilmenu = makemenu(root, 'Soil type:', 3, soilVariable, "Sand", "Sandy loam", "Loamy sand","Loam", "Clay")
rainmenu = makemenu(root, 'Soil moisture dynamics:', 4, rainVariable,  *rainOptions.keys())
durationinput = makeentry(root, 'Duration (days):', 5, durationInp)
moistinput = makeentry(root, 'Init. soil moisture (%):', 6, sInit)
svLbl = tk.Label(root, text = "Save results as")
svLbl.grid(row = 7, column = 0)
saveEntry = tk.Entry(root, textvariable = resultsFile)
saveEntry.grid(row = 7, column = 1)
savebutton = tk.Button(root, text="Select file", command=chooseoutputfile)
savebutton.grid(row=7, column = 2)
runbutton = tk.Button(root, text="Run", command=ok)
runbutton.grid(row=8, column = 1)
root.mainloop()                 #This command will tell the window to come out

df = pd.read_excel(weatherFile.get())

if weatherOption == 'ACTUAL':
    tempC = df['Temperature'] # extracts temperature column in Celcius
    taInp = tempC + 273. # convert to K

    rh = df['Relative Humidity'] # extracts rh column (%)
    psat = 610.78*np.exp(tempC / (tempC + 238.3)*17.2694)  # saturated vapor pressure in Pa
    qaInp = 0.622*rh/100.*psat/P_ATM # needs to be in kg/kg

    # vpd = df['VPD'] #assumed to be in kPa
    # qaInp = (esat - vpd*1000)*.622/po

    qaInp = list(qaInp.values)
    qaInp = np.tile(qaInp, repeatWeather)
    taInp = list(taInp.values)
    taInp = np.tile(taInp, repeatWeather)
else:
    taInp = 0
    qaInp = 0

if lightOption == 'ACTUAL':
    globalsr = df['GHI'] # extracts global solar radiation column from Excel Worksheet in W/m^2
    phiInp = list(globalsr.values)
    phiInp = np.tile(phiInp, repeatWeather)
else:
    phiInp = 0

if rainOption == 'ACTUAL':
    rain = df['Rain']
    rInp = list(rain.values)
else:
    rInp = 0

#Run model
 #Initial Conditions
#s0 = sInit;  #Initial Soil moisture
vwi =  1.*vwt[species] #Initial Capacitance Depth (m3/m2 leaf area) 

if phenologyOn==True:
    pheDict = phe.initialize(species)
    zr = pheDict['zr'][-1]
    lai = pheDict['lai'][-1]
    ared = pheDict['ared'][-1]
else:
    zr = Zr[species]
    lai = LAI[species]
    ared = 1.

#Declare Arrays for Data Storage
tl_a = [tio] # Leaf Temperature 
shf_a = [] # Sensible Heat Flux
ev_a = [] # Transpiration
s_a = [s0] # Soil Moisture
psi_l_a = [] # Leaf Potential
gsv_a = [] # Stomatal Conductance 
an_a = [] # Assimilated Carbon
ci_a = [] # Intercellular CO2 concentration
cm_a = [] # Mesophyll CO2 concentration
cs_a = [ca] # stomatal carbon concentration
phiL=[]
qi_a=[] # specific humidity internal to leaf stomata
gp_a = [gpmax[species]]

#Extraneous arrays for data storage
psi_atma =[]
vpd_a = [] # Vapor Pressure Deficit

# Arrays for CAM model
if pType[species]=="CAM":
    m_a = [M0]
    z_a = [z0]
    cc_a = [] # Cytoplasm Acid Carbon Concentration 
else:
    pass

# Arrays for C4 model
cbs_a = [ca] # Bundle sheath cell CO2 concentration

#Arrays for capacitance model
vw_a = [vwi] #Array for stored water depth

theta_a = [theta_c(tio, ho)] #Initial potential Temperature
qio = .622*esat(tio)/P_ATM #Initial specific humidity, assuming a starting value at saturation*)
# qaNight = .0082041;#Nighttime specific humidity
# qaDay = .0082041 + .0004293;#(*daytime specific humidity*)
# TaDay = 303.15;
# TaNight = 288.15;
phi_a=[]
if weatherOption == "BLM":
    qa_a = [qio] # Specific humidity
    ta_a = [tio] # Air temperature
    hgt_a = [ho] # Boundary layer height
else: #if weatherOption == "AVG": or "STEP" or "CONST"
    qa_a = []
    ta_a = []

for i in range(steps(duration, int(timestepM))):
    #Update weather conditions
    if lightOption == "PAR":
        phi_a.append(phiPar(i))
    elif lightOption=="STEP": 
        phi_a.append(phiStep(i))     
    elif lightOption=="ACTUAL":
        phi_a.append(phiInp[int((i-1)*60/timestepM)*(timestepM/timestepD)]) #return phiInp[int(t)] #int(t) rounds down to nearest hour integer
    if weatherOption == "AVG":
        ta_a.append(taSimple(i))
        qa_a.append(qaSimple(i))
    elif weatherOption == "STEP":
        ta_a.append(taStep(i))
        qa_a.append(qaStep(i))
        #qa.append(qaSimple(i))
    elif weatherOption == "CONST":
        #Ta.append(TaNight)
        ta_a.append(TaStep(i))
        qa.append(qaNight)
    elif weatherOption == "ACTUAL":
        # Ta.append(taInp[i])
        # qa.append(qaInp[i])
        ta_a.append(taInp[int(i*timestepM/timestepD)])
        qa_a.append(qaInp[int(i*timestepM/timestepD)])
         #Ta.append(taInp[int(i*timestepM/60)])
         #qa.append(qaInp[int(i*timestepM/60)])
    elif weatherOption == "VPD":
        ta_a.append(TaStep(i))
        qa_a.append(qaV(i, VPdiff))
    else: #weatherOption == "BLM"
         pass
    
    vpd_a.append(VPD(ta_a[i], qa_a[i]))
    psi_atma.append(psi_atm(ta_a[i], qa_a[i]))

    if pType[species]=="CAM":
        ci_a.append(ciNew(species, cs_a[i], ta_a[i], qa_a[i]))
        cm_a.append(cmNew(species, cs_a[i], ta_a[i], qa_a[i]))
        cc_a.append(ccNew(species, cs_a[i], ta_a[i], qa_a[i], z_a[i], m_a[i])) 
        CAMargs = {"M": m_a[i], "z": z_a[i]}
    else:
        ci_a.append(ciNew(species, cs_a[i], ta_a[i], qa_a[i]))
        cm_a.append(cmNew(species, cs_a[i], ta_a[i], qa_a[i]))
        CAMargs = {}

    if pType[species]=="C4":
        c1 = cbs_a[i]
    else:
        c1 = cm_a[i]

    #Solving the water balance to find the leaf water potential and leaf temperature
    if len(psi_l_a) == 0:
        psi_l_start = -1.
    else:
        psi_l_start = psi_l_a[-1]

    
    fBalargs = {}

    if capOn == True:
        fBalargs['vwa'] = vw_a[i]
    if pType[species] =="CAM":
        fBalargs['M'] = m_a[i]
        fBalargs['z'] = z_a[i]


    psi_lNew, tlNew = fsolve(fBal, (psi_l_start, 290.), args= (sType, species, phi_a[i], ta_a[i], qa_a[i], tl_a[i], c1, s_a[i], lai, gp_a[i], ared, zr, fBalargs))

    psi_l_a.append(psi_lNew)
    tl_a.append(tlNew)
    qi_a.append(qi(tlNew, psi_lNew))
    gp_a.append(gp(species, psi_lNew))
    ev_a.append(evf(species, phi_a[i], ta_a[i], psi_l_a[i], qa_a[i], tl_a[i], c1, lai, ared, **CAMargs))
    shf_a.append(shfNew(species, tl_a[i], ta_a[i], lai))
    
    #update the soil moisture
    if rainOption == 'CONST':
        s_a.append(s0)
    else:
        if capOn ==True:
            sargs = {'vw':vw_a[i], 'gp':gp_a[i], 'psi_l':psi_l_a[i]}
        else:
            sargs = {}
        if rainOption=="ACTUAL":
            rainAmt = rInp[i]
        elif rainOption == 'STOCH':
            if np.random.random() > lambda_r*timestepM/(60.*24):
                rainAmt =  0.
            else:
                gamma = (n[sType]*zr*100.)/alpha
                rainAmt = np.random.exponential(1./gamma)
        else:
            rainAmt = 0.
        
        s_a.append((min(1., sNew(ev_a[i], s_a[i], zr, **sargs) + rainAmt )))

    if weatherOption == "BLM":
        hgt_a.append(hgtNew(phi_a[i], shf_a[i], hgt_a[i]))
        theta_a.append(theta_aNew(phi_a[i], theta_a[i], shf_a[i], hgt_a[i], hgt_a[i+1]))
        ta_a.append(Temp(theta_a[i+1], ho))
        qa_a.append(qaNew(i))
    else:
        pass
        
    an_a.append(an(species, phi_a[i], psi_l_a[i], tl_a[i], c1, ared, **CAMargs))
    cs_a.append(csNew(species, an_a[i]))
    gsv_a.append(gsvNew(species, phi_a[i], ta_a[i], psi_l_a[i], qa_a[i], tl_a[i], c1, ared, **CAMargs))
    cbs_a.append(cbsNew(an_a[i], cm_a[i], psi_l_a[i]))

    #CAM MODEL
    if pType[species] == "CAM":
        m_a.append(mNew(species, phi_a[i], psi_l_a[i], cc_a[i], tl_a[i], z_a[i], m_a[i]))
        z_a.append(zNew(species, phi_a[i], m_a[i], z_a[i], tl_a[i])) 
    else:
        pass

    #CAPACITANCE MODEL
    if capOn == True:
        vw_a.append(vwf(species, vw_a[i], ev_a[i], gp_a[i], psi_l_a[i], lai))
    else:
        pass

    if phenologyOn ==True:
        pheDict = phe.update(species, an_a[i], tl_a[i], ta_a[i], pheDict)
        zr = pheDict['zr'][-1]
        lai = pheDict['lai'][-1]
        ared = pheDict['ared'][-1]

del tl_a[-1]
del vw_a[-1]
#del hgt_a[-1]
del s_a[-1]
del gp_a[-1]
if pType[species] =="CAM":
    del m_a[-1]
    del z_a[-1]
else:
    pass
if phenologyOn==True:
    del pheDict['lai'][-1]
    del pheDict['zr'][-1]
    del pheDict['bmL'][-1]
    del pheDict['bmS'][-1]
    del pheDict['bmO'][-1]
    del pheDict['bmR'][-1]

if pType[species] =="CAM" and capOn ==True:
    output = np.transpose([gsv_a, [psi_s(sType, x) for x in s_a], psi_l_a, psi_atma, gp_a, s_a, [x*24*3.6 for x in ev_a], tl_a, [x*24*3.6/1000 for x in an_a], qa_a, qi_a, [x/vwt[species] for x in vw_a], vpd_a, ci_a, cm_a, cc_a, [x/Mmax[species] for x in m_a], z_a, pheDict['lai'], pheDict['zr'], pheDict['bmL'], pheDict['bmR'], pheDict['bmS'], pheDict['bmO']])
    output = pd.DataFrame(output, columns = ['gsv', 'psi_s', 'psi_ll', 'psi_atma', 'gp', 's', 'Ev', 'Tl', 'An', 'qa', 'qi', 'w', 'VPD', 'Ci', 'Cm', 'Cc', 'M', 'z', 'lai', 'zr', 'bmL', 'bmR', 'bmS', 'bmO'])
elif pType[species] =="CAM":
    output = np.transpose([gsv_a, [psi_s(sType, x) for x in s_a], psi_l_a, psi_atma, gp_a, s_a, [x*24*3.6 for x in ev_a], tl_a, [x*24*3.6/1000 for x in an_a], qa_a, qi_a, vpd_a, ci_a, cm_a, cc_a, m_a, z_a])
    output = pd.DataFrame(output, columns = ['gsv', 'psi_s', 'psi_ll', 'psi_atma', 'gp', 's', 'Ev', 'Tl', 'An', 'qa', 'qi', 'VPD', 'Ci', 'Cm', 'Cc', 'M', 'z'])
elif capOn==True:
    output = np.transpose([gsv_a, [psi_s(sType, x) for x in s_a], psi_l_a, psi_atma, gp_a, s_a, [x*24*3.6 for x in ev_a], tl_a, [x*24*3.6/1000 for x in an_a], qa_a, qi_a, [x/vwt[species] for x in vw_a], vpd_a, ci_a, cm_a])
    output = pd.DataFrame(output, columns = ['gsv', 'psi_s', 'psi_ll', 'psi_atma', 'gp', 's', 'Ev', 'Tl', 'An', 'qa', 'qi', 'w', 'VPD', 'Ci', 'Cm'])
else: 
    # output = np.transpose([gsv_a, [psi_s(sType, x) for x in s_a], psi_l_a, psi_atma, gp_a, s_a, [x*24*3.6 for x in ev_a], tl_a, [x*24*3.6/1000 for x in an_a], qa_a, qi_a, vpd_a, ci_a, cm_a])
    # output = pd.DataFrame(output, columns = ['gsv', 'psi_s', 'psi_ll', 'psi_atma', 'gp', 's', 'Ev', 'Tl', 'An', 'qa', 'qi', 'VPD', 'Ci', 'Cm'])

    output = np.transpose([gsv_a, [psi_s(sType, x) for x in s_a], psi_l_a, psi_atma, gp_a, s_a, [x*24*3.6 for x in ev_a], tl_a, [x*24*3.6/1000 for x in an_a], qa_a, qi_a, vpd_a, ci_a, cm_a, pheDict['lai'], pheDict['zr'], pheDict['bmL'], pheDict['bmR'], pheDict['bmS'], pheDict['bmO']])
    output = pd.DataFrame(output, columns = ['gsv', 'psi_s', 'psi_ll', 'psi_atma', 'gp', 's', 'Ev', 'Tl', 'An', 'qa', 'qi', 'VPD', 'Ci', 'Cm', 'lai', 'zr', 'bmL', 'bmR', 'bmS', 'bmO'])
    #output = pd.DataFrame({'gsv':gsv_a, 'psi_l' : psi_l_a, 'gp': gp_a, 's':s_a, 'ev': ev_a, 'tl': tl_a, 'an': an_a, 'qa': qa_a, 'qi': qi_a,\
    # 'vpd': vpd_a, 'ci':ci_a, 'cm': cm_a, 'lai': pheDict['lai'], 'zr': pheDict['zr'], 'bmL': pheDict['bmL'], 'bmS': pheDict['bmS'], \
    # 'bmO': pheDict['bmO']}, columns = columns)

data =output
#Save data as DataFrame object
data.to_pickle(resultsFile.get())
#Save data as csv file
#data.to_csv(resultsFile.get())


#Display output graphically
timestepM = 30
startDay = 0
endDay = duration
dispDuration = endDay-startDay
daySteps = 60/timestepM*24
timevec = np.linspace(0,duration,duration*daySteps)
timevecHr = np.linspace(0,duration*24,duration*daySteps)

# anp = plt.figure()
# plt.title("Carbon assimilation")
# plt.xlabel("time (h)")
# plt.ylabel("An (mol/m2/d)")
# #plt.plot(timevec, data['Vp'][144*startDay:144*endDay], label = 'Vp')
# plt.plot(timevecHr[0:144*dispDuration], data['An'][144*startDay:144*endDay], label = 'An')
# #plt.xticks([0.,6.,12.,18.,24.])
# #plt.legend()
# anp.show()

# anp = plt.figure()
# plt.title("Carbon assimilation")
# plt.xlabel("time (d)")
# plt.ylabel("An (mol/m2/d)")
# plt.plot(timevec[0:daySteps*dispDuration], data['An'][daySteps*startDay:daySteps*endDay], label = 'An')
# #plt.xticks([0.,6.,12.,18.,24.])
# #plt.legend()
# anp.show()

# anp = plt.figure()
# plt.title("Biomass")
# plt.xlabel("time (d)")
# plt.ylabel("B (umol/m2 leaf area")
# plt.plot(timevec[0:daySteps*dispDuration], data['B'][daySteps*startDay:daySteps*endDay])
# #plt.xticks([0.,6.,12.,18.,24.])
# #plt.legend()
# anp.show()

# anp = plt.figure()
# plt.title("LAI")
# plt.xlabel("time (d)")
# plt.ylabel("LAI")
# plt.plot(timevec[0:daySteps*dispDuration], data['LAI'][daySteps*startDay:daySteps*endDay])
# #plt.xticks([0.,6.,12.,18.,24.])
# #plt.legend()
# anp.show()

# anp = plt.figure()
# plt.title("Zr")
# plt.xlabel("time (d)")
# plt.ylabel("Zr")
# plt.plot(timevec[0:daySteps*dispDuration], data['Zr'][daySteps*startDay:daySteps*endDay])
# #plt.xticks([0.,6.,12.,18.,24.])
# #plt.legend()
# anp.show()


# ep = plt.figure()
# plt.title("Transpiration")
# plt.xlabel("time (d)")
# plt.ylabel("Ev (mm/d)")
# plt.plot(timevec, data['Ev'])
# ep.show()


# stp = plt.figure()
# plt.title("Stomatal conductance to water vapor")
# plt.xlabel("time (d)")
# plt.ylabel("gsv (mm/s)")
# plt.plot(timevec[0:daySteps*dispDuration], data['gsv'][daySteps*startDay:daySteps*endDay])
# #plt.xticks([0.,6.,12.,18.,24.])
# stp.show()

# stp = plt.figure()
# plt.title("VPD")
# plt.xlabel("time (d)")
# plt.ylabel("VPD (Pa)")
# plt.plot(timevec[0:daySteps*dispDuration], data['VPD'][daySteps*startDay:daySteps*endDay], label = 'VPD')
# plt.hold(True)
# plt.plot(timevec[0:daySteps*dispDuration], [x*10 for x in srInp[daySteps*startDay:daySteps*endDay]] , label = 'phi')
# plt.legend
# stp.show()

# sp = plt.figure()
# plt.title("soil potential")
# plt.xlabel("time (d)")
# plt.ylabel("psi_s (MPa)")
# plt.plot(timevec[0:daySteps*dispDuration], data['psi_s'][daySteps*startDay:daySteps*endDay])
# sp.show()

# sp = plt.figure()
# plt.title("soil moisture")
# plt.xlabel("time (d)")
# plt.ylabel("s (%)")
# plt.plot(timevec[0:daySteps*dispDuration], data['s'][daySteps*startDay:daySteps*endDay])
# sp.show()

# psip = plt.figure()
# plt.title("leaf water potential")
# plt.xlabel("time (d)")
# plt.ylabel("psil (MPa)")
# plt.plot(timevec[0:daySteps*dispDuration], data['psi_ll'][daySteps*startDay:daySteps*endDay])
# psip.show()

# if capOn==True:
#     wp = plt.figure()
#     plt.title("stored water")
#     plt.xlabel("time (d)")
#     plt.ylabel("w (%)")
#     plt.plot(timevec[0:daySteps*dispDuration], data['w'][daySteps*startDay:daySteps*endDay])
#     wp.show()
# else:
#     pass