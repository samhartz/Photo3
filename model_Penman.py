from scipy.optimize import fsolve
from sympy import *
import Tkinter as tk
import numpy as np
import pandas as pd
from math import exp, pi, sqrt
import matplotlib.pyplot as plt

heatCap = False
cpw = 4184. #specific heat of water (J/kg/K)
transpirationType = "qi" #qi #PM

#Soil parameters
psi_ss = {'Sandy loam': -.7*10**-3, 'Loamy sand': -.17*10.**-3.,'Loam': -1.43*10.**-3., 'Clay': -1.82*10**-3.} #(MPa)
b = {'Sandy loam': 4.9, 'Loamy sand':4.38, 'Loam': 5.39, 'Clay': 11.4}#(unitless)
Ks = {'Sandy loam': 80.,'Loamy sand':100., 'Loam': 20., 'Clay': 1.} #Saturated hydraulic condutivity (cm/day)
n = {'Sandy loam': .43, 'Loamy sand':.42,'Loam': .45, 'Clay': .5} #Porosity
sh = {'Sandy loam':.14, 'Loamy sand':.08,'Loam': .19, 'Clay': 0.47} #Hygroscopic point

def run(species, sType, capOn, weatherOption, lightOption, rainOption, Duration, taInp, qaInp, phiInp, rInp, sInit):

    #using same Zr for all 3 crops
    Zr=0.5
    #Zr = {'T. aest': 0.75, 'P. menz': 0.65, 'O. ficu': 0.1, 'A. tequ': 0.3, 'S. bico': 0.5} # Rooting depth (m)

    #Plant photosynthetic parameters

    pType = {'T. aest': 'C3', 'P. menz': 'C3', 'O. ficu': 'CAM', 'A. tequ': 'CAM', 'S. bico': 'C4'} #Photosynthetic type
    LAI = {'T. aest': 5., 'P. menz': 8.4, 'O. ficu': 3., 'A. tequ': 6., 'S. bico': 5.} #Leaf area index (m2/m2)

    #Vcmax0 = {'T. aest': 107.4,   'P. menz': 57.7,   'O. ficu': 13., 'A. tequ': 19.5, 'S. bico': 39.} #Max. carboxylation rate (\[Mu]mol/(m^2s^2))
    Vcmax0 = {'T. aest': 83.,   'P. menz': 57.7,   'O. ficu': 13., 'A. tequ': 19.5, 'S. bico': 39.} #Max. carboxylation rate (\[Mu]mol/(m^2s^2))
    Hkc =    {'T. aest': 59430.,  'P. menz': 59430., 'O. ficu': 54930., 'A. tequ': 59430., 'S. bico': 59430.} #Activation Energy for Kc (J/mol)
    Hko =    {'T. aest': 36000.,  'P. menz': 36000., 'O. ficu': 36000., 'A. tequ': 36000., 'S. bico': 36000.} #Activation Energy for Ko (J/mol)
    HkR =    {'T. aest': 53000.,  'P. menz': 53000., 'O. ficu': 53000., 'A. tequ': 53000., 'S. bico': 53000.}  #Activation Energy for Rd (J/mol)
    Rd =     {'T. aest': 4.93,    'P. menz': .32,    'O. ficu': .32, 'A. tequ': .32, 'S. bico': .32} #Standard Dark respiration at 25 C (\[Mu]mol/(m^2s))
    HaV =    {'T. aest': 62000.,  'P. menz': 61210., 'O. ficu': 72000., 'A. tequ': 72000., 'S. bico': 72000.} #Activation Energy for Vc,max (J/mol)
    HdV =    {'T. aest': 202900., 'P. menz': 200000.,'O. ficu': 200000., 'A. tequ': 200000., 'S. bico': 200000.} #Deactivation Energy for Vc,max (J/mol)
    HaJ =    {'T. aest': 50000.,  'P. menz': 37870., 'O. ficu': 50000., 'A. tequ': 50000., 'S. bico': 50000.} #Activation Energy for Jmax (J/mol)
    HdJ =    {'T. aest': 200000., 'P. menz': 200000.,'O. ficu': 200000., 'A. tequ': 200000., 'S. bico': 200000.} #Deactivation Energy for Jmax (J/mol)
    #Jmax0 =  {'T. aest': 184.9,   'P. menz': 98.5,   'O. ficu': 2.*13., 'A. tequ': 2.*19.5, 'S. bico': 180.} #Max. e- transport rate (\[Mu]mol/(m^2s))
    Jmax0 =  {'T. aest': 132.,   'P. menz': 98.5,   'O. ficu': 2.*13., 'A. tequ': 2.*19.5, 'S. bico': 180.} #Max. e- transport rate (\[Mu]mol/(m^2s))

    Dxc02 = {'T. aest': .0021,  'P. menz': .0021,  'O. ficu': .0021,  'A. tequ': .0021,  'S. bico': .0021} #VPD sensitivity (kg/kg)
    a1 =    {'T. aest': .15625, 'P. menz': .15625, 'O. ficu': .15625, 'A. tequ': .15625, 'S. bico': .15625} #VPD sensitivity

    Dxnew = {'T. aest': 300.,  'P. menz': 300.,  'O. ficu': 300.,  'A. tequ': 300.,  'S. bico': 300.}
    a1new = {'T. aest': 15., 'P. menz': 15., 'O. ficu': 0.6*15., 'A. tequ': 0.6*15., 'S. bico': 0.5*15.}

    Kc0 = {'T. aest': 302., 'P. menz': 302., 'O. ficu': 302., 'A. tequ': 302., 'S. bico': 302.} #Michaelis constant for C02 at To (\[Mu]mol/mol)
    Ko0 = {'T. aest': 256., 'P. menz': 256., 'O. ficu': 256., 'A. tequ': 256., 'S. bico': 256.} #Michaelis constant for 02 at To (mmol/mol)
    oi =  {'T. aest': .209, 'P. menz': .209, 'O. ficu': .209, 'A. tequ': .209, 'S. bico': .209} #Oxygen Concentration (mol/mol )
    SvC = {'T. aest': 649., 'P. menz': 649., 'O. ficu': 649., 'A. tequ': 649., 'S. bico': 649.} #Entropy term for carboxylation (J/mol)
    SvQ = {'T. aest': 646., 'P. menz': 646., 'O. ficu': 646., 'A. tequ': 646., 'S. bico': 646.} #Entropy term for e-transport (J/mol)

    To = 293.2; #Reference Temperature for photosynthetic parameters (K)
    kappa_2 = .3; #Quantum yield of photosynthesis (mol CO2/mol photon)
    gamma_o = {'T. aest': 34.6, 'P. menz': 34.6, 'O. ficu': 34.6, 'A. tequ': 34.6, 'S. bico': 10.}  #C02 compensation point at To (umol/mol)
    gamma_1 = .0451; #(1/K)
    gamma_2 = .000347; #(1/K^2)

    #Plant hydraulic parameters
    gcut = {'T. aest': 0.15*2., 'P. menz': 0.007, 'O. ficu': 0., 'A. tequ': 0., 'S. bico': 0.1802} #Cuticular conducatnce, per unit leaf area (mm/s)
    ga = {'T. aest': 61., 'P. menz': 324., 'O. ficu': 324., 'A. tequ': 61., 'S. bico': 61.} #Atomospheric  Conductance, per unit ground area (mm/s)
    gmgsratio = {'C3': 1.65, 'C4': 2.65, 'CAM': 1.}; #gm/gs (mol/(m^2-s)) 
    Rc = {'C3': 0.7, 'C4': 0.4, 'CAM': 0.5}; #Ratio of ci:ca
    RAIW = {'T. aest': 5.6, 'P. menz': 10., 'O. ficu': 3., 'A. tequ': 3., 'S. bico': 5.6} #Root area index (m2/m2)
    froot = {'T. aest': .25, 'P. menz': .25, 'O. ficu': .25, 'A. tequ': .25, 'S. bico': .25} #Ratio of root to total biomass
    a = 8.; #Parameter accounting for root growth

    psi_la1O = {'T. aest': -.7, 'P. menz': -.5, 'O. ficu': -.5, 'A. tequ': -.5, 'S. bico': -.5} 
    psi_laoO = {'T. aest': -2., 'P. menz': -3., 'O. ficu': -3., 'A. tequ': -3., 'S. bico': -1.8} 

    gwmax = {'T. aest': 0., 'P. menz': 0.05, 'O. ficu': .002, 'A. tequ': .002, 'S. bico': 0.} #Conductance b/t psi_x and psi_w, per unit root area (um/(MPa s))
    gpmax = {'T. aest': 11.7, 'P. menz': 0.056, 'O. ficu': .04, 'A. tequ': .04, 'S. bico': 0.13} #Plant conductance, per unit leaf area(um/(MPa s))
    Vwt = {'T. aest': .000001, 'P. menz': .027/LAI['P. menz'], 'O. ficu': .00415, 'A. tequ': .00415, 'S. bico': .000001} #Max. storage depth (m3/m2 leaf area)
    cap = {'T. aest': 0.15, 'P. menz': 0.15, 'O. ficu': 0.27, 'A. tequ': 0.27, 'S. bico': 0.15} #Hydraulic capacitance (MPa-1)

    f = .5;
    beta_r = 1.25; #Root respiration coefficient

    #Initial Conditions
    s0 = sInit;  #Initial Soil moisture
    Vwi =  .75*Vwt[species] #Initial Capacitance Depth (m3/m2 leaf area) 

    #C4 Constants from "Modeling of C3 and C4 photosynthesis under water-stressed conditions" Vico 2008
    gbs = .013 #Conductance between bundle sheath and mesophyll (mol m^-2s^-1)
    Vpmax0 = 120. #Maximum PEP carboxylase under well-watered conditions (umol/(m^2s))
    Vpr = 80. #PEP regeneration rate (umol/(m^2s))
    Kp = 80. #Michaelis-Menten coefficient for C4 species (umol/(mol))

    #Parameters for CAM model from Bartlett 2014
    Mmax = {'O. ficu': 190000000., 'A. tequ': 130000000.} #max concentration of malic acid (umol/m^3)
    #Mmax = {'O. ficu': 19., 'A. tequ': 130000000.} #max concentration of malic acid (umol/m^3)
    Ammax = {'O. ficu': 13.5, 'A. tequ': 11.1} #rate of malic acid storage flux (umol/(m^2 s)
    tr = 90.; #Relaxation time (m)
    c0 = 3000. #parameter for decarboxylation of malic acid (umol/mol)
    alpha_1 = 1/100. ;#Half reaction concentration for malic acid (umol/m^3)
    alpha_2 = 1/7.; #Half reaction concentration for malic acid (umol/m^3)
    alpha_3 = 1/7.; 
    k = .003 ;
    Topt1 = 288.65; #(K)
    Vcm = 0.0027 ; #Value controlling relative storage of malate (m)
    #Vcm = Vcm*10.**7
    mu = .5; #Circadian oscillator constant
    beta = 2.764; #Circadian oscillator constant
    c1 = .365; #Circadian oscillator constant
    c2 = .55; #Circadian oscillator constant
    c3 = 10.; #Circadian oscillator constant
    z0 = .55 ; #Initial value of z
    M0 = 0. #Initial Malic Acid Carbon Concentration (umol/m^3)

    #TimeStep Variables
    timestepM = 10.; #Change in time at each step (min)

    #Rainfall parameters
    alpha = 1.5; #Average Precipitaion Amount (cm)
    lambda_r = .3; #Frequency of Rainfall events (1/d)
    gamma = (n[sType]*Zr*100.)/alpha; #Normalized Depth of Rainfall
    EvMax = 1. #Maximum soil evaporation (mm/d)

    #General Constants
    lambda_w = 2.5*10**6; #Latent heat of water vaporization (J/kg )
    lambda_l = 550.*10**-9; #Wavelength of light (m)
    h = 6.63*10**-34; #Plank's constant (J s)
    cc = 3.00*10**8 ;#Speed of light (m/s)
    Ep = (h*cc)/lambda_l ;#Energy of photon (J)
    g = 9.8; #Gravity (m/s^2)
    rho_w = 998.; #Water density (kg/m^3)
    cp = 1012.; #Specific Heat of Air (J/(kg K))
    Ra = 290.; #Specific Gas Constant for Air (J/(kg K))
    NA = 6.022*10**23; #Avogadro's Constant (1/mol)
    R = 8.314; #Universal Gas Constant (J/(mol K))
    vw = 18.02/(rho_w*1000.) #Molar volume of water (mol/m3)

    #Atmospheric parameters
    po = 101.325*10**3; #Atmospheric pressure (Pa)
    rho = 1.27 ;#Air Density (kg/m^3)
    ca = 350.; #Atmopsheric CO2 concentration (ppm)
    gamma_w = (po*cp)/(.622*lambda_w);#Psychrometric constant for Penman Equation (J/(K-m^3))
    gammaTheta  = 4.78; #(K/km)
    B = 0.; #Empirical factor for boundary layer growth
    asat = 613.75;
    bsat = 17.502;
    csat = 240.97;
    rho_v = 0.87; #Density of water vapor (kg/m3)

    #Parameters for windspeed
    kvc = 0.41 #von Karman constant
    uz = 2.  # windspeed (m/s)
    zwind = 2. # height at which windspeed is measured (m)

    #Initial atmospheric conditions
    ho = 120.; #Initial Height of Boundary Layer (m)
    Tio = 280.; #Initial Temperature (K)

    #Declare Arrays for Data Storage
    ta = [];
    hgt = [ho];#Boundary layer height
    Tla = [Tio];#Leaf Temperature 
    H = [];#Sensible Heat Flux
    Ev = [];#Transpiration
    qsa = [];
    s = [s0];#Soil Moisture
    psi_ll = []; #Leaf Potential
    gsv = [];#Stomatal Conductance Pressure Deficit
    Ana = [];#Assimilated Carbon
    Asva =[]
    Avca = []
    Asca = []
    Rdca = []
    Rdva = []
    Ci = []; #Intercellular CO2 concentration
    Cm = []; #Mesophyll CO2 concentration
    Cs = [ca];#(*stomatal carbon concentration*)
    qa = []; #(*define array for specific humidity*)
    phiL=[];
    qia=[];
    gpa = [gpmax[species]];

    #Extraneous arrays for data storage
    psi_atma =[]
    VPDa = []; #Vapor Pressure Deficit

    # Arrays for CAM model
    if pType[species]=="CAM":
        Ma = [M0]
        za = [z0]
        Cc = []; #(*Cytoplasm Acid Carbon Concentration *)
        fDcO2 = []; #(*fdCO2 monitors the value of the air dryness function*)
    else:
        pass

    Th = 302.65; #High temperature for CAM model
    Tw = 283.15; #Low temperature for CAM model

    #Arrays for C4 model
    Cbs = [ca] #Bundle sheath cell CO2 concentration

    #Arrays for capacitance model
    psi_wa = []; #Array for potential of Storage
    Vwa = [Vwi]; #Array for stored water depth
    psi_xa = []; #Intermediate Potential at Storage Branch
    qwa = [];
    qsa = [];

    #Time Conversions
    def Steps(Duration, TimeStep):
        """Change Duration of Simulation to to number of timesteps according to timestep value"""
        return(Duration*24*60)/TimeStep
    def tcm(t):
        """TimeStep conversion, step counter to minutes"""
        return timestepM*t
    def tch(t):
        """TimeStep conversion, step counter to hours"""
        return tcm(t+1.)/60.
    def td(days):
        """Number of timesteps (days)"""
        return (days*60.*24.)/timestepM
    dt = timestepM*60. #no. of seconds in timestep, used to advance differential equations

    #SOLAR RADIATION
    phi_max = 500. #Maximum solar radiation (W/m^2);
    delta = 12. #day length (h)
    to = 6. #(h);
    t1=6.;
    t2=18.;
    def phi(t):
        """Leaf available energy, takes input of t in terms of hours"""
        if lightOption == "PAR":
            return max((4.*phi_max)/delta**2.*(-(t%24.)**2. + (delta + 2.*to)*(t%24.) - to*(to + delta)), 0.)
        elif lightOption=="STEP": #if lightOption=="STEP"
            if 0. <= (t%24.) < t1:
                return 0.
            elif t1 <= (t%24.) <= t2:
                return phi_max
            else: #{0, t2 <= Mod[t, 24] < 24}
                return 0.       
        elif lightOption=="ACTUAL":
                #return phiInp[int(t)] #int(t) rounds down to nearest hour integer
                return phiInp[int((t-1)*60/timestepM)] 
     
    #Atmospheric boundary layer
    def theta(h):
        """Potential temp at top of boundary layer (K), h input in m"""
        return gammaTheta*h/1000. + 293.6
    def q(h):
        """Relative Humidity at Top of Boundary layer (kg/kg), h input in m"""
        return -.00285*h/1000 + .01166
    def p(h):
        """Pressure Change with Altitude"""
        return po*(1. - 2.25577*10**-5*h)**5.25588
    def theta_c(Ta, h):
        """Potential Temperature Calculation from Atmopheric Temperature"""
        return Ta*(po/p(h))**(Ra/cp)
    def Temp(theta, h):
        """Temperature Calculation from Potential Temperature"""
        return theta*(p(h)/po)**(Ra/cp)
    theta_a = [theta_c(Tio, ho)] #Initial potential Temperature
    def esat(T):
        """Saturated vapor pressure (Pa)"""
        return asat*exp((bsat*(T - 273.))/(csat + T - 273.))
    def deltaS(T):
        return (bsat*(csat+T-273.)-bsat*(T-273.))/(csat+T-273.)**2.*esat(T)

    qio = .622*esat(Tio)/po #Initial specific humidity, assuming a starting value at saturation*)
    qaNight = .0082041;#Nighttime specific humidity
    qaDay = .0082041 + .0004293;#(*daytime specific humidity*)
    TaDay = 303.15;
    TaNight = 288.15;
    
    def TaSimple(t):
        return TaNight + (TaDay - TaNight)*phi(tch(t))/phi_max
    def qaSimple(t):
        return qaNight + (qaDay - qaNight)*phi(tch(t))/phi_max
    def TaStep(t):
        if phi(tch(t))>0:
            return TaDay
        else:
            return TaNight
    def qaStep(t):
        if phi(tch(t))>0:
            return qaDay
        else:
            return qaNight
        
    if weatherOption == "BLM":
        qa = [qio] #(*define array for specific humidity*)
        Ta = [Tio]#(*Temperature Atmosphere*)
    else: #if weatherOption == "AVG": or "STEP" or "CONST"
        qa = []
        Ta = []

    def gaCalc(h):
        """calculate atmospheric conductance (mm/s) from plant height (m)"""
        return 1000.*kvc**2.*uz/log((zwind-0.64*h)/(0.13*h))**2.
    def VPD(T, q):
        """Vapor pressure deficit (Pa)"""
        return esat(T) - (q*po)/.622
    def Drh(T, q):
        """Specific humidity deficit (kg/kg)"""
        return (esat(T)*.622)/po - q 
    def psi_atm(T, q):
        """Atmospheric water potential (MPa)"""
        return R*T/vw*log(q*po/.622/esat(T))/1000000.
    def qi(T, psi):
        """Specific humidity internal to leaf (kg/kg)"""
        return .622*esat(T)/po*exp(psi*1000000.*vw/R/T)
    def hgtNew(t):
        """Boundary layer height (m)"""
        if phi(tch(t)) > 0.:
           return ((1. + 2.*B)*H[t]*1000.)/(rho*cp*hgt[t]*gammaTheta)*dt + hgt[t]
        else:
           return ho
    def theta_aNew(t):
        """Change in temperature"""
        if phi(tch(t)) > 0.:
            return theta_a[t] + dt*H[t]/(rho*cp*hgt[t]) + (theta(hgt[t]) - theta_a[t])/hgt[t]*(hgt[t + 1] - hgt[t])
        else:
            return theta_c(Tio, ho)#K
    def qaNew(t): 
        """Specific Humidity (kg/kg)"""
        if phi(tch(t)) > 0.:
            return  qa[t] + (rho_w*Ev[t]/10.**6*dt)/(rho*hgt[t]) + (q(hgt[t]) - qa[t])/hgt[t]*(hgt[t + 1] - hgt[t])
        else:
            return qio
    def rain(i):
        """Rainfall (unitless)"""
        if rainOption == 'ACTUAL':
            return rInp[i];
        elif rainOption == 'STOCH':
            if np.random.random() > lambda_r*timestepM/(60.*24):
                return 0.
            else:
                return np.random.exponential(1./gamma)
        else:
            return 0.
    def HNew(Tl, Ta):
        """Sensible heat flux (W/m^2), per unit ground area"""
        return cp*rho*ga[species]*(Tl-Ta)/1000.*LAI[species]
    def HNewPM(t):
        """Sensible Heat (W/m^2)"""
        return phi(tch(t)) - lambda_w*rho_w*Ev[t]/10.**6 

    def TlNew(Ta, H):
        """Leaf temperature (K)"""
        if pType[species] == "CAM":
            return Ta + (H*10.)/(cp*rho*ga[species]) 
        else:
            return Ta + (H*1000.)/(cp*rho*ga[species])

    #Soil and plant water fluxes
    def L(s):
        """Leakage (um/s) """
        return .11574*Ks[sType]*s**(2.*b[sType] + 3.)                                   
    def psi_s(s):
        """Soil Potential (MPa)"""
        return psi_ss[sType]*(s**-b[sType])  
    def RAI(s):
        """Root area index"""
        return RAIW[species]*s**-a
    def gsr(s):
        """Soil-Root Conductance, per unit ground area (um/(s-MPa))"""
        return (L(s)*sqrt(RAI(s))*1000000.)/(float(pi)*g*rho_w*Zr)
    def gp(psi_l):
        """Plant conductance, per unit leaf area (um/(s-MPa))"""
        #return gpmax[species]
        if psi_l<-10:
            return 0.
        else:
            return gpmax[species]*exp(-(-psi_l/2.)**2.)
    def grsp(s, gp):
        """Soil-Root-Plant Conductance, per unit ground area (um/(s-MPa))"""
        return (LAI[species]*gsr(s)*gp)/(gsr(s) + LAI[species]*gp)
    def grsfp(s, t):
        """Soil-root-plant fraction conductance, per unit ground area (um/(s-MPa))"""
        return (LAI[species]*gsr(s)*gpa[t]/f)/(gsr(s) +  LAI[species]*gpa[t]/f) 
    def gw(Vw):
        """Xylem-storage conductance, per unit leaf area (um/(MPa-s))"""
        return gwmax[species]*(Vw/Vwt[species])**4. 
    def EVMx(s): 
        """Soil evaporation rate, per unit ground area (mm/day)"""
        if s > sh[sType]:
            return EvMax*(s - sh[sType])/(1. - sh[sType])
        else:
            return 0.
    def sNew(i): 
        """Soil moisture"""
        if capOn == True:
            return (dt/(n[sType]*Zr*10.**6)*(-Ev[i] +LAI[species]*gw(Vwa[i])*((psi_wf(Vwa[i]))-Ev[i]*(1.-f)/(LAI[species]*gpa[i])-psi_ll[i]) - (EVMx(s[i])*1000.)/(24.*60*60)- L(s[i]))) + s[i]
        else:
            return (dt/(n[sType]*Zr*10.**6)*(-Ev[i] - (EVMx(s[i])*1000.)/(24.*60*60)- L(s[i]))) + s[i]

    #DEFINITIONS FOR CAPACITANCE
    def psi_wf(Vw): 
        return (1./cap[species])*Vw/Vwt[species] - (1./cap[species])
    def Vwf(i):
        """Stored water volume, per unit leaf area (m3/m2)"""
        return min(Vwa[i] - gw(Vwa[i])*((psi_wf(Vwa[i])) - (Ev[i]*(1. - f))/(LAI[species]*gpa[i]) - psi_ll[i])*dt/10.**6, Vwt[species])
    def psi_xf(i):
        return Ev[i]*(1. - f)/(LAI[species]*gpa[i]) + psi_ll[i]
    def qwf(i):
        """Stored water flux, per unit ground area"""
        return (Vwa[i] - Vwf(i))*LAI[species]*10.**6/dt
    def qsf(i):
        """Soil water flux, per unit ground area"""
        return Ev[i] - qwf(i)
    def Vcmax(Tl):
        """Maximum carboxylation rate (umol/(m^2s))"""
        return Vcmax0[species]*exp(HaV[species]/(R*To)*(1. - To/Tl))/(1. + exp((SvC[species]*Tl - HdV[species])/(R*Tl)))
    def Gamma(Tl):
        """CO2 compensation point (umol/mol)"""
        return gamma_o[species]*(1. + gamma_1*(Tl - To) + gamma_2*(Tl - To)**2.);
    def Jmax(Tl):
        """Max. e- transport rate (umol/(m^2s))"""
        return Jmax0[species]*exp(HaJ[species]/(R*To)*(1. - To/Tl))/(1. + exp((SvQ[species]*Tl - HdJ[species])/(R*Tl))) 
    def J(phi, Tl):
        """Electron transport rate (umol/(m^2s))"""
        return min((phi*10.**6)/(Ep*NA)*kappa_2*.5, Jmax(Tl)) 
    def Ko(Tl):
        """Michaelis-menten coefficient for O2"""
        return Ko0[species]*exp(Hko[species]/(R*To)*(1. - To/Tl))
    def Kc(Tl):
        """Michaelis-menten coefficient for CO2"""
        return Kc0[species]*exp(Hkc[species]/(R*To)*(1. - To/Tl))
    def Ac(ci, Tl):
        """Rubisco-limited photosynthetic rate (umol/(m^2s^1))"""
        return Vcmax(Tl)*(ci - Gamma(Tl))/(ci + Kc(Tl)*(1. + (oi[species]*1000.)/Ko(Tl)))
    def Aq(phi, ci, Tl):
        """Light-limited photosynthetic rate (umol/(m^2s^1))"""
        return (J(phi, Tl)*(ci - Gamma(Tl)))/(4.*(ci + 2.*Gamma(Tl)))
    def AphiciTl(phi, ci, Tl):
        """Net photosynthetic demand for CO2 (umol/(m^2s^1))"""
        return min(Ac(ci, Tl), Aq(phi, ci, Tl))
    def Apsilc02(psi_l):  
        """Vulnerability curve for water potential"""
        if psi_l < psi_laoO[species] :
            return 0.
        elif psi_laoO[species] <= psi_l <= psi_la1O[species] :
            return (psi_l - psi_laoO[species] )/(psi_la1O[species]  - psi_laoO[species])
        else: #if \[Psi]l > \[Psi]la1O 
            return 1.
        
        
    #C3 and CAM models diverge:
    def An(phi, psi_l, Tl, ci, t): 
        """Photosynthetic rate, per unit leaf area (umol/(m^2s))"""
        if pType[species]=="CAM":
            M=Ma[t]
            z=za[t]
            return Apsilc02(psi_l)*(AphiciTl(phi, ci, Tl) - Rdc(phi, Tl))*(1. - fTM(z, M)) + Asv(phi, Tl, psi_l, z, M) #(*Flux Asc + Asv*)
            #return Apsilc02(psi_l)*(AphiciTl(phi, ci, Tl) )*(1. - fTM(z, M)) + Asv(phi, Tl, psi_l, z, M) 
        else: #same function for both C3 and C4. no respiration 
            return Apsilc02(psi_l)*AphiciTl(phi, ci, Tl)
            #return Apsilc02(psi_l)*(AphiciTl(phi(tch(t)), ci, Tl) - Rdf(Tl))) # with dark respiration
    def Asc(phi, psi_l, Tl, ci, t):
        """Flux Asc from stomata to Calvin cycle"""
        M=Ma[t]
        z=za[t]
        return Apsilc02(psi_l)*(AphiciTl(phi, ci, Tl) - Rdc(phi, Tl))*(1. - fTM(z, M)) #(*Flux Asc + Asv*) 
    def Rroot(Tl):
        """Root respiration (umol/(m^2s))"""
        return beta_r*fA(Tl)*fTI(Tl)*froot[species]*B[t]
    def Rgrowth(Tl, t, psi_l, ci):
        """Growth respiration (umol/(m^2s))"""
        return 0.2*An(phi(tch(t)), psi_l, Tl, ci, t)
    def Bm(t):
        """Biomass (umol/m2)"""
        return dt*(An-Rroot(Tl)-Rgrowth(Tl,t, psi_l, ci))+Ba[t]
    def Rdf(Tl):
        """Dark respiration flux (umol/(m^2s))"""
        return Rd[species]*exp(HkR[species]/(R*To)*(1. - To/Tl))
    def Rdv(phi, Tl):
        """Flux of dark respiration to vacuole (umol/(m^2s))"""
        return Rdf(Tl)*exp(-phi)
    def Rdc(phi, Tl):
        """Flux of dark respiration to calvin cycle (umol/(m^2s))"""
        return Rdf(Tl)*(1. - exp(-phi))
    def fT(z):
        """Circadian order function"""
        return exp(-(z/mu)**c3)
    def fTM(z, M):
        """Carbon circadian control"""
        return (1. - fT(z))*M/(alpha_1*Mmax[species] + M)
    def Asv(phi, Tl, psi_l, z, M):
        """Flux from stomata to vacuole (umol/(m^2s))"""
        if phi>0 and M < .005:
            return 0.
        else:
            if Mmax[species]*((Th - Tl)/(Th - Tw)*(1. - alpha_3) + alpha_3) > M and (1. - k*(Tl - Topt1)**2.) >0:
                return (Ammax[species]*(1. - k*(Tl - Topt1)**2.) - Rdv(phi, Tl))*fT(z)*(Mmax[species]*((Th - Tl)/(Th - Tw)*(1. - alpha_3) + alpha_3) - M)/(alpha_2*Mmax[species]*((Th - Tl)/(Th - Tw)*(1. - alpha_3) + alpha_3) + (Mmax[species]*((Th - Tl)/(Th - Tw)*(1. - alpha_3) + alpha_3) - M))*Apsilc02(psi_l)
            else:
                return 0.
    def Avc(phi, cc, Tl, z, M):
        """Flux from vacuole to calvin cycle (umol/(m^2s))"""
        return (AphiciTl(phi, cc, Tl) - Rdc(phi, Tl))*fTM(z, M)
    def MT(z, Tl, phi): 
        """Malic acid equilibrium value"""
        if phi>0.:
            return Mmax[species]*(c1*((Th - Tl)/(Th - Tw) + 1.)*(beta*(z - mu))**3. - beta*(Th - Tl)/(Th - Tw)*(z - mu) + c2*(Th - Tl)/(Th - Tw)) 
        else:
            return Mmax[species]*(c1*((Th - Tl)/(Th - Tw) + 1.)*(beta*(z - mu))**3. - beta*(Th - Tl)/(Th - Tw)*(z - mu) + c2*(Th - Tl)/(Th - Tw)+ 50.*z**6.) 
    def zNew(phi, M, z, Tl):
        return max(0, dt*(M - MT(z, Tl, phi))/(Mmax[species]*60.*tr) + z)
    def fD(vpd):
        """Stomatal response to vapor pressure deficit"""
        return 3/13./sqrt(vpd/1000.)
    def gsN(phi, Ta, psi_l, qa, Tl, ci, t): 
        """Stomatal conductance to CO2, per unit leaf area (mol/m2/s)"""
        if An(phi, psi_l, Tl, ci, t) < 0.:
            return 0.
        else:
            return a1new[species]*An(phi, psi_l, Tl, ci, t)/ca*fD(VPD(Ta, qa))
    def gwN(phi, Ta, psi_l, qa, Tl, ci, t): 
        """Stomatal conductance to water, per unit leaf area (mol/m2/sec)"""
        #return gsN(phi, Ta, psi_l, qa, Tl, ci, t)*(1.6*(1. + gmgsratio[pType[species]]))/(1.6 + gmgsratio[pType[species]]) + (gcut[species]*po/(1000.*R*Ta))
        return gsN(phi, Ta, psi_l, qa, Tl, ci, t)*1.6 + (gcut[species]*po/(1000.*R*Ta))
    def Evf(phi, Ta, psi_l, qa, Tl, ci, t):
        """Transpiration, per unit ground area (um/sec)"""
        return LAI[species]*(1./(gwN(phi, Ta, psi_l, qa, Tl, ci, t)*R*Ta/po*1000000.)+1./(ga[species]*1000.))**(-1.)*rho/rho_w*(qi(Tl, psi_l)-qa)
    def Evfpmc(phi, Ta, psi_l, qa, Tl, ci, t): 
        """Transpiration, Penman-Monteith Equation (um/sec)"""
        return ((lambda_w*gamma_w*ga[species]/1000.*rho*Drh(Ta, qa) + deltaS(Ta)*phi)*R*Ta/po*gwN(phi, Ta, psi_l, qa, Tl, ci, t)*1000000.*LAI[species])/(rho_w*lambda_w*(gamma_w*(ga[species]/1000. + (R*Ta)/po*gwN(phi, Ta, psi_l, qa, Tl, ci, t)*LAI[species]) + (R*Ta)/po*gwN(phi, Ta, psi_l, qa, Tl, ci, t)*LAI[species]*deltaS(Ta)))
    def CsNew(i):
        """CO2 concentration at leaf surface (ppm)"""
        return ca - Ana[i]/ga[species]
    def CiNew(i):
        """CO2 concentration in mesophyll cytosol (ppm)"""
        #return Rc[pType[species]]*ca  
        return Cs[i]*(1.-1./(a1new[species]*fD(VPD(Ta[i], qa[i])))) 
    # def CiN(i):
    #     """CO2 concentration in stomatal cavity (ppm), Norman model"""
    #     return Rc[pType[species]]*ca   
    def CmNew(i):
        return CiNew(i)
    #     """CO2 concentration in mesophyll (ppm)"""
    #     return CiNew(i) - ((1.+VPD(Ta[i], qa[i])/Dxnew[species])*Cs[i])/(a1new[species]*gmgsratio[pType[species]])   
    def CcNew(i, z, M):
        """CO2 concentration in mesophyll cytosol resulting from malic acid decarboxylation (ppm)"""
        return CmNew(i) + fTM(z, M)*c0
    def MNew(i, psi_l, cc, Tl, z, M): 
        """Malic acid concentration"""
        return max(((dt/ Vcm)*(Asv(phi(tch(i)), Tl, psi_l, z, M) - Avc(phi(tch(i)), cc, Tl, z, M) + Rdv(phi(tch(i)), Tl))) + M, 0.)
    def gsvNew(phi, Ta, psi_l, qa, Tl, ci, t):
        """Stomatal conductance to water vapor, per unit leaf area (mm/s)"""
        return (R*Ta)/po*1000.*gwN(phi, Ta, psi_l, qa, Tl, ci, t)

    #DEFINITIONS FOR C4
    def Vp(ci):
        return min((ci*Vpmax0)/(ci + Kp), Vpr)
    def CbsNew(i, psi_l):
        """CO2 concentration in bundle sheath cell (ppm)"""
        return ((Vp(Cm[i]) - AphiciTl(phi(tch(i)), Cbs[i], Tla[i]))*Apsilc02(psi_l))/gbs + Cm[i]


    for i in range(Steps(Duration, int(timestepM))):
        #time=i*dt/3600./24.; #time in days
        #ta.append(time)
        phiL.append(phi(tch(i)))
        if weatherOption == "AVG":
            Ta.append(TaSimple(i))
            qa.append(qaSimple(i))
        elif weatherOption == "STEP":
            Ta.append(TaStep(i))
            qa.append(qaStep(i))
            #qa.append(qaSimple(i))
        elif weatherOption == "CONST":
            #Ta.append(TaNight)
            Ta.append(TaStep(i))
            qa.append(qaNight)
        elif weatherOption == "ACTUAL":
            Ta.append(taInp[i])
            qa.append(qaInp[i])
             #Ta.append(taInp[int(i*timestepM/60)])
             #qa.append(qaInp[int(i*timestepM/60)])
        else: #weatherOption == "BLM"
             pass
        
        VPDa.append(VPD(Ta[i], qa[i]))
        psi_atma.append(psi_atm(Ta[i], qa[i]))

        if pType[species]=="CAM":
            Ci.append(CiNew(i))
            Cm.append(CmNew(i))
            Cc.append(CcNew(i, za[i], Ma[i])) 
        else:
            Ci.append(CiNew(i))
            Cm.append(CmNew(i))

        if pType[species]=="C4":
            C1 = Cbs[i]
        else:
            C1 = Cm[i]

        if transpirationType == "PM":
            def evBal(psi_l, i):
                if pType[species]=="C4":
                    C1 = Cbs[i]
                else:
                    C1 = Ci[i]
                if capOn == True:
                    return Evfpmc(phi(tch(i)), Ta[i], psi_l, qa[i], Tla[i], C1, i) - \
                           (grsfp(s[i], i)*(psi_s(s[i]) - psi_l) + \
                            LAI[species]*gw(Vwa[i])*(psi_wf(Vwa[i]) - psi_l))/(1. + (grsfp(s[i], i)*(1. - f))/(LAI[species]*gpa[i]) + (gw(Vwa[i])*(1. - f))/gpa[i])
                else:
                    return Evfpmc(phi(tch(i)), Ta[i], psi_l, qa[i], Tla[i], C1, i) - grsp(s[i], gpa[i])*(psi_s(s[i]) - psi_l)
            psi_lNew = fsolve(evBal, -3., args=(i))
            psi_ll.append(psi_lNew[-1])
            gpa.append(gp(psi_lNew))
            qia.append(qi(Tla[i], psi_ll[i]))
            Ev.append(Evfpmc(phi(tch(i)), Ta[i], psi_ll[i], qa[i], Tla[i], C1, i))
        else:
            def fBal(p, i):
                psi_l, Tl =p
                #if heatCap ==True:
                #    return (phi(tch(i)) - HNew(Tl, Ta[i]) -lambda_w*rho_w*Evf(phi(tch(i)), Ta[i], psi_l, qa[i], Tl, C1, i)/1000000. - cpw*Vwt[species]*LAI[species]*1000./dt*(Tl-Tla[i]), Evf(phi(tch(i)), Ta[i], psi_l, qa[i], Tl, C1, i) - grsp(s[i], gpa[i])*(psi_s(s[i]) - psi_l))
                if capOn==True:
                    return (Evf(phi(tch(i)), Ta[i], psi_l, qa[i], Tl, C1, i)\
                           -(grsfp(s[i], i)*(psi_s(s[i]) - psi_l) + LAI[species]*gw(Vwa[i])*(psi_wf(Vwa[i]) - psi_l))/(1. + (grsfp(s[i], i)*(1. - f))/(LAI[species]*gpa[i]) + (gw(Vwa[i])*(1. - f))/gpa[i]),\
                            phi(tch(i)) - HNew(Tl, Ta[i]) -lambda_w*rho_w*Evf(phi(tch(i)), Ta[i], psi_l, qa[i], Tl, C1, i)/1000000.)
                else:
                    return (phi(tch(i)) - HNew(Tl, Ta[i]) -lambda_w*rho_w*Evf(phi(tch(i)), Ta[i], psi_l, qa[i], Tl, C1, i)/1000000., Evf(phi(tch(i)), Ta[i], psi_l, qa[i], Tl, C1, i) - grsp(s[i], gpa[i])*(psi_s(s[i]) - psi_l))

            if len(psi_ll) == 0:
                psi_lNew, TlNew = fsolve(fBal, (-1., 305.), args=(i))
            else:
                psi_lNew, TlNew = fsolve(fBal, (psi_ll[-1], 305.), args=(i))

            psi_ll.append(psi_lNew)
            Tla.append(TlNew)
            qia.append(qi(TlNew, psi_lNew))
            gpa.append(gp(psi_lNew))
            Ev.append(Evf(phi(tch(i)), Ta[i], psi_ll[i], qa[i], Tla[i], C1, i))
        

        if rainOption == 'CONST':
            s.append(s0);
        else:
            s.append((min(1., sNew(i) + rain(i) )))

        
        if transpirationType =="PM":
            H.append(HNewPM(i))
            Tla.append(TlNew(Ta[i], H[i]))
        else: 
            H.append(HNew(Tla[i], Ta[i]))

        hgt.append(hgtNew(i))
        theta_a.append(theta_aNew(i))

        if weatherOption == "BLM":
            Ta.append(Temp(theta_a[i+1], ho))
            qa.append(qaNew(i))
        else:
            pass
            
        Ana.append(An(phi(tch(i)), psi_ll[i], Tla[i], C1, i))
        Asva.append(Asv(phi(tch(i)), Tla[i], psi_ll[i], za[i], Ma[i]))
        Avca.append(Avc(phi(tch(i)), C1, Tla[i], za[i], Ma[i]))
        Asca.append(Asc(phi(tch(i)), psi_ll[i], Tla[i], C1, i))
        Rdca.append(Rdc(phi(tch(i)), Tla[i]))
        Rdva.append(Rdv(phi(tch(i)), Tla[i]))
        Cs.append(CsNew(i))
        gsv.append(gsvNew(phi(tch(i)), Ta[i], psi_ll[i], qa[i], Tla[i], C1, i))
        Cbs.append(CbsNew(i, psi_ll[i]))

        #CAM MODEL
        if pType[species] == "CAM":
            Ma.append(MNew(i, psi_ll[i], Cc[i], Tla[i], za[i], Ma[i]))
            za.append(zNew(phi(tch(i)), Ma[i], za[i], Tla[i])) 
        else:
            pass

        #CAPACITANCE MODEL
        if capOn == True:
            psi_wa.append(psi_wf(Vwa[i])) 
            psi_xa.append(psi_xf(i)) 
            Vwa.append(Vwf(i))
            qwa.append(qwf(i)) 
            qsa.append(qsf(i))
        else:
            pass

    del Tla[-1]
    del Vwa[-1]
    del hgt[-1]
    del s[-1]
    del gpa[-1]
    if pType[species] =="CAM":
        del Ma[-1]
        del za[-1]
    else:
        pass

    if pType[species] =="CAM" and capOn ==True:
        output = np.transpose([gsv, [psi_s(x) for x in s], psi_ll, psi_atma, gpa, s, [x*24*3.6 for x in Ev], Tla, [x*24*3.6/1000 for x in Ana], qa, qia, [x/Vwt[species] for x in Vwa], VPDa, Ci, Cm, Cc, Ma, za, Asva, Avca, Asca, Rdca, Rdva])
        output = pd.DataFrame(output, columns = ['gsv', 'psi_s', 'psi_ll', 'psi_atma', 'gp', 's', 'Ev', 'Tl', 'An', 'qa', 'qi', 'w', 'VPD', 'Ci', 'Cm', 'Cc', 'M', 'z', 'Asv', 'Avc', 'Asc', 'Rdc', 'Rdv'])
    elif pType[species] =="CAM":
        output = np.transpose([gsv, [psi_s(x) for x in s], psi_ll, psi_atma, gpa, s, [x*24*3.6 for x in Ev], Tla, [x*24*3.6/1000 for x in Ana], qa, qia, VPDa, Ci, Cm, Cc, Ma, za])
        output = pd.DataFrame(output, columns = ['gsv', 'psi_s', 'psi_ll', 'psi_atma', 'gp', 's', 'Ev', 'Tl', 'An', 'qa', 'qi', 'VPD', 'Ci', 'Cm', 'Cc', 'M', 'z'])
    elif capOn==True:
        output = np.transpose([gsv, [psi_s(x) for x in s], psi_ll, psi_atma, gpa, s, [x*24*3.6 for x in Ev], Tla, [x*24*3.6/1000 for x in Ana], qa, qia, [x/Vwt[species] for x in Vwa], VPDa, Ci, Cm])
        output = pd.DataFrame(output, columns = ['gsv', 'psi_s', 'psi_ll', 'psi_atma', 'gp', 's', 'Ev', 'Tl', 'An', 'qa', 'qi', 'w', 'VPD', 'Ci', 'Cm'])
    else: 
        output = np.transpose([gsv, [psi_s(x) for x in s], psi_ll, psi_atma, gpa, s, [x*24*3.6 for x in Ev], Tla, [x*24*3.6/1000 for x in Ana], qa, qia, VPDa, Ci, Cm])
        output = pd.DataFrame(output, columns = ['gsv', 'psi_s', 'psi_ll', 'psi_atma', 'gp', 's', 'Ev', 'Tl', 'An', 'qa', 'qi', 'VPD', 'Ci', 'Cm'])

    #outputs s (%), Ev (mm/d), An (mol/m2/d), gs (mm/s) 
    return output
