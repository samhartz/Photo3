#Soil parameters
psi_ss = {'Sand': -.34*10**-3, 'Sandy loam': -.7*10**-3, 'Loamy sand': -.17*10.**-3.,'Loam': -1.43*10.**-3., 'Clay': -1.82*10**-3.} # (MPa)
b = {'Sand': 4.05,'Sandy loam': 4.9, 'Loamy sand':4.38, 'Loam': 5.39, 'Clay': 11.4} # (-)
Ks = {'Sand': 200, 'Sandy loam': 80.,'Loamy sand':100., 'Loam': 20., 'Clay': 1.} # Saturated hydraulic condutivity (cm/day)
n = {'Sand': .35, 'Sandy loam': .43, 'Loamy sand':.42,'Loam': .45, 'Clay': .5} # Porosity
sh = {'Sand': .08, 'Sandy loam':.14, 'Loamy sand':.08,'Loam': .19, 'Clay': 0.47} # Hygroscopic point

Zr = {'T. aest': 0.75, 'P. menz': 0.65, 'O. ficu': 0.3, 'A. tequ': 0.3, 'S. bico': 0.5} # Rooting depth (m)

#Plant photosynthetic parameters

pType = {'T. aest': 'C3', 'P. menz': 'C3', 'O. ficu': 'CAM', 'A. tequ': 'CAM', 'S. bico': 'C4'} # Photosynthetic type
LAI = {'T. aest': 5., 'P. menz': 8.4, 'O. ficu': 4.4, 'A. tequ': 6., 'S. bico': 5.} # Leaf area index (m2/m2)

Vcmax0 = {'T. aest': 83.,   'P. menz': 57.7,   'O. ficu': 13., 'A. tequ': 19.5, 'S. bico': 39.} # Max. carboxylation rate (umol/(m^2s^2))
Hkc =    {'T. aest': 59430.,  'P. menz': 59430., 'O. ficu': 54930., 'A. tequ': 59430., 'S. bico': 59430.} # Activation Energy for Kc (J/mol)
Hko =    {'T. aest': 36000.,  'P. menz': 36000., 'O. ficu': 36000., 'A. tequ': 36000., 'S. bico': 36000.} # Activation Energy for Ko (J/mol)
HkR =    {'T. aest': 53000.,  'P. menz': 53000., 'O. ficu': 53000., 'A. tequ': 53000., 'S. bico': 53000.}  # Activation Energy for Rd (J/mol)
Rd0 =     {'T. aest': 4.93,    'P. menz': .32,    'O. ficu': .32, 'A. tequ': .32, 'S. bico': .32} # Standard Dark respiration at 25 C (umol/(m^2s))
HaV =    {'T. aest': 62000.,  'P. menz': 61210., 'O. ficu': 72000., 'A. tequ': 72000., 'S. bico': 72000.} # Activation Energy for Vc,max (J/mol)
HdV =    {'T. aest': 202900., 'P. menz': 200000.,'O. ficu': 200000., 'A. tequ': 200000., 'S. bico': 200000.} # Deactivation Energy for Vc,max (J/mol)
HaJ =    {'T. aest': 50000.,  'P. menz': 37870., 'O. ficu': 50000., 'A. tequ': 50000., 'S. bico': 50000.} # Activation Energy for Jmax (J/mol)
HdJ =    {'T. aest': 200000., 'P. menz': 200000.,'O. ficu': 200000., 'A. tequ': 200000., 'S. bico': 200000.} # Deactivation Energy for Jmax (J/mol)
Jmax0 =  {'T. aest': 132.,   'P. menz': 98.5,   'O. ficu': 2.*Vcmax0['O. ficu'], 'A. tequ': 2.*19.5, 'S. bico': 180.} # Max. e- transport rate (umol/(m^2s))

a1new = {'T. aest': 15., 'P. menz': 15., 'O. ficu': 0.6*15., 'A. tequ': 0.6*15., 'S. bico': 0.5*15.} # Parameter relating gsv to an (-)

Kc0 = {'T. aest': 302., 'P. menz': 302., 'O. ficu': 302., 'A. tequ': 302., 'S. bico': 302.} # Michaelis constant for C02 at TO (umol/mol)
Ko0 = {'T. aest': 256., 'P. menz': 256., 'O. ficu': 256., 'A. tequ': 256., 'S. bico': 256.} # Michaelis constant for 02 at TO (mmol/mol)
oi =  {'T. aest': .209, 'P. menz': .209, 'O. ficu': .209, 'A. tequ': .209, 'S. bico': .209} # Oxygen Concentration (mol/mol)
SvC = {'T. aest': 649., 'P. menz': 649., 'O. ficu': 649., 'A. tequ': 649., 'S. bico': 649.} # Entropy term for carboxylation (J/mol)
SvQ = {'T. aest': 646., 'P. menz': 646., 'O. ficu': 646., 'A. tequ': 646., 'S. bico': 646.} # Entropy term for e-transport (J/mol)

TO = 293.2 # Reference Temperature for photosynthetic parameters (K)
KAPPA_2 = .3 # Quantum yield of photosynthesis (mol CO2/mol photon)
gamma_o = {'T. aest': 34.6, 'P. menz': 34.6, 'O. ficu': 34.6, 'A. tequ': 34.6, 'S. bico': 10.}  # C02 compensation point at To (umol/mol)
GAMMA_1 = .0451 # Parameter for temp dependence of CO2 compensation point (1/K)
GAMMA_2 = .000347 # Parameter for temp dependence of CO2 compensation point (1/K^2)

# Plant hydraulic parameters
gcut = {'T. aest': 0.15*2., 'P. menz': 0.007, 'O. ficu': 0., 'A. tequ': 0., 'S. bico': 0.1802} # Cuticular conducatnce, per unit leaf area (mm/s)
ga = {'T. aest': 61., 'P. menz': 324., 'O. ficu': 324., 'A. tequ': 61., 'S. bico': 61.} # Atomospheric  Conductance, per unit ground area (mm/s)
gmgsratio = {'C3': 1.65, 'C4': 2.65, 'CAM': 1.}; #gm/gs (mol/(m^2-s)) 
Rc = {'C3': 0.7, 'C4': 0.4, 'CAM': 0.5}; #Ratio of ci:ca
RAIW = {'T. aest': 5.6, 'P. menz': 10., 'O. ficu': 3., 'A. tequ': 3., 'S. bico': 5.6} # Root area index (m2/m2)
A_ROOT = 8.; # Parameter accounting for root growth
gwmax = {'T. aest': 0., 'P. menz': 0.05, 'O. ficu': .02, 'A. tequ': .002, 'S. bico': 0.} # Conductance b/t psi_x and psi_w, per unit root area (um/(MPa s))
gpmax = {'T. aest': 11.7, 'P. menz': 0.056, 'O. ficu': .4, 'A. tequ': .04, 'S. bico': 0.13} # Plant conductance, per unit leaf area (um/(MPa s))
vwt = {'T. aest': .000001, 'P. menz': .027/LAI['P. menz'], 'O. ficu': .0113, 'A. tequ': .00415, 'S. bico': .000001} # Max. storage depth (m3/m2 leaf area) 
cap = {'T. aest': 0.15, 'P. menz': 0.15, 'O. ficu': 0.83, 'A. tequ': 0.27, 'S. bico': 0.15} # Hydraulic capacitance (MPa-1)
F_CAP = .5 # Parameter for plant hydraulic capacitance scheme

psi_la1O = {'T. aest': -.7, 'P. menz': -.5, 'O. ficu': -.5, 'A. tequ': -.5, 'S. bico': -.5} # Parameter describing vulnerability to leaf water potential (MPa)
psi_laoO = {'T. aest': -2., 'P. menz': -3., 'O. ficu': -3., 'A. tequ': -3., 'S. bico': -1.8} # Parameter describing vulnerability to leaf water potential (MPa)

# Parameters for C4 model
GBS = .013 #Conductance between bundle sheath and mesophyll (mol m^-2s^-1)
VPMAX0 = 120. #Maximum PEP carboxylase under well-watered conditions (umol/(m^2s))
VPR = 80. #PEP regeneration rate (umol/(m^2s))
KP = 80. #Michaelis-Menten coefficient for C4 species (umol/(mol))

# Parameters for CAM model 
Mmax = {'O. ficu': 190000000., 'A. tequ': 130000000.} # max concentration of malic acid (umol/m^3)
Ammax = {'O. ficu': 13.5, 'A. tequ': 11.1} # rate of malic acid storage flux (umol/(m^2 s)
tr = 90.; # Relaxation time for circadian oscillator (min)
c0 = 3000. # parameter for decarboxylation of malic acid (umol/mol)
ALPHA_1 = 1/100.
ALPHA_2 = 1/7. 
k = .003 
Topt1 = 288.65 # (K)
VCM = 0.0027 # Value controlling relative storage of malate (m)
MU = .5 # Circadian oscillator constant
BETA = 2.764 # Circadian oscillator constant
CIRC_1 = .365 # Circadian oscillator constant
CIRC_2 = .55 # Circadian oscillator constant
CIRC_3 = 10. # Circadian oscillator constant
z0 = .55  # Initial value of z (-)
M0 = 0. # Initial Malic Acid Carbon Concentration (umol/m^3)
TH = 302.65 # High temperature for CAM model (K)
TW = 283.15 # Low temperature for CAM model (K)

# TimeStep Variables
timestepM = 30 # Model change in time at each step (min)
timestepD = 30 # timestep of input data 
dt = timestepM*60. # no. of seconds in timestep, used to advance differential equations

# Rainfall parameters
alpha = 1.5 # Average Precipitaion Amount (cm)
lambda_r = .3 # Frequency of Rainfall events (1/d)
#gamma = (n[sType]*zr*100.)/alpha; # Normalized Depth of Rainfall
EVMAX = 1. # Maximum soil evaporation (mm/d)

# General Constants
LAMBDA_W = 2.5*10**6 # Latent heat of water vaporization (J/kg )
LAMBDA_L = 550.*10**-9 # Wavelength of light (m)
H = 6.63*10**-34 # Planck's constant (J s)
CC = 3.00*10**8 # Speed of light (m/s)
EP = (H*CC)/LAMBDA_L # Energy of photon (J)
g = 9.8 # Gravity (m/s^2)
RHO_W = 998. # Water density (kg/m^3)
CP_A = 1012. # Specific Heat of Air (J/(kg K))
R_A = 290. # Specific Gas Constant for Air (J/(kg K))
NA = 6.022*10**23 # Avogadro's Constant (1/mol)
R = 8.314; # Universal Gas Constant (J/(mol K))
VW = 18.02/(RHO_W*1000.) # Molar volume of water (mol/m3)

# Atmospheric parameters
P_ATM = 101.325*10**3 # Atmospheric pressure (Pa)
RHO_A = 1.27 # Air Density (kg/m^3)
ca = 350. # Atmopsheric CO2 concentration (ppm)
#gamma_w = (P_ATM*CP_A)/(.622*LAMBDA_W);#Psychrometric constant for Penman Equation (J/(K-m^3))
GAMMA_THETA  = 4.78 # (K/km)
B = 0. # Empirical factor for boundary layer growth
A_SAT = 613.75
B_SAT = 17.502
C_SAT = 240.97
RHO_V = 0.87 # Density of water vapor (kg/m3)
CP_W = 4184. # specific heat of water (J/kg/K)

# Initial atmospheric conditions
ho = 120. # Initial Height of Boundary Layer (m)
tio = 280. # Initial Temperature (K)

qaNight = .00802 # Nighttime specific humidity
qaDay = .00739
VPdiff = 2500.

TaDay = 300.
TaNight = 288.