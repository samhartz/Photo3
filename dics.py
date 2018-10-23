timestepM = 30 # Model change in time at each step (min)
timestepD = 30 # timestep of input data 
dt = timestepM*60. # no. of seconds in timestep, used to advance differential equations

#c
# Rainfall parameters
alpha = 1.5 # Average Precipitaion Amount (cm)
lambda_r = .3 # Frequency of Rainfall events (1/d)
#gamma = (n[sType]*zr*100.)/alpha; # Normalized Depth of Rainfall

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
#ca = 350. # Atmopsheric CO2 concentration (ppm)
#gamma_w = (P_ATM*CP_A)/(.622*LAMBDA_W);#Psychrometric constant for Penman Equation (J/(K-m^3))
GAMMA_THETA  = 4.78 # (K/km)
B = 0. # Empirical factor for boundary layer growth
A_SAT = 613.75
B_SAT = 17.502
C_SAT = 240.97
RHO_V = 0.87 # Density of water vapor (kg/m3)
CP_W = 4184. # specific heat of water (J/kg/K)

