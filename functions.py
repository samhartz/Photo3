from math import exp, pi, sqrt, log
from dics import *


def steps(duration, timeStep):
    """Change Duration of Simulation to to number of timesteps according to timestep value"""
    return (duration*24*60)//timeStep

def VPD(ta, qa):
    """Vapor pressure deficit (Pa)"""
    return esat(ta) - (qa*P_ATM)/.622

def esat(ta):
    """Saturated vapor pressure (Pa)"""
    return A_SAT*exp((B_SAT*(ta - 273.))/(C_SAT + ta - 273.))
def qaRh(rh, ta):
    """Specific humidity (kg/kg), input of rh in %, ta in K"""
    return 0.622*rh/100.*esat(ta)/P_ATM # needs to be in kg/kg

