from math import exp, pi, sqrt, log
from scipy.optimize import fsolve
from sympy import *
import numpy as np
from dics import *
from functions import *

class Soil(object):
	EVMAX = 3. # Maximum bare soil evaporation rate (mm/d)
	def __init__(self, stype, dynamics, zr, s):
		self.PSI_SS = stype.PSI_SS # Saturation suction (MPa)
		self.B = stype.B # Exponent of soil retention curve (-)
		self.KS = stype.KS # Saturated hydraulic conductivity (cm/d)
		self.N = stype.N # Soil porosity (-)
		self.SH = stype.SH # Hygroscopic point of soil (-)
		self.ZR = zr # Rooting depth (m)
		self.s = s # Soil moisture (-)
		self.s_a = [] # Soil moisture array (-)
		self.psi_s_a = [] # Soil water potential array (MPa)
		self.dynamics = dynamics # Rules according to which soil moisture evolves
		self.rain_amt = 0 # Rain amount (mm)
		self.sm_inp = s
	def update(self, dt, zr, qs):
		self.s = self.dynamics.snew(self, dt, zr, qs)
		#self.s = (dt/(self.N*zr*10.**6)*(-qs - (self.evap(self.s)*1000.)/(24.*60*60)- self.leak(self.s))) + self.s
		self.s_a.append(self.s)
		self.psi_s_a.append(self.psi_s(self.s))
	def output(self):
		return {'s': self.s_a, 'psi_s': self.psi_s_a}
	def leak(self, s):
		"""Leakage (um/s) """    
		try:
			ans = .11574*self.KS*s**(2.*self.B + 3.)         
		except OverflowError:
			ans = 0.
		return ans                          
	def psi_s(self, s):
		"""Soil Potential (MPa)"""
		return self.PSI_SS*(s**-self.B)  
	def evap(self, s): 
		"""Soil evaporation rate, per unit ground area (mm/day)"""
		if s > self.SH:
			return self.EVMAX*(s - self.SH)/(1. - self.SH)
		else:
			return 0.

class DrydownSoil(object):
	def __init__(self):
		pass
	def snew(self, soil, dt, zr, qs):
		return (dt/(soil.N*zr*10.**6)*(-qs - (soil.evap(soil.s)*1000.)/(24.*60*60)- soil.leak(soil.s))) + soil.s

class ConstantSoil(object):
	def __init__(self):
		pass
	def snew(self, soil, dt, zr, qs):
		return soil.s

class StochasticSoil(object):
	"""takes alpha in cm, lda in 1/d"""
	def __init__(self, alpha, lda):
		self.alpha = alpha # Average rainfall depth (cm)
		self.lambda_r = lda # Average rainfall interarrival rate (days-1)
	def rain(self, dt, gamma):
		if np.random.random() > self.lambda_r*dt/(3600.*24.):
			return 0.
		else:
			return np.random.exponential(1./gamma)
	def sLoss(self, soil, dt, zr, qs):
		return (dt/(soil.N*zr*10.**6)*(-qs - (soil.evap(soil.s)*1000.)/(24.*60*60)- soil.leak(soil.s))) + soil.s
	def snew(self, soil, dt, zr, qs):
		gamma = (soil.N*zr*100.)/self.alpha; #Normalized Depth of Rainfall
		return min(1., self.sLoss(soil, dt, zr, qs) + self.rain(dt, gamma))

class RainSoil(object):
	def __init__(self):
		pass
	def sLoss(self, soil, dt, zr, qs):
		return (dt/(soil.N*zr*10.**6)*(-qs - (soil.evap(soil.s)*1000.)/(24.*60*60)- soil.leak(soil.s))) + soil.s
	def snew(self, soil, dt, zr, qs):
		"""Takes input of rainfall in mm"""
		return min(1., self.sLoss(soil, dt, zr, qs) + soil.rain_amt/(soil.N*zr*1000.))

class SetSoil(object):
	def __init__(self):
		pass
	def sLoss(self, soil, dt, zr, qs):
		return (dt/(soil.N*zr*10.**6)*(-qs - (soil.evap(soil.s)*1000.)/(24.*60*60)- soil.leak(soil.s))) + soil.s
	def snew(self, soil, dt, zr, qs):
		return soil.sm_inp

class SaltySoil(Soil):
	TS = 293. # Soil water temp (K)
	IV = 2. # Van't hoff coefficient for NaCl
	E = 0.95 # Plant filtration efficiency (-)
	def __init__(self, stype, zr, s, cs):
		Soil.__init__(self, stype, zr, s)
		self.cs = cs # Salt concentration in soil (mol/m3)
		self.MS = cs*self.ZR*self.N*s # mass of salt in soil (mol/m2)
		self.cs_a = []
	def update(self, dt, zr, qs):
		self.s = (dt/(self.N*zr*10.**6)*(-qs - (self.evap(self.s)*1000.)/(24.*60*60)- self.leak(self.s))) + self.s
		self.cs = self.MS/(self.s*self.N*self.ZR) # salt concentration in soil (mol/m3)
		self.s_a.append(self.s)
		self.cs_a.append(self.cs)
	def output(self):
		return {'s': self.s_a, 'cs': self.cs_a}
	def psi_s(self, s):
		return self.PSI_SS*(s**-self.B) - E*self.cs*R*self.IV*self.TS*10.**(-6.)

class Loam(object):
	PSI_SS = -1.43*10.**-3. # Saturation suction (MPa)
	B = 5.39 # Exponent of soil retention curve (-)
	KS = 20. # Saturated hydraulic conductivity (cm/d)
	N = .45 # Soil porosity (-)
	SH = .19 # Hygroscopic point of soil (-)
	def __init__(self):
		pass
		
class Sand(object):
	PSI_SS = -.34*10**-3 # Saturation suction (MPa)
	B = 4.05 # Exponent of soil retention curve (-)
	KS = 200. # Saturated hydraulic conductivity (cm/d)
	N = .35 # Soil porosity (-)
	SH = .08 # Hygroscopic point of soil (-)
	def __init__(self):
		pass

class SandyLoam(object):
	PSI_SS = -.7*10**-3 # Saturation suction (MPa)
	B = 4.9 # Exponent of soil retention curve (-)
	KS = 80. # Saturated hydraulic conductivity (cm/d)
	N = .43 # Soil porosity (-)
	SH = .14 # Hygroscopic point of soil (-)
	def __init__(self):
		pass

class LoamySand(object):
	PSI_SS = -.17*10**-3 # Saturation suction (MPa)
	B = 4.38 # Exponent of soil retention curve (-)
	KS = 100. # Saturated hydraulic conductivity (cm/d)
	N = .42 # Soil porosity (-)
	SH = .08 # Hygroscopic point of soil (-)
	def __init__(self):
		pass

class Clay(object):
	PSI_SS = -1.82*10**-3 # Saturation suction (MPa)
	B = 11.4 # Exponent of soil retention curve (-)
	KS = 1. # Saturated hydraulic conductivity (cm/d)
	N = .5 # Soil porosity (-)
	SH = .47 # Hygroscopic point of soil (-)
	def __init__(self):
		pass
