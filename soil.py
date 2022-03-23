from math import exp, pi, sqrt, log
from scipy.optimize import fsolve
from sympy import *
import numpy as np
from dics import *
from functions import *

class Soil(object):
	EVMAX = 3.
	def __init__(self, stype, dynamics, zr, s):
		self.PSI_SS = stype.PSI_SS
		self.B = stype.B
		self.KS = stype.KS
		self.N = stype.N
		self.SH = stype.SH
		self.ZR = zr
		self.s = s
		self.s_a = []
		self.psi_s_a = []
		self.dynamics = dynamics
		self.rain_amt = 0 # this rain amount is in mm!!
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
		self.alpha = alpha
		self.lambda_r = lda
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
	TS = 293. # soil water temp (K)
	IV = 2. # van't hoff coefficient for NaCl
	E = 0.95
	def __init__(self, stype, zr, s, cs):
		Soil.__init__(self, stype, zr, s)
		self.cs = cs # salt concentration in soil, mol/m3
		self.MS = cs*self.ZR*self.N*s # mass of salt in soil, mol/m2
		self.cs_a = []
	def update(self, dt, zr, qs):
		self.s = (dt/(self.N*zr*10.**6)*(-qs - (self.evap(self.s)*1000.)/(24.*60*60)- self.leak(self.s))) + self.s
		self.cs = self.MS/(self.s*self.N*self.ZR) # salt concentration in soil, mol/m3
		self.s_a.append(self.s)
		self.cs_a.append(self.cs)
	def output(self):
		return {'s': self.s_a, 'cs': self.cs_a}
	def psi_s(self, s):
		return self.PSI_SS*(s**-self.B) - E*self.cs*R*self.IV*self.TS*10.**(-6.)

class Loam(object):
	PSI_SS = -1.43*10.**-3.
	B = 5.39
	KS = 20.
	N = .45
	SH = .19
	def __init__(self):
		pass

class Sand(object):
	PSI_SS = -.34*10**-3
	B = 4.05
	KS = 200.
	N = .35
	SH = .08
	def __init__(self):
		pass

class SandyLoam(object):
	PSI_SS = -.7*10**-3
	B = 4.9
	KS = 80.
	N = .43
	SH = .14
	def __init__(self):
		pass

class LoamySand(object):
	PSI_SS = -.17*10**-3
	B = 4.38
	KS = 100.
	N = .42
	SH = .08
	def __init__(self):
		pass

class Clay(object):
	PSI_SS = -1.82*10**-3
	B = 11.4
	KS = 1.
	N = .5
	SH = .47
	def __init__(self):
		pass
