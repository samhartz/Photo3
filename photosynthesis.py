from math import exp, pi, sqrt, log
from scipy.optimize import fsolve
from sympy import *
import numpy as np
from dics import *
from functions import *

class Photo(object):

	TO = 293.2 # Reference Temperature for photosynthetic parameters (K)
	GAMMA_1 = .0451 # Parameter for temp dependence of CO2 compensation point (1/K)
	GAMMA_2 = .000347 # Parameter for temp dependence of CO2 compensation point (1/K^2)
	KC0 = 302. # Michaelis constant for C02 at TO (umol/mol)
	KO0 = 256. # Michaelis constant for 02 at TO (mmol/mol)
	OI = .209  # Oxygen Concentration (mol/mol)
	SVC = 649. # Entropy term for carboxylation (J/mol)
	SVQ = 646. # Entropy term for e-transport (J/mol)
	HKC =  59430. # Activation Energy for Kc (J/mol)
	HKO =  36000. # Activation Energy for Ko (J/mol)
	HKR =  53000. # Activation Energy for Rd (J/mol)
	HDJ = 200000. # Deactivation Energy for Jmax (J/mol)
	HAJ = 50000. # Activation Energy for Jmax (J/mol)
	RD0 = .32 # Standard Dark respiration at 25 C (umol/(m^2s))
	HAV =  72000.  # Activation Energy for Vc,max (J/mol)
	HDV =  200000. # Deactivation Energy for Vc,max (J/mol)

	def __init__(self, ptype, species):
		self.ptype = ptype
		if hasattr(species, 'RD0'): self.RD0 = species.RD0
		if hasattr(species, 'HAV'): self.HAV = species.HAV
		if hasattr(species, 'HDV'): self.HDV = species.HDV
		self.VCMAX0 = species.VCMAX0
		self.JMAX0 = species.JMAX0
		self.PSILA0 = species.PSILA0
		self.PSILA1 = species.PSILA1
		self.ared = 1.
		self.light_atten = 1.

	def a_c(self, ci, tl, ared):
		"""Rubisco-limited photosynthetic rate (umol/(m^2s^1))"""
		return self.v_cmax(tl, ared)*(ci - self.gamma(tl))/(ci + self.k_c(tl)*(1. + (self.OI*1000.)/self.k_o(tl)))
	def v_cmax(self, tl, ared):
		"""Maximum carboxylation rate (umol/(m^2s))"""
		return ared*self.VCMAX0*exp(self.HAV/(R*self.TO)*(1. - self.TO/tl))/(1. + exp((self.SVC*tl - self.HDV)/(R*tl)))
	def k_o(self, tl):
		"""Michaelis-menten coefficient for O2"""
		return self.KO0*exp(self.HKO/(R*self.TO)*(1. - self.TO/tl))
	def k_c(self, tl):
		"""Michaelis-menten coefficient for CO2"""
		return self.KC0*exp(self.HKC/(R*self.TO)*(1. - self.TO/tl))
	def a_q(self, phi, ci, tl):
		"""Light-limited photosynthetic rate (umol/(m^2s^1))"""
		return (self.j(phi*self.light_atten, tl)*(ci - self.gamma(tl)))/(4.*(ci + 2.*self.gamma(tl)))
	def gamma(self, tl):
		"""CO2 compensation point (umol/mol)"""
		return self.GAMMA_0*(1. + self.GAMMA_1*(tl - self.TO) + self.GAMMA_2*(tl - self.TO)**2.);
	def jmax(self, tl):
		"""Max. e- transport rate (umol/(m^2s))"""
		return self.JMAX0*exp(self.HAJ/(R*self.TO)*(1. - self.TO/tl))/(1. + exp((self.SVQ*tl - self.HDJ)/(R*tl))) 
	def j(self, phi, tl):
		"""Electron transport rate (umol/(m^2s))"""
		return min((phi*10.**6)/(EP*NA)*self.KAPPA_2*.5, self.jmax(tl)) 
	def jpar(self, phi, tl):
		"""Electron transport rate (umol/(m^2s), based off of PAR, not total solar radiatoion)"""
		return min(phi*self.KAPPA_2, self.jmax(tl)) 
	def a_phiciTl(self, phi, ci, tl, ared):
		"""Net photosynthetic demand for CO2 (umol/(m^2s^1))"""
		return max(min(self.a_c(ci, tl, 1.), self.a_q(phi, ci, tl)),0)*ared
	def a_psilc02(self, psi_l):  
		"""Vulnerability curve for water potential (-)"""
		if psi_l < self.PSILA0:
			return 0.
		elif self.PSILA0 <= psi_l <= self.PSILA1 :
			return (psi_l - self.PSILA0)/(self.PSILA1  - self.PSILA0)
		else: 
			return 1.
	def r_d(self, tl):
		"""Dark respiration flux (umol/(m^2s))"""
		return self.RD0*exp(self.HKR/(R*self.TO)*(1. - self.TO/tl))
	def csNew(self, an):
		"""CO2 concentration at leaf surface (ppm)"""
		return ca - an/self.GA
	def ciNew(self, cs, ta, qa):
		"""CO2 concentration in mesophyll cytosol (ppm)""" 
		# Dx = .0068 # .0068 #kg/kg
		# Drh = esat(ta)*.622/101325 - qa
		# ca = 350.
		# a1temp = 6.4 # c3 value

		# return cs*(1 - (1+Drh/Dx)/a1temp)

		return cs*(1.-1./(self.A1*self.fD(VPD(ta, qa))))  
	def cmNew(self, cs, ta, qa):
		"""CO2 concentration in mesophyll (ppm)"""
		return self.ciNew(cs, ta, qa) 
	def fD(self, vpd):
		"""Stomatal response to vapor pressure deficit (-)"""

		if vpd < 0.1:
			return 1.
		else:
			return 3/13./sqrt(vpd/1000.)
	def gsc(self, phi, ta, psi_l, qa, tl, cx, ared, **kwargs):
		"""Stomatal conductance to CO2, per unit leaf area (mol/m2/s)"""
		if self.an(phi, psi_l, tl, cx, ared, **kwargs) < 0.:
			return 0.
		else:
			return self.A1*self.an(phi, psi_l, tl, cx, ared, **kwargs)/self.ca*self.fD(VPD(ta, qa))

	# def gsc(self, phi, ta, psi_l, qa, tl, cx, ared, **kwargs):
	#     """Stomatal conductance to CO2, per unit leaf area (mol/m2/s)"""
	#     Dx = .0068 # .0068 #kg/kg
	#     Drh = esat(ta)*.622/101325 - qa
	#     ca = 350.
	#     a1temp = 6.4 # c3 value
	#     if self.an(phi, psi_l, tl, cx, ared, **kwargs) < 0.:
	#         return 0.
	#     else:
	#     	return a1temp*self.an(phi, psi_l, tl, cx, ared, **kwargs)/(ca*(1+Drh/Dx))
	#         # return a1temp*self.an(phi, psi_l, tl, cx, ared, **kwargs)/(self.cs*(1+Drh/Dx))


class C3(Photo):
	A1 = 15.
	GAMMA_0 = 34.6
	RC = 0.7
	GMGSRATIO = 1.65

	KAPPA_2 = .3 # Quantum yield of photosynthesis (mol CO2/mol photon) 
	def __init__(self, species, atm):
		Photo.__init__(self, "C3", species)
		self.ca = atm.ca
		self.cs = atm.ca
		self.ci = self.ciNew(atm.cs, atm.ta, atm.qa)
		self.cm = self.cmNew(atm.cs, atm.ta, atm.qa)
		self.cx = self.cm
		self.a_a = []
	def update(self, atm, psi_l, tl, dt):
		self.ci = self.ciNew(self.cs, atm.ta, atm.qa)
		self.cm = self.cmNew(self.cs, atm.ta, atm.qa)
		self.a = self.an(atm.phi, psi_l, tl, self.cm, self.ared)
		self.a_a.append(self.a)
	def output(self):
		return {'a': self.a_a}
	def an(self, phi, psi_l, tl, ci, ared): 
		"""Photosynthetic rate, per unit leaf area (umol/(m^2s))"""
		return self.a_psilc02(psi_l)*self.a_phiciTl(phi, ci, tl, ared)

class C4(Photo):
	A1 = 0.5*15.
	GAMMA_0 = 10.
	RC = 0.4
	GMGSRATIO = 2.65
	GBS = .013 # Conductance between bundle sheath and mesophyll (mol m^-2s^-1)
	VPMAX0 = 120. # Maximum PEP carboxylase under well-watered conditions (umol/(m^2s))
	VPR = 80. # PEP regeneration rate (umol/(m^2s))
	KP = 80. # Michaelis-Menten coefficient for C4 species (umol/(mol))

	KAPPA_2 = .3 # Quantum yield of photosynthesis (mol CO2/mol photon) 
	def __init__(self, species, atm):
		Photo.__init__(self, "C4", species)
		self.ca = atm.ca
		self.cs = atm.ca
		self.ci = self.ciNew(atm.cs, atm.ta, atm.qa)
		self.cm = self.cmNew(atm.cs, atm.ta, atm.qa)
		self.cbs = self.cm
		self.cx = self.cbs
		self.a_a = []

	def update(self, atm, psi_l, tl, dt):
		self.ci = self.ciNew(self.cs, atm.ta, atm.qa)
		self.cm = self.cmNew(self.cs, atm.ta, atm.qa)
		self.cbs = self.cbsNew(self.an(atm.phi, psi_l, tl, self.cbs, self.ared), self.cm, psi_l)
		self.cx = self.cbs
		self.a = self.an(atm.phi, psi_l, tl, self.cbs, self.ared)
		self.a_a.append(self.a)

	def output(self):
		return {'a': self.a_a}

	def an(self, phi, psi_l, tl, ci, ared): 
		"""Photosynthetic rate, per unit leaf area (umol/(m^2s))"""
		return self.a_psilc02(psi_l)*self.a_phiciTl(phi, self.cbs, tl, ared)
	def v_p(self, psi_l, ci):
		"""CO2 concentrating flux (umol/m2/s)"""
		return min((ci*self.VPMAX0)/(ci + self.KP), self.VPR)
	def cbsNew(self, an, cm, psi_l):
		"""CO2 concentration in bundle sheath cell (ppm)"""
		return (self.v_p(psi_l, cm) - an)/self.GBS + cm

class CAM(Photo):
	A1 = .8*15. # 0.6*15. # this is the value consistent with the Leuning model
	GAMMA_0 = 34.6
	RC = 0.5
	GMGSRATIO = 1.
	TR = 90.; # Relaxation time for circadian oscillator (min)
	C0 = 3000. # parameter for decarboxylation of malic acid (umol/mol)
	ALPHA_1 = 1/100.
	ALPHA_2 = 1/7. 
	K = .003 
	TOPT = 288.65 # (K)
	VCM = 0.0027 # Value controlling relative storage of malate (m)
	MU = .5 # Circadian oscillator constant
	BETA = 2.764 # Circadian oscillator constant
	CIRC_1 = .365 # Circadian oscillator constant
	CIRC_2 = .55 # Circadian oscillator constant
	CIRC_3 = 10. # Circadian oscillator constant
	Z0 = .55 # Initial value of z (-)
	M0 = 0. # 1000. # 0. Initial Malic Acid Carbon Concentration (umol/m^3)
	TH = 302.65 # 302.65 High temperature for CAM model (K)
	TW = 283.15 # 283.15 Low temperature for CAM model (K)

	KAPPA_2 = .1 # Quantum yield of photosynthesis (mol CO2/mol photon) (note that this overrides the value of 0.3 for typical photosynthesis)
	# KAPPA_2 = .3 # Quantum yield of photosynthesis (mol CO2/mol photon) (note that this overrides the value of 0.3 for typical photosynthesis)

	def __init__(self, species, atm):
		Photo.__init__(self, "CAM", species)
		self.MMAX = species.MMAX
		self.AMMAX = species.AMMAX
		self.ca = atm.ca
		self.cs = atm.ca
		self.ci = self.ciNew(atm.cs, atm.ta, atm.qa)
		self.cm = self.cmNew(atm.cs, atm.ta, atm.qa)
		self.z = self.Z0
		self.m = self.M0
		self.cc = self.ccNew(self.cs, atm.ta, atm.qa, self.z, self.m)
		self.cx = self.cc
		self.a_a = []
		self.z_a = []
		self.m_a = []
		self.asc_a = []
		self.asv_a = []
		self.avc_a = []
		self.aphicitl_a = []
		self.aphicctl_a = []
		self.ci_a = []
		self.cc_a = []
		self.cs_a = []
		
	def update(self, atm, psi_l, tl, dt):
		self.ci = self.ciNew(self.cs, atm.ta, atm.qa)
		self.cm = self.cmNew(self.cs, atm.ta, atm.qa)
		self.cc = self.ccNew(self.cs, atm.ta, atm.qa, self.z, self.m)
		self.cx = self.cc
		self.z = self.zNew(atm.phi, self.m, self.z, tl, dt) 
		self.m = self.mNew(atm.phi, psi_l, self.cc, tl, self.z, self.m, self.ared, dt)
		# self.a = self.an(atm.phi, psi_l, tl, self.cc, self.ared)
		self.a = self.an(atm.phi, psi_l, tl, self.ci, self.ared)
		self.z_a.append(self.z)
		self.m_a.append(self.m)
		self.a_a.append(self.a)
		self.asc_a.append(self.a_sc(atm.phi, psi_l, tl, self.ci, self.z, self.m, self.ared))
		self.aphicitl_a.append(self.a_phiciTl(atm.phi, self.ci, tl, self.ared))
		self.aphicctl_a.append(self.a_phiciTl(atm.phi, self.cc, tl, self.ared))
		self.ci_a.append(self.ci)
		self.cc_a.append(self.cc)

	def output(self):
		return {'a': self.a_a, 'z': self.z_a, 'm': self.m_a, 'asc': self.asc_a, 'aphicitl': self.aphicitl_a, 'aphicctl': self.aphicctl_a, 'ci': self.ci_a, 'cc': self.cc_a}

	def a_sc(self, phi, psi_l, tl, ci, z, m, ared):
		"""Flux from stomata to Calvin cycle (umol/(m^2s))"""
		## return self.a_psilc02(psi_l)*(self.a_phiciTl(phi, ci, tl, ared) - self.r_dc(phi, tl))*(1. - self.f_c(z, m))
		return max(0, self.a_psilc02(psi_l)*(self.a_phiciTl(phi, ci, tl, ared) - self.r_dc(phi, tl))*(1. - self.f_c(z, m)))
	def r_dv(self, phi, tl):
		"""Flux of dark respiration to vacuole (umol/(m^2s))"""
		return self.r_d(tl)*exp(-phi)
	def r_dc(self, phi, tl):
		"""Flux of dark respiration to calvin cycle (umol/(m^2s))"""
		return self.r_d(tl)*(1. - exp(-phi))
	def f_o(self, z):
		"""Circadian order function (-)"""
		return exp(-(z/self.MU)**self.CIRC_3)
	def f_m(self, z, m, tl):
		"""Malic acid storage function"""
		return self.f_o(z)*(self.m_s(tl) - m)/(self.ALPHA_2*self.m_s(tl) + (self.m_s(tl) - m))
	def m_s(self, tl):
		return self.MMAX*((self.TH - tl)/(self.TH - self.TW)*(1. - self.ALPHA_2) + self.ALPHA_2)
	def f_c(self, z, m):
		"""Carbon circadian control function"""
		return (1. - self.f_o(z))*m/(self.ALPHA_1*self.MMAX + m)
	def a_sv(self, phi, tl, psi_l, z, m):
		"""Flux from stomata to vacuole (umol/(m^2s))"""

		## if self.MMAX*((self.TH - tl)/(self.TH - self.TW)*(1. - self.ALPHA_2) + self.ALPHA_2) > m:
		## 	return (self.AMMAX*(1. - self.K*(tl - self.TOPT)**2.) - self.r_dv(phi, tl))*self.f_m(z, m, tl)*self.a_psilc02(psi_l)
		## else:
		## 	return 0.

		if self.MMAX*((self.TH - tl)/(self.TH - self.TW)*(1. - self.ALPHA_2) + self.ALPHA_2) > m and (1. - self.K*(tl - self.TOPT)**2.) >0:
			return (self.AMMAX*(1. - self.K*(tl - self.TOPT)**2.) - self.r_dv(phi, tl))*self.f_m(z, m, tl)*self.a_psilc02(psi_l)
		else:
			return 0.
	def a_vc(self, phi, cc, tl, z, m, ared):
		"""Flux from vacuole to calvin cycle (umol/(m^2s))"""
		return (self.a_phiciTl(phi, cc, tl, ared) - self.r_dc(phi, tl))*self.f_c(z, m)
	def m_e(self, z, m, tl, phi): 
		"""Malic acid equilibrium value"""
		# return self.MMAX*(self.CIRC_1*((self.TH - tl)/(self.TH - self.TW) + 1.)*(self.BETA*(z - self.MU))**3. - self.BETA*(self.TH - tl)/(self.TH - self.TW)*(z - self.MU) + \
		#         self.CIRC_2*(self.TH - tl)/(self.TH - self.TW)) #-(1- self.f_o(z))*(1-m/(m+self.ALPHA_1*self.MMAX))) 


		if phi>0.:
			# return self.MMAX*(self.CIRC_1*((self.TH - tl)/(self.TH - self.TW) + 1.)*(self.BETA*(z - self.MU))**3. - self.BETA*(self.TH - tl)/(self.TH - self.TW)*(z - self.MU) + \
		 #        self.CIRC_2*(self.TH - tl)/(self.TH - self.TW) )
			return self.MMAX*(self.CIRC_1*((self.TH - tl)/(self.TH - self.TW) + 1.)*(self.BETA*(z - self.MU))**3. - self.BETA*(self.TH - tl)/(self.TH - self.TW)*(z - self.MU) + \
			   self.CIRC_2*(self.TH - tl)/(self.TH - self.TW) -(1- self.f_o(z))*(1-m/(m+self.ALPHA_1*self.MMAX))) 
			
		else:
			return self.MMAX*(self.CIRC_1*((self.TH - tl)/(self.TH - self.TW) + 1.)*(self.BETA*(z - self.MU))**3. - self.BETA*(self.TH - tl)/(self.TH - self.TW)*(z - self.MU) + \
				self.CIRC_2*(self.TH - tl)/(self.TH - self.TW)+ (1-self.f_o(z))) 

	def zNew(self, phi, m, z, tl, dt):
		return max(0, dt*(m - self.m_e(z, m, tl, phi))/(self.MMAX*60.*self.TR) + z)
	def an(self, phi, psi_l, tl, ci, ared): 
		"""Photosynthetic rate, per unit leaf area (umol/(m^2s))"""
		return self.a_sc(phi, psi_l, tl, ci, self.z, self.m, ared) + self.a_sv(phi, tl, psi_l, self.z, self.m) 
	def ccNew(self, cs, ta, qa, z, m):
		"""CO2 concentration in mesophyll cytosol resulting from malic acid decarboxylation (ppm)"""
		return self.cmNew(cs, ta, qa) + self.f_c(z, m)*self.C0
	def mNew(self, phi, psi_l, cc, tl, z, m, ared, dt): 
		"""Malic acid concentration"""
		return max(((dt/ self.VCM)*(self.a_sv(phi, tl, psi_l, z, m) - self.a_vc(phi, cc, tl, z, m, ared) + self.r_dv(phi, tl))) + m, 0.)
