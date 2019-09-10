from math import exp, pi, sqrt, log
from scipy.optimize import fsolve
from sympy import *
import numpy as np
from dics import *
from functions import *

class Simulation(object):
	def __init__(self, species_cls, atm_cls, soil_cls, photo_cls, hydro_cls):
		self.species = species_cls
		self.atm = atm_cls
		self.soil = soil_cls
		self.photo = photo_cls
		self.hydro = hydro_cls
	def update(self, dt, phi, ta, qa):
		self.atm.update(phi, ta, qa) 
		self.photo.update(self.atm, self.hydro.psi_l, self.hydro.tl, dt)
		self.hydro.update(self.atm, self.soil, self.photo, dt)
		self.soil.update(dt, self.species.ZR, self.hydro.qs)
	def output(self):
		out = {}
		out.update(self.photo.output())
		out.update(self.hydro.output())
		out.update(self.soil.output())
		return out

class Soil(object):
	EVMAX = 1.
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
	def snew():
		return min(1., self.sLoss(soil, dt, zr, qs) + self.rain(dt, gamma))

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


class Atmosphere(object):
	ca = 400.
	def __init__(self, phi, ta, qa):
		self.phi = phi
		self.ta = ta
		self.qa = qa
		self.cs = self.ca
	def update(self, phi, ta, qa):
		self.phi = phi
		self.ta = ta
		self.qa = qa

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
	HAV =  72000. # Activation Energy for Vc,max (J/mol)
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
	    return max(min(self.a_c(ci, tl, ared), self.a_q(phi, ci, tl)),0)
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
	A1 = 0.6*15.
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
	Z0 = .55  # Initial value of z (-)
	M0 = 0. # Initial Malic Acid Carbon Concentration (umol/m^3)
	TH = 302.65 # High temperature for CAM model (K)
	TW = 283.15 # Low temperature for CAM model (K)

	KAPPA_2 = .1 # Quantum yield of photosynthesis (mol CO2/mol photon) (note that this overrides the value of 0.3 for typical photosynthesis)

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
		
	def update(self, atm, psi_l, tl, dt):
		self.ci = self.ciNew(self.cs, atm.ta, atm.qa)
		self.cm = self.cmNew(self.cs, atm.ta, atm.qa)
		self.cc = self.ccNew(self.cs, atm.ta, atm.qa, self.z, self.m)
		self.cx = self.cc
		self.z = self.zNew(atm.phi, self.m, self.z, tl, dt) 
		self.m = self.mNew(atm.phi, psi_l, self.cc, tl, self.z, self.m, dt)
		self.a = self.an(atm.phi, psi_l, tl, self.cc, self.ared)
		self.z_a.append(self.z)
		self.m_a.append(self.m)
		self.a_a.append(self.a)

	def output(self):
		return {'a': self.a_a, 'z': self.z_a, 'm': self.m_a}

	def a_sc(self, phi, psi_l, tl, ci, z, m, ared):
	    """Flux from stomata to Calvin cycle (umol/(m^2s))"""
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
	    return self.f_o(z)*(self.MMAX*((self.TH - tl)/(self.TH - self.TW)*(1. - self.ALPHA_2) + self.ALPHA_2) - m)/(self.ALPHA_2*self.MMAX*((self.TH - tl)/\
	        (self.TH - self.TW)*(1. - self.ALPHA_2) + self.ALPHA_2) + (self.MMAX*((self.TH - tl)/(self.TH - self.TW)*(1. - self.ALPHA_2) + self.ALPHA_2) - m))
	def f_c(self, z, m):
	    """Carbon circadian control function"""
	    return (1. - self.f_o(z))*m/(self.ALPHA_1*self.MMAX + m)
	def a_sv(self, phi, tl, psi_l, z, m):
		"""Flux from stomata to vacuole (umol/(m^2s))"""
		if self.MMAX*((self.TH - tl)/(self.TH - self.TW)*(1. - self.ALPHA_2) + self.ALPHA_2) > m and (1. - self.K*(tl - self.TOPT)**2.) >0:
			return (self.AMMAX*(1. - self.K*(tl - self.TOPT)**2.) - self.r_dv(phi, tl))*self.f_m(z, m, tl)*self.a_psilc02(psi_l)
		else:
			return 0.
	def a_vc(self, phi, cc, tl, z, m):
	    """Flux from vacuole to calvin cycle (umol/(m^2s))"""
	    return (self.a_phiciTl(phi, cc, tl, 1.) - self.r_dc(phi, tl))*self.f_c(z, m)
	def m_e(self, z, m, tl, phi): 
	    """Malic acid equilibrium value"""
	    if phi>0.:
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
	def mNew(self, phi, psi_l, cc, tl, z, m, dt): 
		"""Malic acid concentration"""
		return max(((dt/ self.VCM)*(self.a_sv(phi, tl, psi_l, z, m) - self.a_vc(phi, cc, tl, z, m) + self.r_dv(phi, tl))) + m, 0.)


class Hydro(object):
	A_ROOT = 8.
	def __init__(self, species):
		self.GPMAX = species.GPMAX
		self.GA = species.GA
		self.gp = species.GPMAX
		self.GCUT = species.GCUT
		self.RAIW = species.RAIW
		self.zr = species.ZR
		self.lai = species.LAI
		self.psi_l_a = []
		self.gp_a = []
		self.gsv_a = []
		self.tl_a = []
		self.ev_a = []
	def rai(self, s):
		"""Root area index (-)"""
		return self.RAIW*s**-self.A_ROOT
	def evf(self, photo, phi, ta, psi_l, qa, tl, ci, lai, ared, **kwargs):
	    """Transpiration, per unit ground area (um/sec)"""
	    if self.gsw(photo, phi, ta, psi_l, qa, tl, ci, ared, **kwargs) < 0.00001:
	    	return 0.
	    else:
	    	return max(lai*(1./(self.gsw(photo, phi, ta, psi_l, qa, tl, ci, ared, **kwargs)*R*ta/P_ATM*1000000.)+1./(self.GA*1000.))**(-1.)\
	    	*RHO_A/RHO_W*(self.qi(tl, psi_l)-qa), 0.)
	def qi(self, tl, psi_l):
	    """Specific humidity internal to leaf (kg/kg)"""
	    try: 
	        ans =  .622*esat(tl)/P_ATM*exp(psi_l*1000000.*VW/R/tl)
	    except OverflowError:
	        ans = 0.
	    return ans

	    #return .622*esat(tl)/P_ATM*exp(psi_l*1000000.*VW/R/tl)
	def gpf(self, psi_l):
	    """Plant conductance, per unit leaf area (um/(s-MPa))"""
	    if psi_l<-10:
	        return 0.
	    else:
	        return self.GPMAX*exp(-(-psi_l/2.)**2.)
	def shf(self, tl, ta, lai):
		"""Sensible heat flux (W/m^2), per unit ground area"""
		return CP_A*RHO_A*self.GA*(tl-ta)/1000.*lai
	def gsr(self, soil, s, zr):
	    """Soil-Root Conductance, per unit ground area (um/(s-MPa))"""
	    return (soil.leak(s)*sqrt(self.rai(s))*1000000.)/(float(pi)*g*RHO_W*zr)
	def gsw(self, photo, phi, ta, psi_l, qa, tl, ci, ared): 
	    """Stomatal conductance to water, per unit leaf area (mol/m2/sec)"""
	    #return gsN(phi, Ta, psi_l, qa, Tl, ci, t)*(1.6*(1. + GMGSRATIO[pType[species]]))/(1.6 + GMGSRATIO[pType[species]]) + (gcut[species]*po/(1000.*R*Ta))
	    return photo.gsc(phi, ta, psi_l, qa, tl, ci, ared)*1.6 + (self.GCUT*P_ATM/(1000.*R*ta))

class HydroNC(Hydro):
	def __init__(self, species, atm, soil, photo, vwi):
		Hydro.__init__(self, species)
		self.psi_l, self.tl = fsolve(self.fBal, (-1., 290.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr))
		self.gp = self.gpf(self.psi_l)
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, photo.ared)
	def update(self, atm, soil, photo, dt):
		self.psi_l, self.tl = fsolve(self.fBal, (-1., 290.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr))
		self.gp = self.gpf(self.psi_l)
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, photo.ared)
		self.qs = self.ev
		self.psi_l_a.append(self.psi_l)
		self.gp_a.append(self.gp)
		self.gsv_a.append(self.gsw(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, photo.ared))
		self.tl_a.append(self.tl) 
		self.ev_a.append(self.ev)
	def output(self):
		return {'psi_l': self.psi_l_a, 'gp': self.gp_a, 'gsv': self.gsv_a, 'tl': self.tl_a, 'ev': self.ev_a}
	def gsrp(self, soil, s, gp, lai, zr):
	    """Soil-Root-Plant Conductance, per unit ground area (um/(s-MPa))"""
	    return (lai*self.gsr(soil, s, zr)*gp)/(self.gsr(soil, s, zr) + lai*gp)
	def qsf(self, photo, phi, ta, psi_l, qa, tl, ci, lai):
		return self.evf(photo, phi, ta, psi_l, qa, tl, ci, lai, photo.ared)
	def fBal(self, p, soil, photo, phi, ta, qa, c1, s, lai, gp, ared, zr):
	    psi_l, tl =p

	    if lai < 1.: # assumes only a portion of solar radiation is absorbed by crops
	        return (phi*lai - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)/1000000., \
	            self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared) - self.gsrp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l))
	    else:
	        return (phi - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)/1000000., \
	            self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared) - self.gsrp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l))


class HydroCap(Hydro):
	F_CAP = 0.5
	def __init__(self, species, atm, soil, photo, vwi):
		Hydro.__init__(self, species)
		self.GWMAX = species.GWMAX
		self.VWT = species.VWT
		self.CAP = species.CAP
		self.vw = vwi*self.VWT
		self.psi_l, self.tl = fsolve(self.fBal, (-1., 290.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr, self.vw))
		self.gp = self.gpf(self.psi_l)
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, photo.ared)
		self.vw_a = []

	def update(self, atm, soil, photo, dt):
		self.psi_l, self.tl = fsolve(self.fBal, (-1., 290.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr, self.vw))
		self.gp = self.gpf(self.psi_l)
		self.vw = self.vwf(self.vw, self.ev, self.gp, self.psi_l, self.lai, dt)
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, photo.ared)
		self.qs = self.qsf(self.vw, self.ev, self.gp, self.psi_l, self.lai, dt)
		self.psi_l_a.append(self.psi_l)
		self.gp_a.append(self.gp)
		self.gsv_a.append(self.gsw(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, photo.ared))
		self.tl_a.append(self.tl) 
		self.ev_a.append(self.ev)
		self.vw_a.append(self.vw)

	def output(self):
		return {'psi_l': self.psi_l_a, 'gp': self.gp_a, 'gsv': self.gsv_a, 'tl': self.tl_a, 'ev': self.ev_a, 'vw': self.vw_a}

	def psi_wf(self, vw): 
	    """Water potential of stored water (MPa)"""
	    return (1./self.CAP)*vw/self.VWT - (1./self.CAP)
	def vwf(self, vw, ev, gp, psi_l, lai, dt):
	    """Stored water volume, per unit leaf area (m3/m2)"""
	    return min(vw - self.gwf(self.psi_wf(vw))*(self.psi_wf(vw) - (ev*(1. - self.F_CAP))/(lai*gp) - psi_l)*dt/10.**6, self.VWT)
	def psi_xf(self, ev, gp, psi_l):
	    """Water potential at connection node x (MPa)"""
	    return ev*(1. - self.F_CAP)/(lai*gp) + psi_l
	def qwf(self, vw, ev, gp, psi_l, lai, dt):
	    """Stored water flux, per unit ground area"""
	    return (vw - self.vwf(vw, ev, gp, psi_l, lai, dt))*lai*10.**6/dt
	def qsf(self, vw, ev, gp, psi_l, lai, dt):
	    """Soil water flux, per unit ground area"""
	    return ev - self.qwf(vw, ev, gp, psi_l, lai, dt)
	def gwf(self, psi_w):
		"""Xylem-storage conductance, per unit leaf area (um/(MPa-s))"""
		return self.GWMAX*exp(-(-psi_w/2.)**2.)
	    #return GWMAX[species]*(vw/VWT[species])**4. 
	def gsrfp(self, soil, s, gp, lai, zr):
	    """Soil-root-plant fraction conductance, per unit ground area (um/(s-MPa))"""
	    return (lai*self.gsr(soil, s, zr)*gp/self.F_CAP)/(self.gsr(soil, s, zr) +  lai*gp/self.F_CAP)
	def fBal(self, params, soil, photo, phi, ta, qa, c1, s, lai, gp, ared, zr, vw):
	    psi_l, tl = params
	    psi_w = self.psi_wf(vw)
	    if lai < 1.: # assumes only a portion of solar radiation is absorbed by crops
		    return (phi*lai - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)/1000000.,  \
		    	self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)\
		    	-(self.gsrfp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l) + lai*self.gwf(psi_w)*(psi_w - psi_l))/ \
		    	(1. + (self.gsrfp(soil, s, gp, lai, zr)*(1. - self.F_CAP))/(lai*gp) + (self.gwf(psi_w)*(1. - self.F_CAP))/gp))
	    else:
	    	return (phi - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)/1000000., 
	    		self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)
				-(self.gsrfp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l) + lai*self.gwf(psi_w)*(psi_w - psi_l))/
				(1. + (self.gsrfp(soil, s, gp, lai, zr)*(1. - self.F_CAP))/(lai*gp) + (self.gwf(psi_w)*(1. - self.F_CAP))/gp))

class HydroLeafCap(Hydro):
	F_CAP = 0.5
	def __init__(self, species, atm, soil, photo, vwi, vli, lcap, vlt, gleaf):
		Hydro.__init__(self, species)
		self.GWMAX = species.GWMAX
		self.VWT = species.VWT
		self.CAP = species.CAP
		self.LCAP = lcap
		self.GLMAX = gleaf
		self.VLT = vlt
		self.vl = vli*self.VLT
		self.vw = vwi*self.VWT
		self.psi_l = self.psi_hf(self.vl)
		self.psi_w = self.psi_wf(self.vw)

		self.psi_2, self.tl = fsolve(self.fBal, (-1., 290.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr, self.vw, self.vl))
		self.psi_1 = self.psi_1f(self.psi_w, self.psi_2, self.gp, soil, self.lai)
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, photo.ared, self.psi_2)
		self.vw_a = []
		self.vl_a = []
		self.qs_a = []
		self.qw_a = []
		self.ql_a = []
		self.psi_1_a = []
		self.psi_2_a = []
		self.psi_w_a = []

	def update(self, atm, soil, photo, dt):
		# self.vw = self.vwf(self.vw, self.psi_1, dt)
		# self.vl = self.vlf(self.vl, self.psi_l, self.psi_2, dt)
		# self.psi_w = self.psi_wf(self.vw)
		# self.psi_l = self.psi_lf(self.vl)
		# self.psi_2, self.tl = fsolve(self.fBal, (-1., 290.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr, self.vw, self.vl))
		# self.psi_1 = self.psi_1f(self.psi_l, self.psi_w, self.psi_2, self.gp, soil)
		self.gp = self.gpf(self.psi_l)
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, photo.ared, self.psi_2)
		self.qs = self.qsf(self.vw, self.vl, self.ev, self.psi_l, self.psi_1, self.psi_2, self.lai, dt)
		self.psi_l_a.append(self.psi_l)
		self.gp_a.append(self.gp)
		self.gsv_a.append(self.gsw(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, photo.ared))
		self.tl_a.append(self.tl) 
		self.ev_a.append(self.ev)
		self.vw_a.append(self.vw)
		self.vl_a.append(self.vl)
		self.qs_a.append(self.qs)
		self.psi_1_a.append(self.psi_1)
		self.psi_2_a.append(self.psi_2)
		self.psi_w_a.append(self.psi_w)
		self.qw_a.append(self.qwf(self.vw, self.psi_1, self.lai, dt))
		self.ql_a.append(self.qlf(self.vl, self.psi_l, self.psi_2, self.lai, dt))
		self.vw = self.vwf(self.vw, self.psi_1, dt)
		self.vl = self.vlf(self.vl, self.psi_l, self.psi_2, dt)
		self.psi_w = self.psi_wf(self.vw)
		self.psi_l = self.psi_hf(self.vl)
		self.psi_2, self.tl = fsolve(self.fBal, (-1., 290.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr, self.vw, self.vl))
		self.psi_1 = self.psi_1f(self.psi_w, self.psi_2, self.gp, soil, self.lai)

	def output(self):
		return {'psi_l': self.psi_l_a, 'gp': self.gp_a, 'gsv': self.gsv_a, 'tl': self.tl_a, 'ev': self.ev_a, 'vw': self.vw_a, 'vl' :self.vl_a, 'qs':self.qs_a, 'qw':self.qw_a, 'qs':self.qs_a, 'psi_1':self.psi_1_a, 'psi_2':self.psi_2_a, 'psi_w':self.psi_w_a}

	def psi_wf(self, vw): 
	    """Water potential of stored water (MPa)"""
	    return (1./self.CAP)*vw/self.VWT - (1./self.CAP)
	def psi_hf(self, vl):
		"""Water potential of stored leaf water (MPa)"""
		return (1./self.LCAP)*vl/self.VLT - (1./self.LCAP)
	def gl_f(self, psi_l):
		"""Leaf hydraulic conductance (um/MPa/s)"""
		return self.GLMAX*exp(-(-psi_l/2.)**2.)
	def psi_1f(self, psi_w, psi_2, gp, soil, lai):
		"""Water potential at connection node 1 (lower node) (MPa)"""
		return (psi_2*self.gfp(gp)*lai + psi_w*self.gwf(psi_w)*lai + soil.psi_s(soil.s)*self.gsrfp(soil, soil.s, gp, self.lai, self.zr))/(self.gfp(gp)*lai + self.gwf(psi_w)*lai + self.gsrfp(soil, soil.s, gp, self.lai, self.zr))
	def vwf(self, vw, psi_1, dt):
	    """Stored water volume, per unit leaf area (m3/m2)"""
	    return max(min(vw - self.gwf(self.psi_wf(vw))*(self.psi_wf(vw)-psi_1)*dt/10.**6, self.VWT), 0.)
	def vlf(self, vl, psi_l, psi_2, dt):
		return max(min(vl - self.gl_f(psi_l)*(psi_l-psi_2)*dt/10.**6, self.VLT), 0.)
	def qwf(self, vw, psi_1, lai, dt):
	    """Stored water flux, per unit ground area"""
	    return (vw - self.vwf(vw, psi_1, dt))*lai*10.**6/dt
	def qlf(self, vl, psi_l, psi_2, lai, dt):
	    """Stored water flux, per unit ground area"""
	    return (vl - self.vlf(vl, psi_l, psi_2, dt))*lai*10.**6/dt
	def qsf(self, vw, vl, ev, psi_l, psi_1, psi_2, lai, dt):
	    """Soil water flux, per unit ground area"""
	    return ev - self.qwf(vw, psi_1, lai, dt) - self.qlf(vl, psi_l, psi_2, lai, dt)
	def gwf(self, psi_w):
		"""Xylem-storage conductance, per unit leaf area (um/(MPa-s))"""
		return self.GWMAX*exp(-(-psi_w/2.)**2.)
	    #return GWMAX[species]*(vw/VWT[species])**4. 
	def gsrfp(self, soil, s, gp, lai, zr):
	    """Soil-root-plant fraction conductance, per unit ground area (um/(s-MPa))"""
	    return (lai*self.gsr(soil, s, zr)*gp/self.F_CAP)/(self.gsr(soil, s, zr) +  lai*gp/self.F_CAP)
	def gfp(self, gp):
		return gp/(1. - self.F_CAP)
	def fBal(self, params, soil, photo, phi, ta, qa, c1, s, lai, gp, ared, zr, vw, vl):
	    psi_2, tl = params
	    psi_w = self.psi_wf(vw)
	    psi_l = self.psi_hf(vl)
	    return (phi - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared, psi_2)/1000000., 

		    self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared, psi_2)
		     -(psi_l-psi_2)*self.gl_f(psi_l)*lai - (self.psi_1f(psi_w, psi_2, gp, soil, lai)-psi_2)*self.gfp(gp)*lai)

	def evf(self, photo, phi, ta, psi_l, qa, tl, ci, lai, ared, psi_2):
	    """Transpiration, per unit ground area (um/sec)"""
	    if self.gsw(photo, phi, ta, psi_l, qa, tl, ci, ared) < 0.00001:
	    	return 0.
	    else:
	    	return max(lai*(1./(self.gsw(photo, phi, ta, psi_l, qa, tl, ci, ared)*R*ta/P_ATM*1000000.)+1./(self.GA*1000.))**(-1.)\
	    	*RHO_A/RHO_W*(self.qi(tl, psi_2)-qa), 0.)

class HydroLeafCapOld(Hydro):
	F_CAP = 0.5
	def __init__(self, species, atm, soil, photo, vwi, vli, lcap, vlt, gleaf):
		Hydro.__init__(self, species)
		self.GWMAX = species.GWMAX
		self.VWT = species.VWT
		self.CAP = species.CAP
		self.LCAP = lcap
		self.GLMAX = gleaf
		self.VLT = vlt
		self.vl = vli*self.VLT
		self.vw = vwi*self.VWT
		self.psi_h = self.psi_hf(self.vl)
		self.psi_w = self.psi_wf(self.vw)
		self.psi_l, self.tl = fsolve(self.fBal, (-1., 290.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr, self.vw, self.vl))
		self.psi_1 = self.psi_1f(self.psi_w, self.psi_l, self.gp, soil, self.lai)
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, photo.ared)
		self.vw_a = []
		self.vl_a = []
		self.qs_a = []
		self.qw_a = []
		self.psi_1_a = []
		self.psi_h_a = []
		self.psi_w_a = []

	def update(self, atm, soil, photo, dt):
		# self.vw = self.vwf(self.vw, self.psi_1, dt)
		# self.vl = self.vlf(self.vl, self.psi_l, self.psi_2, dt)
		# self.psi_w = self.psi_wf(self.vw)
		# self.psi_l = self.psi_lf(self.vl)
		# self.psi_2, self.tl = fsolve(self.fBal, (-1., 290.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr, self.vw, self.vl))
		# self.psi_1 = self.psi_1f(self.psi_l, self.psi_w, self.psi_2, self.gp, soil)
		self.gp = self.gpf(self.psi_l)
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, photo.ared)
		self.qs = self.qsf(self.vw, self.vl, self.ev, self.psi_l, self.psi_1, self.psi_h, self.lai, dt)
		self.psi_l_a.append(self.psi_l)
		self.gp_a.append(self.gp)
		self.gsv_a.append(self.gsw(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, photo.ared))
		self.tl_a.append(self.tl) 
		self.ev_a.append(self.ev)
		self.vw_a.append(self.vw)
		self.vl_a.append(self.vl)
		self.qs_a.append(self.qs)
		self.psi_1_a.append(self.psi_1)
		self.psi_h_a.append(self.psi_h)
		self.psi_w_a.append(self.psi_w)
		self.qw_a.append(self.qwf(self.vw, self.psi_1, self.lai, dt))
		self.vw = self.vwf(self.vw, self.psi_1, dt)
		self.vl = self.vlf(self.vl, self.psi_h, self.psi_l, dt)
		self.psi_w = self.psi_wf(self.vw)
		self.psi_h = self.psi_hf(self.vl)
		self.psi_l, self.tl = fsolve(self.fBal, (-1., 290.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr, self.vw, self.vl))
		self.psi_1 = self.psi_1f(self.psi_w, self.psi_l, self.gp, soil, self.lai)

	def output(self):
		return {'psi_l': self.psi_l_a, 'gp': self.gp_a, 'gsv': self.gsv_a, 'tl': self.tl_a, 'ev': self.ev_a, 'vw': self.vw_a, 'vl' :self.vl_a, 'qs':self.qs_a, 'qw':self.qw_a, 'psi_1':self.psi_1_a, 'psi_h':self.psi_h_a, 'psi_w':self.psi_w_a}

	def psi_wf(self, vw): 
	    """Water potential of stored water (MPa)"""
	    return (1./self.CAP)*vw/self.VWT - (1./self.CAP)
	def psi_hf(self, vl):
		"""Water potential of stored leaf water (MPa)"""
		return (1./self.LCAP)*vl/self.VLT - (1./self.LCAP)
	def gl_f(self, psi_l):
		"""Leaf hydraulic conductance (um/MPa/s)"""
		return self.GLMAX*exp(-(-psi_l/2.)**2.)
	def psi_1f(self, psi_w, psi_l, gp, soil, lai):
		"""Water potential at connection node 1 (lower node) (MPa)"""
		return (psi_l*self.gfp(gp)*lai + psi_w*self.gwf(psi_w)*lai + soil.psi_s(soil.s)*self.gsrfp(soil, soil.s, gp, self.lai, self.zr))/(self.gfp(gp)*lai + self.gwf(psi_w)*lai + self.gsrfp(soil, soil.s, gp, self.lai, self.zr))
	def vwf(self, vw, psi_1, dt):
	    """Stored water volume, per unit leaf area (m3/m2)"""
	    return max(min(vw - self.gwf(self.psi_wf(vw))*(self.psi_wf(vw)-psi_1)*dt/10.**6, self.VWT), 0.)
	def vlf(self, vl, psi_h, psi_l, dt):
		return max(min(vl - self.gl_f(psi_l)*(psi_h-psi_l)*dt/10.**6, self.VLT), 0.)
	def qwf(self, vw, psi_1, lai, dt):
	    """Stored water flux, per unit ground area"""
	    return (vw - self.vwf(vw, psi_1, dt))*lai*10.**6/dt
	def qlf(self, vl, psi_h, psi_l, lai, dt):
	    """Stored water flux, per unit ground area"""
	    return (vl - self.vlf(vl, psi_h, psi_l, dt))*lai*10.**6/dt
	def qsf(self, vw, vl, ev, psi_l, psi_1, psi_h, lai, dt):
	    """Soil water flux, per unit ground area"""
	    return ev - self.qwf(vw, psi_1, lai, dt) - self.qlf(vl, psi_h, psi_l, lai, dt)
	def gwf(self, psi_w):
		"""Xylem-storage conductance, per unit leaf area (um/(MPa-s))"""
		return self.GWMAX*exp(-(-psi_w/2.)**2.)
	    #return GWMAX[species]*(vw/VWT[species])**4. 
	def gsrfp(self, soil, s, gp, lai, zr):
	    """Soil-root-plant fraction conductance, per unit ground area (um/(s-MPa))"""
	    return (lai*self.gsr(soil, s, zr)*gp/self.F_CAP)/(self.gsr(soil, s, zr) +  lai*gp/self.F_CAP)
	def gfp(self, gp):
		return gp/(1. - self.F_CAP)
	def fBal(self, params, soil, photo, phi, ta, qa, c1, s, lai, gp, ared, zr, vw, vl):
	    psi_l, tl = params
	    psi_w = self.psi_wf(vw)
	    psi_h = self.psi_hf(vl)
	    return (phi - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)/1000000., 
		    self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)
		     -(psi_h-psi_l)*self.gl_f(psi_l)*lai - ((psi_l*self.gfp(gp)*lai + psi_w*self.gwf(psi_w)*lai + soil.psi_s(soil.s)*self.gsrfp(soil, soil.s, gp, self.lai, self.zr))\
		     /(self.gfp(gp)*lai+self.gwf(psi_w)*lai+self.gsrfp(soil, soil.s, gp, self.lai, self.zr))-psi_l)*self.gfp(gp)*lai)

class Halophyte(Hydro):
	F_CAP = 0.5
	CW = 150. # default salt concentration in plant, mol/m3
	E = 0.95 # filtration efficiency, unitless
	TS = 293. # soil water temp (K)
	IV = 2. # van't hoff coefficient for NaCl
	CS = 150. # salt concentration in soil, mol/m3
	D1 = .028 # parameter for turgor pressure, for agave
	D2 = 8. # parameter for turgor pressure, for agave

	def __init__(self, species, atm, soil, photo, vwi):
		Hydro.__init__(self, species)
		self.GWMAX = species.GWMAX
		self.VWT = species.VWT
		self.vw = vwi*self.VWT
		self.MW = self.CW*self.vw
		self.CAP = species.CAP
		self.cw = self.MW/self.vw
		self.psi_l, self.tl = fsolve(self.fBal, (-1., 290.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cm, soil.s, self.lai, self.gp, 1., self.zr, self.psi_wf(self.vw, self.cw)))
		self.gp = self.gpf(self.psi_l)
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, 1.)
		self.vw_a = []
		self.cw_a = []
		

	def update(self, atm, soil, photo, dt):
		self.psi_l, self.tl = fsolve(self.fBal, (-1., 290.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cm, soil.s, self.lai, self.gp, 1., self.zr, self.psi_wf(self.vw, self.cw)))
		self.gp = self.gpf(self.psi_l)
		self.vw = self.vwf(self.vw, self.ev, self.gp, self.psi_l, self.lai, self.cw, dt)
		self.cw = self.MW/self.vw
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, 1.)
		self.qs = self.qsf(self.vw, self.ev, self.gp, self.psi_l, self.lai, self.cw, dt)
		self.psi_l_a.append(self.psi_l)
		self.gp_a.append(self.gp)
		self.gsv_a.append(self.gsw(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, 1.))
		self.tl_a.append(self.tl) 
		self.ev_a.append(self.ev)
		self.vw_a.append(self.vw)
		self.cw_a.append(self.cw)

	def output(self):
		return {'psi_l': self.psi_l_a, 'gp': self.gp_a, 'gsv': self.gsv_a, 'tl': self.tl_a, 'ev': self.ev_a, 'vw': self.vw_a, 'cw': self.cw_a}

	def psi_wf(self, vw, cw): 
	    """Water potential of stored water (MPa)"""
	    TL = 293.
	    return (vw/self.VWT-self.D1)**self.D2 - cw*R*self.IV*TL*10.**(-6.)
	def qwf(self, vw, ev, gp, psi_l, lai, cw, dt):
	    """Stored water flux, per unit ground area"""
	    return (vw - self.vwf(vw, ev, gp, psi_l, lai, cw, dt))*lai*10.**6/dt
	def qsf(self, vw, ev, gp, psi_l, lai, cw, dt):
	    """Soil water flux, per unit ground area"""
	    return ev - self.qwf(vw, ev, gp, psi_l, lai, cw, dt)
	def gwf(self, psi_w):
	    """Xylem-storage conductance, per unit leaf area (um/(MPa-s))"""
	    return self.GWMAX*exp(-(-psi_w/2.)**2.)
	    #return GWMAX[species]*(vw/VWT[species])**4. 
	def gsrfp(self, soil, s, gp, lai, zr):
	    """Soil-root-plant fraction conductance, per unit ground area (um/(s-MPa))"""
	    return (lai*self.gsr(soil, s, zr)*gp/self.F_CAP)/(self.gsr(soil, s, zr) +  lai*gp/self.F_CAP)
	def vwf(self, vw, ev, gp, psi_l, lai, cw, dt):
	    """Stored water volume, per unit leaf area (m3/m2)"""
	    psi_w = self.psi_wf(vw, cw)
	    return min(vw - self.gwf(psi_w)*(psi_w - (ev*(1. - self.F_CAP))/(lai*gp) - psi_l)*dt/10.**6, self.VWT)
	def fBal(self, params, soil, photo, phi, ta, qa, c1, s, lai, gp, ared, zr, psi_w):
	    psi_l, tl = params
	    if lai < 1.: # assumes only a portion of solar radiation is absorbed by crops
    		return (phi*lai - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)/1000000.,  \
    			self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)\
    			-(self.gsrfp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l) + lai*self.gwf(psi_w)*(psi_w - psi_l))/ \
    			(1. + (self.gsrfp(soil, s, gp, lai, zr)*(1. - self.F_CAP))/(lai*gp) + (self.gwf(psi_w)*(1. - self.F_CAP))/gp))
	    else:
	    	return (phi - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)/1000000., 
			self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)
			-(self.gsrfp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l) + lai*self.gwf(psi_w)*(psi_w - psi_l))/
			(1. + (self.gsrfp(soil, s, gp, lai, zr)*(1. - self.F_CAP))/(lai*gp) + (self.gwf(psi_w)*(1. - self.F_CAP))/gp))

class Amari(object):
	NAME = 'A. mari'
	PTYPE = C3

	ZR = 0.5
	LAI = 4.
	GCUT = 0.
	GA = 61.
	RAIW = 10.
	GPMAX = 2.

	GWMAX = 0.01
	VWT = .002
	CAP = 0.15

	RD0 = 4.93
	HAV = 62000.
	HDV = 202900.
	VCMAX0 = 83.
	JMAX0 = 132.
	PSILA0 = -3.5
	PSILA1 = -0.05

class Taest(object):
	NAME = 'T. aest'
	PTYPE = C3

	ZR = 0.75
	LAI = 5.
	GCUT = 0.3
	GA = 61.
	RAIW = 5.6
	GPMAX = 11.7

	GWMAX = 0.
	VWT = .000001
	CAP = 0.15

	RD0 = 4.93
	HAV = 62000.
	HDV = 202900.
	VCMAX0 = 83.
	JMAX0 = 132.
	PSILA0 = -2.
	PSILA1 = -0.7

class Zmays(object):
	NAME = 'Z. mays'
	PTYPE = C4

	ZR = 0.5 # most of these are old values for sorghum
	LAI = 5.
	GCUT = 0.1802
	GA = 61.
	RAIW = 5.6
	GPMAX = 0.13

	GWMAX = 0.
	VWT = .000001
	CAP = 0.15

	VCMAX0 = 39.
	JMAX0 = 180.
	PSILA0 = -1.9
	PSILA1 = -0.3

class Sbico(object):
	NAME = 'S. bico'
	PTYPE = C4

	ZR = 0.5
	LAI = 5.
	GCUT = 0.1802
	GA = 61.
	RAIW = 5.6
	GPMAX = 0.13

	GWMAX = 0.
	VWT = .000001
	CAP = 0.15

	VCMAX0 = 39.
	JMAX0 = 180.
	PSILA0 = -1.8
	PSILA1 = -0.5


class Oficu(object):
	NAME = 'O. ficu'
	PTYPE = CAM

	ZR = 0.3
	LAI = 3.5
	GCUT = 0.
	GA = 324.
	RAIW = 3.
	GPMAX = .4

	GWMAX = .02
	VWT = .0113
	CAP = 0.83

	VCMAX0 = 20. # old value 13
	JMAX0 = 40. # old value 26
	PSILA0 = -3.
	PSILA1 = -0.5

	MMAX = 190000000. # max concentration of malic acid (umol/m^3)
	AMMAX = 13.5  # rate of malic acid storage flux (umol/(m^2 s)

class Atequ(object):
	NAME = 'A. tequ'
	PTYPE = CAM

	ZR = 0.3
	LAI = 6.
	GCUT = 0.
	GA = 61.
	RAIW = 3.
	GPMAX = .04

	GWMAX = .002
	VWT = .00415
	CAP = 0.27

	VCMAX0 = 19.5
	JMAX0 = 39.
	MMAX = 130000000. # max concentration of malic acid (umol/m^3)
	AMMAX =  11.1 # rate of malic acid storage flux (umol/(m^2 s)
	PSILA0 = -3.
	PSILA1 = -0.5
	capOn = True

class Clusia(object):
	NAME = 'Clusia'
	PTYPE = CAM

	ZR = 0.4
	# ZR = 0.1
	LAI = 2.
	GCUT = 0.
	GA = 95.
	RAIW = 6.
	GPMAX = 2.

	GWMAX = 0.01
	VWT = .00058 # max. storage depth (m3/m2 leaf area)
	CAP = 0.1


	VCMAX0 = 3.15
	JMAX0 = 6.3
	MMAX = 92000000.
	AMMAX = 1.1
	PSILA0 = -1.5
	PSILA1 = -.5

	VLT = .0011 # water storage depth in leaves (m3/m2 leaf area)
	LCAP = 0.1 # hydraulic capacitance of leaf tissue (MPa)

class Pmenz(object):
	NAME = 'P. menz'
	PTYPE = C3

	ZR = 0.65
	LAI = 8.4
	GCUT = .007
	GA = 324.
	RAIW = 10.
	GWMAX = .005
	GPMAX = 0.056
	VWT = 0.27/LAI
	CAP = 0.15

	VCMAX0 = 57.7
	JMAX0 = 98.5
	PSILA0 = -3.
	PSILA1 = -0.5
	capOn = True

class FacCAM(Photo):
	# A1 = 0.6*15.
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
	Z0 = .55  # Initial value of z (-)
	M0 = 0. # Initial Malic Acid Carbon Concentration (umol/m^3)
	TH = 302.65 # High temperature for CAM model (K)
	TW = 283.15 # Low temperature for CAM model (K)
	def __init__(self, species, atm):
		Photo.__init__(self, "FacCAM", species)
		self.pmode = "C3"
		self.A1 = 15.
		self.MUPPER = species.MMAX
		self.AMMAX = species.AMMAX
		self.ca = atm.ca
		self.cs = atm.ca
		self.ci = self.ciNew(atm.cs, atm.ta, atm.qa)
		self.cm = self.cmNew(atm.cs, atm.ta, atm.qa)
		self.z = self.Z0
		self.m = self.M0
		# self.cc = self.ccNew(self.cs, atm.ta, atm.qa, self.z, self.m)
		self.cx = self.cm
		self.a_a = []
		
	def update(self, atm, psi_l, tl, dt):
		if psi_l < -1.:
			self.pmode = "CAM"
			self.A1 = 0.6*15.
			self.mmax = self.MUPPER
			if psi_l < -2.5:
				self.mmax = self.MUPPER
		else:
			pass
		self.ci = self.ciNew(self.cs, atm.ta, atm.qa)
		self.cm = self.cmNew(self.cs, atm.ta, atm.qa)
		if self.pmode == "CAM":
			self.cc = self.ccNew(self.cs, atm.ta, atm.qa, self.z, self.m)
			self.z = self.zNew(atm.phi, self.m, self.z, tl, dt) 
			self.m = self.mNew(atm.phi, psi_l, self.cc, tl, self.z, self.m, dt)
		else:
			pass
		self.cx = self.cxNew()
		self.a = self.an(atm.phi, psi_l, tl, self.cx, self.ared)
		self.a_a.append(self.a)

	def output(self):
		return {'a': self.a_a}

	def cxNew(self):
		if self.pmode == "CAM":
			return self.cc
		else:
			return self.cm
	def an(self, phi, psi_l, tl, ci, ared): 
	    """Photosynthetic rate, per unit leaf area (umol/(m^2s))"""
	    if self.pmode == "CAM":
	    	return self.a_sc(phi, psi_l, tl, ci, self.z, self.m, ared) + self.a_sv(phi, tl, psi_l, self.z, self.m) 
	    else:
	    	return self.a_psilc02(psi_l)*self.a_phiciTl(phi, ci, tl, ared)
		
	def a_sc(self, phi, psi_l, tl, ci, z, m, ared):
	    """Flux from stomata to Calvin cycle (umol/(m^2s))"""
	    return max(0, self.a_psilc02(psi_l)*(self.a_phiciTl(phi, ci, tl, ared) - self.r_dc(phi, tl))*(1. - self.f_c(z, m)))
	def r_dv(self, phi, tl):
	    """Flux of dark respiration to vacuole (umol/(m^2s))"""
	    # return self.r_d(tl)*exp(-phi)
	    return 0.
	def r_dc(self, phi, tl):
	    """Flux of dark respiration to calvin cycle (umol/(m^2s))"""
	    # return self.r_d(tl)*(1. - exp(-phi))
	    return 0.
	def f_o(self, z):
	    """Circadian order function (-)"""
	    return exp(-(z/self.MU)**self.CIRC_3)
	def f_m(self, z, m, tl):
	    """Malic acid storage function"""
	    return self.f_o(z)*(self.mmax*((self.TH - tl)/(self.TH - self.TW)*(1. - self.ALPHA_2) + self.ALPHA_2) - m)/(self.ALPHA_2*self.mmax*((self.TH - tl)/\
	        (self.TH - self.TW)*(1. - self.ALPHA_2) + self.ALPHA_2) + (self.mmax*((self.TH - tl)/(self.TH - self.TW)*(1. - self.ALPHA_2) + self.ALPHA_2) - m))
	def f_c(self, z, m):
	    """Carbon circadian control function"""
	    return (1. - self.f_o(z))*m/(self.ALPHA_1*self.mmax + m)
	def a_sv(self, phi, tl, psi_l, z, m):
		"""Flux from stomata to vacuole (umol/(m^2s))"""
		if self.mmax*((self.TH - tl)/(self.TH - self.TW)*(1. - self.ALPHA_2) + self.ALPHA_2) > m and (1. - self.K*(tl - self.TOPT)**2.) >0:
			return (self.AMMAX*(1. - self.K*(tl - self.TOPT)**2.) - self.r_dv(phi, tl))*self.f_m(z, m, tl)*self.a_psilc02(psi_l)
		else:
			return 0.
	def a_vc(self, phi, cc, tl, z, m):
	    """Flux from vacuole to calvin cycle (umol/(m^2s))"""
	    return (self.a_phiciTl(phi, cc, tl, 1.) - self.r_dc(phi, tl))*self.f_c(z, m)
	def m_e(self, z, m, tl, phi): 
	    """Malic acid equilibrium value"""
	    if phi>0.:
	        return self.mmax*(self.CIRC_1*((self.TH - tl)/(self.TH - self.TW) + 1.)*(self.BETA*(z - self.MU))**3. - self.BETA*(self.TH - tl)/(self.TH - self.TW)*(z - self.MU) + \
	            self.CIRC_2*(self.TH - tl)/(self.TH - self.TW) -(1- self.f_o(z))*(1-m/(m+self.ALPHA_1*self.mmax))) 
	    else:
	        return self.mmax*(self.CIRC_1*((self.TH - tl)/(self.TH - self.TW) + 1.)*(self.BETA*(z - self.MU))**3. - self.BETA*(self.TH - tl)/(self.TH - self.TW)*(z - self.MU) + \
	            self.CIRC_2*(self.TH - tl)/(self.TH - self.TW)+ (1-self.f_o(z))) 
	def zNew(self, phi, m, z, tl, dt):
	    return max(0, dt*(m - self.m_e(z, m, tl, phi))/(self.mmax*60.*self.TR) + z)
	def ccNew(self, cs, ta, qa, z, m):
		"""CO2 concentration in mesophyll cytosol resulting from malic acid decarboxylation (ppm)"""
		return self.cmNew(cs, ta, qa) + self.f_c(z, m)*self.C0
	def mNew(self, phi, psi_l, cc, tl, z, m, dt): 
		"""Malic acid concentration"""
		return max(((dt/ self.VCM)*(self.a_sv(phi, tl, psi_l, z, m) - self.a_vc(phi, cc, tl, z, m) + self.r_dv(phi, tl))) + m, 0.)

