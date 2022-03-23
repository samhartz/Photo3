from math import exp, pi, sqrt, log
from scipy.optimize import fsolve
from sympy import *
import numpy as np
from dics import *
from functions import *

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

	def evfPen(self, photo, phi, ta, psi_l, qa, tl, ci,  lai, ared):
		"""Penman-Monteith transpiration (um/sec)"""
		GAMMA_W = (P_ATM*CP_A)/(.622*LAMBDA_W)
		def delta_s(ta):
			return esat(ta)*(C_SAT*B_SAT)/(C_SAT + ta -273)**2
		def drh(ta, qa):
			return VPD(ta, qa)*.622/P_ATM
		GMGSRATIO = 1.
		gsCAM = self.gsw(photo, phi, ta, psi_l, qa, tl, ci, ared)*(1.6*(1.+GMGSRATIO)/(1.6+GMGSRATIO))

		return ((LAMBDA_W*GAMMA_W*self.GA/1000.*RHO_A*drh(ta, qa) + delta_s(ta)*phi)*(R*ta/P_ATM)*gsCAM*1000000.*lai)/ \
		(RHO_W*LAMBDA_W*(GAMMA_W*(self.GA/1000. + (R*ta/P_ATM)*gsCAM*lai) + (R*ta/P_ATM)*gsCAM*lai*delta_s(ta)))
		
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
		# if psi_l<-10:
		# 	return 0.
		# else:
		# 	return self.GPMAX*exp(-(-psi_l/2.)**2.)
		# return self.GPMAX
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
		if self.qi(self.tl, self.psi_l) < atm.qa:
			self.psi_l = psi_i(atm.ta, atm.qa)
			self.tl = fsolve(self.fBal_psil_known, (290.), args= (self.psi_l, soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr))
		self.gp = self.gpf(self.psi_l)
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, photo.ared)
	def update(self, atm, soil, photo, dt):
		self.psi_l, self.tl = fsolve(self.fBal, (-1., 290.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr))
		if self.qi(self.tl, self.psi_l) < atm.qa:
			self.psi_l = psi_i(atm.ta, atm.qa)
			self.tl = fsolve(self.fBal_psil_known, (290.), args= (self.psi_l, soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr))
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
	def fBal_psil_known(self, p, psi_l, soil, photo, phi, ta, qa, c1, s, lai, gp, ared, zr):
		tl = p

		if lai < 1.: # assumes only a portion of solar radiation is absorbed by crops
			return phi - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*self.gsrp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l)/1000000.
		else:
			return phi - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*self.gsrp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l)/1000000.


class HydroCap(Hydro):
	F_CAP = 0.5


# # Below is the protocol for solving with Penman-Monteith (no solving energy balance!)
# 	def __init__(self, species, atm, soil, photo, vwi):
# 		Hydro.__init__(self, species)
# 		self.GWMAX = species.GWMAX
# 		self.VWT = species.VWT
# 		self.CAP = species.CAP
# 		self.vw = vwi*self.VWT
# 		self.tl = atm.ta
# 		self.psi_l = fsolve(self.fBalPen, (-1.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr, self.vw, self.tl))
# 		self.gp = self.gpf(self.psi_l)
# 		self.ev = self.evfPen(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, photo.ared)
# 		self.tl = self.tlPen(atm.ta, atm.phi, self.ev)
# 		self.vw_a = []

# 	def update(self, atm, soil, photo, dt):
# 		self.psi_l = fsolve(self.fBalPen, (-1.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr, self.vw, self.tl))
# 		self.gp = self.gpf(self.psi_l)
# 		self.vw = self.vwf(self.vw, self.ev, self.gp, self.psi_l, self.lai, dt)
# 		self.ev = self.evfPen(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, photo.ared)
# 		self.tl = self.tlPen(atm.ta, atm.phi, self.ev)
# 		self.qs = self.qsf(self.vw, self.ev, self.gp, self.psi_l, self.lai, dt)
# 		self.psi_l_a.append(self.psi_l)
# 		self.gp_a.append(self.gp)
# 		self.gsv_a.append(self.gsw(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, photo.ared))
# 		self.tl_a.append(self.tl) 
# 		self.ev_a.append(self.ev)
# 		self.vw_a.append(self.vw)

	def __init__(self, species, atm, soil, photo, vwi):
		Hydro.__init__(self, species)
		self.GWMAX = species.GWMAX
		self.VWT = species.VWT
		self.CAP = species.CAP
		self.vw = vwi*self.VWT
		self.psi_l, self.tl = fsolve(self.fBal, (-1., 290.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr, self.vw))
		if self.qi(self.tl, self.psi_l) < atm.qa:
			self.psi_l = psi_i(atm.ta, atm.qa)
			self.tl = fsolve(self.fBal_psil_known, (290.), args= (self.psi_l, soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr, self.vw))
		self.gp = self.gpf(self.psi_l)
		self.ev = self.evf(photo, atm.phi, atm.ta, self.psi_l, atm.qa, self.tl, photo.cx, self.lai, photo.ared)
		self.vw_a = []

	def update(self, atm, soil, photo, dt):
		self.psi_l, self.tl = fsolve(self.fBal, (-1., 290.), args= (soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr, self.vw))
		if self.qi(self.tl, self.psi_l) < atm.qa:
			self.psi_l = psi_i(atm.ta, atm.qa)
			self.tl = fsolve(self.fBal_psil_known, (290.), args= (self.psi_l, soil, photo, atm.phi, atm.ta, atm.qa, photo.cx, soil.s, self.lai, self.gp, photo.ared, self.zr, self.vw))
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
		if gp < 0.000001:
			return vw
		else:
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
		# return self.GWMAX
		return self.GWMAX*exp(-(-psi_w/2.)**2.)
		#return GWMAX[species]*(vw/VWT[species])**4. 
	def gsrfp(self, soil, s, gp, lai, zr):
		"""Soil-root-plant fraction conductance, per unit ground area (um/(s-MPa))"""
		return (lai*self.gsr(soil, s, zr)*gp/self.F_CAP)/(self.gsr(soil, s, zr) +  lai*gp/self.F_CAP)
	def fBal(self, params, soil, photo, phi, ta, qa, c1, s, lai, gp, ared, zr, vw):
		psi_l, tl = params
		psi_w = self.psi_wf(vw)

		if gp == 0.:
			return (psi_l - psi_i(ta, qa), \
				phi*lai - self.shf(tl, ta, lai))
		elif lai < 1.: # assumes only a portion of solar radiation is absorbed by crops
			return (phi*lai - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)/1000000.,  \
				self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)\
				-(self.gsrfp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l) + lai*self.gwf(psi_w)*(psi_w - psi_l))/ \
				(1. + (self.gsrfp(soil, s, gp, lai, zr)*(1. - self.F_CAP))/(lai*gp) + (self.gwf(psi_w)*(1. - self.F_CAP))/gp))
		else:
			# energy balance, phi = shf + lambda rho evf
			# water balance, evf = (gsrfp(psi_s - psi_l) + gw(psi_w-psi_l))/(1 + etc....)
			return (phi - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)/1000000., \
				self.evf(photo, phi, ta, psi_l, qa, tl, c1, lai, ared)\
				-(self.gsrfp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l) + lai*self.gwf(psi_w)*(psi_w - psi_l))/ \
				(1. + (self.gsrfp(soil, s, gp, lai, zr)*(1. - self.F_CAP))/(lai*gp) + (self.gwf(psi_w)*(1. - self.F_CAP))/gp))
	def fBal_psil_known(self, p, psi_l, soil, photo, phi, ta, qa, c1, s, lai, gp, ared, zr, vw):
		tl = p
		psi_w = self.psi_wf(vw)
		if lai < 1.: # assumes only a portion of solar radiation is absorbed by crops
			return phi - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*(self.gsrfp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l) + lai*self.gwf(psi_w)*(psi_w - psi_l))/ \
				(1. + (self.gsrfp(soil, s, gp, lai, zr)*(1. - self.F_CAP))/(lai*gp) + (self.gwf(psi_w)*(1. - self.F_CAP))/gp)/1000000.
		else:
			return phi - self.shf(tl, ta, lai) -LAMBDA_W*RHO_W*(self.gsrfp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l) + lai*self.gwf(psi_w)*(psi_w - psi_l))/ \
				(1. + (self.gsrfp(soil, s, gp, lai, zr)*(1. - self.F_CAP))/(lai*gp) + (self.gwf(psi_w)*(1. - self.F_CAP))/gp)/1000000.
	def tlPen(self, ta, phi, ev):
		H = phi - LAMBDA_W*RHO_W*ev/1000000.
		return ta + H*1000/(CP_A*RHO_A*self.GA)

	def fBalPen(self, param, soil, photo, phi, ta, qa, c1, s, lai, gp, ared, zr, vw, tl):
		psi_l = param
		psi_w = self.psi_wf(vw)
		return self.evfPen(photo, phi, ta, psi_l, qa, tl, c1, lai, ared) \
		- (self.gsrfp(soil, s, gp, lai, zr)*(soil.psi_s(s) - psi_l) + lai*self.gwf(psi_w)*(psi_w - psi_l))/ \
		(1. + (self.gsrfp(soil, s, gp, lai, zr)*(1. - self.F_CAP))/(lai*gp) + (self.gwf(psi_w)*(1. - self.F_CAP))/gp)
