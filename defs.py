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
