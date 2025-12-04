from math import exp, pi, sqrt, log
from scipy.optimize import fsolve
from sympy import *
import numpy as np
from dics import *
from functions import *
from photosynthesis import *

class Oficu(object):
	"""Opuntia ficus-indica (prickly pear)"""
	NAME = 'O. ficu' # species abbreviation (first letter genus, first four letters species)
	PTYPE = CAM # photosynthetic pathway (C3, C4, or CAM)

	# Plant hydraulic paramters
	ZR = 0.3 # rooting depth (m)
	LAI = 3.5 # leaf area index (-)
	GCUT = 0. # cuticular conductance to water vapor (mm/s)
	#GA = 324.
	GA = 30. # atmospheric conductance per unit ground area (mm/s)
	RAIW = 3. # well-watered root area index (-)
	GPMAX = .4 # maximum plant stem hydraulic conductance (um/MPa/s)

	# Plant water storage parameters (only needed for plant hydraulics with capacitance option)
	GWMAX = .02 # max. conductance between water storage tissue and plant xylem (um/MPa/s)
	VWT = .0113 # maximum water storage depth (m)
	CAP = 0.83 # hydraulic capacitance (MPa-1)

	# Photosynthetic parameters
	VCMAX0 = 18. # maximum carboxylation capacity (umol/m2/sec)
	JMAX0 = 36. # maximum electron transport capacity (umol/m2/sec)
	PSILA0 = -3. # leaf water potential at point of full stomatal closure (MPa)
	PSILA1 = -0.5 # leaf water potential at onset of stomatal closure (MPa)
	
	# CAM-specific parameters (only needed for CAM photosynthetic species)
	MMAX = 230000000.  # max concentration of malic acid (umol/m3) 
	AMMAX = 14. # rate of malic acid storage flux (umol/m2/s) 

class Pmenz(object):
	"""Pseudotsuga menziesii (Douglas fir)"""
	NAME = 'P. menz' # species abbreviation (first letter genus, first four letters species)
	PTYPE = C3 # photosynthetic pathway (C3, C4, or CAM)

	# Plant hydraulic paramters
	ZR = 0.65 # rooting depth (m)
	LAI = 8.4 # leaf area index (-)
	GCUT = .007 # cuticular conductance to water vapor (mm/s)
	GA = 324. # atmospheric conductance per unit ground area (mm/s)
	RAIW = 10. # well-watered root area index (-)
	GPMAX = 0.056 # maximum plant stem hydraulic conductance (um/MPa/s)

	# Plant water storage parameters (only needed for plant hydraulics with capacitance option)
	capOn = True
	GWMAX = .005 # max. conductance between water storage tissue and plant xylem (um/MPa/s)
	VWT = 0.27/LAI # maximum water storage depth (m)
	CAP = 0.15 # hydraulic capacitance (MPa-1)

	# Photosynthetic parameters
	VCMAX0 = 57.7 # maximum carboxylation capacity (umol/m2/sec)
	JMAX0 = 98.5 # maximum electron transport capacity (umol/m2/sec)
	PSILA0 = -3. # leaf water potential at point of full stomatal closure (MPa)
	PSILA1 = -0.5 # leaf water potential at onset of stomatal closure (MPa)
	
	
class Sbico(object):
	"""Sorghum bicolor"""
	NAME = 'S. bico'
	PTYPE = C4

	# Plant hydraulic paramters
	ZR = 0.5 # rooting depth (m)
	LAI = 5. # leaf area index (-)
	GCUT = 0.1802 # cuticular conductance to water vapor (mm/s)
	GA = 61. # atmospheric conductance per unit ground area (mm/s)
	RAIW = 5.6 # well-watered root area index (-)
	GPMAX = 0.13 # maximum plant stem hydraulic conductance (um/MPa/s)

	# Plant water storage parameters (only needed for plant hydraulics with capacitance option)
	GWMAX = 0. # max. conductance between water storage tissue and plant xylem (um/MPa/s)
	VWT = .000001 # maximum water storage depth (m)
	CAP = 0.15 # hydraulic capacitance (MPa-1)

	# Photosynthetic parameters
	VCMAX0 = 39. # maximum carboxylation capacity (umol/m2/sec)
	JMAX0 = 180. # maximum electron transport capacity (umol/m2/sec)
	PSILA0 = -1.8 # leaf water potential at point of full stomatal closure (MPa)
	PSILA1 = -0.5 # leaf water potential at onset of stomatal closure (MPa)

class Taest(object):
	"""Triticum aestivum (winter wheat)"""
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
