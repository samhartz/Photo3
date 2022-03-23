from math import exp, pi, sqrt, log
from scipy.optimize import fsolve
from sympy import *
import numpy as np
from dics import *
from functions import *
from photosynthesis import *

class Oficu(object):
	"""Opuntia ficus-indica (prickly pear)"""
	NAME = 'O. ficu'
	PTYPE = CAM

	ZR = 0.3
	LAI = 3.5
	GCUT = 0.
	#GA = 324.
	GA = 30.
	RAIW = 3.
	GPMAX = .4

	GWMAX = .02
	VWT = .0113
	CAP = 0.83

	VCMAX0 = 18. #20. # old value 13
	JMAX0 = 36. #40. # old value 26
	PSILA0 = -3.
	PSILA1 = -0.5

	MMAX = 230000000. # 190000000. # max concentration of malic acid (umol/m^3) 230
	AMMAX = 14. # 13.5  # rate of malic acid storage flux (umol/(m^2 s) 14

class Pmenz(object):
	"""Pseudotsuga menziesii (Douglas fir)"""
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
	
class Sbico(object):
	"""Sorghum bicolor"""
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
