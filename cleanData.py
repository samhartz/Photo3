from scipy.optimize import fsolve
from sympy import *
import Tkinter as tk
import numpy as np
import pandas as pd
from math import exp, pi, sqrt
import matplotlib.pyplot as plt
from datetime import date
import xlrd

#Code to format weather data from NSRDB data files (download from https://maps.nrel.gov/nsrdb-viewer//?aL=UdPEX9)
# data must be formatted as:
# Year Month  Day   Hour   Minute  GHI Temperature Relative Humidity

# data cannot contain gaps larger than 23 hours
# Excel sheet containing data must not additional text beyond headers

xl = pd.ExcelFile("sample_data/TempleApril2015.xlsx")
data = xl.parse("Sheet1")

# convert time array to simple vector of datenumbers
timevecs = data[['Year', 'Month', 'Day', 'Hour', 'Minute']] 
timestamps = pd.to_datetime(timevecs)
newdf = data[['GHI', 'Temperature', 'Relative Humidity']].set_index(timestamps)
converted = newdf.asfreq('30Min') #May change time increment here to match model timestepM
converted = converted.interpolate()
converted['GHI'].plot()
plt.show()

writer = pd.ExcelWriter('sample_data/TempleApril2015Interp30.xlsx')
converted.to_excel(writer,'Sheet1')
writer.save()
