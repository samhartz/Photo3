try:
    import tkinter as tk
except ImportError:
    import Tkinter as tk
# import tkinter as tk
try:
    import tkFileDialog as filedialog
except ImportError:
    import tkinter.filedialog as filedialog
from defs import *

root = tk.Tk()  # This creates a window, but it won't show up

# Default settings
weatherOption = "ACTUAL"  # options: ACTUAL, AVG, BLM, STEP, CONST, VPD
lightOption = "ACTUAL" # options: ACTUAL, STEP, PAR
durationInp = tk.IntVar(value = 10) # days (must be an integer)
sInit = tk.DoubleVar(value = 0.4)
soilVariable = tk.StringVar(value="Sandy loam") # default value
capVariable = tk.StringVar(value = "Yes")
spVariable = tk.StringVar(value = "Opuntia ficus-indica") # default value
rainVariable = tk.StringVar(value = "Drydown") # options: "Drydown", "Constant", "Input rainfall", "Stochastic", "Input VWC"
weatherFile = tk.StringVar(value="sample_data\TempleApril2015Interp30.xlsx") # default value for the location of the weather data
repeatWeather = 1 # option to repeat the data in weatherFile x number of times...default value is 1 for no repeat
resultsFile = tk.StringVar(value = 'sample_output\guitest') # default value for the location where results are saved

capOptions = {"No": HydroNC, "Yes": HydroCap}
speciesOptions = {"Triticum aestivum": Taest, "Sorghum bicolor": Sbico, "Opuntia ficus-indica": Oficu, "Agave tequilana": Atequ}
soilOptions = {"Sand": Sand, "Sandy loam": SandyLoam, "Loamy sand": LoamySand, "Loam": Loam, "Clay": Clay}
rainOptions = {"Drydown": DrydownSoil(), "Constant": ConstantSoil(), "Stochastic rainfall": StochasticSoil(1.5, .3)}

# Take inputs from GUI
def makemenu(root, caption, rowno, default, *options):
    lbl = tk.Label(root, text = caption)
    lbl.grid(row = rowno, column = 0)
    menu = tk.OptionMenu(root, default, *options)
    menu.grid(row = rowno, column = 1)

def makeentry(root, caption, rowno, default):
    lbl = tk.Label(root, text = caption)
    lbl.grid(row = rowno, column = 0)
    entry = tk.Entry(root, textvariable = default)
    entry.grid(row = rowno, column = 1)

def ok():
    global hydrology, species, sType, locs, duration, s0, soilDynamics, weatherLoc, resultsLoc
    duration = durationInp.get()
    s0 = sInit.get()
    hydrology = capOptions[capVariable.get()]
    species = speciesOptions[spVariable.get()]
    sType = soilOptions[soilVariable.get()]
    soilDynamics = rainOptions[rainVariable.get()]
    weatherLoc = weatherFile.get()
    resultsLoc = resultsFile.get()
    root.destroy()

def chooseinputfile():
    global weatherFile
    root.filename = filedialog.askopenfilename(title = "Select file",filetypes = ( ("Xlsx files","*.xlsx"), ("All files","*")))
    weatherFile.set(root.filename)
    print(weatherFile.get())

def chooseoutputfile():
    global resultsFile
    root.filename = filedialog.asksaveasfilename(title = "Select file",filetypes = (  ("Xlsx files","*.xlsx"), ("All files","*")))
    resultsFile.set(root.filename)
    print(resultsFile.get())

speciesmenu =  makemenu(root, 'Plant species:', 0, spVariable,  *speciesOptions.keys())
capmenu = makemenu(root, 'Plant water storage:', 1, capVariable,  *capOptions.keys())
wLbl = tk.Label(root, text = "Weather data file")
wLbl.grid(row = 2, column = 0)
weatherEntry = tk.Entry(root, textvariable = weatherFile)
weatherEntry.grid(row = 2, column = 1)
filebutton = tk.Button(root, text="Select file", command=chooseinputfile)
filebutton.grid(row=2, column = 2)
soilmenu = makemenu(root, 'Soil type:', 3, soilVariable, *soilOptions.keys())
rainmenu = makemenu(root, 'Soil moisture dynamics:', 4, rainVariable,  *rainOptions.keys())
durationinput = makeentry(root, 'Duration (days):', 5, durationInp)
moistinput = makeentry(root, 'Init. soil moisture (%):', 6, sInit)
svLbl = tk.Label(root, text = "Save results as")
svLbl.grid(row = 7, column = 0)
saveEntry = tk.Entry(root, textvariable = resultsFile)
saveEntry.grid(row = 7, column = 1)
savebutton = tk.Button(root, text="Select file", command=chooseoutputfile)
savebutton.grid(row=7, column = 2)
runbutton = tk.Button(root, text="Run", command=ok)
runbutton.grid(row=8, column = 1)
root.mainloop()                 # This command will tell the window to come out
