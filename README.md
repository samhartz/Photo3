# Photo3
The Photo3 model describes C3, C4, and CAM photosynthesis in a consistent manner using a model which is built on the Farquhar et al. model for carbon assimilation. The model incorporates soil and atmospheric conditions through a representation of the soil-plant-atmosphere continuum. Given soil moisture, air temperature, humidity and solar radiation, the model calculates net carbon assimilation and transpiration, among other variables of interest. Photo3 is currently parameterized for three representative species, one from each photosynthetic type: *Triticum aestivum* (C3), *Sorghum bicolor* (C4), and *Opuntia ficus-indica* (CAM).

# Model Execution

To run the model, simply download the essential files (main_gui.py, defs.py, functions.py, dics.py, gui.py) to the same folder and run the main file (main_gui.py). A graphical user interface (GUI) will display allowing the user to choose a plant species, soil type, soil moisture, duration of the simulation, and a data file containing weather inputs (solar radiation, temperature, and humidity). Results will be exported to the file location chosen in the user interface, and/or may be viewed directly from the command prompt or IDE. The sample_data folder contains sample weather inputs for a location in Temple, TX, and the sample_output folder contains results from a few simulations in this area. 

# Model Structure

The model is structured in an object-oriented fashion, using mixins to combine different photosynthetic, hydraulic, soil, and atmosphere sub-components. The main_gui.py script is the engine which runs the model, calling on the gui.py file to create a GUI which interprets the user's input to the model and creating a Simulation() object which is updated at each timestep according to the model inputs. The defs.py file contains the definitions of the classes and their functions, the dics.py file contains the model global variables, and the functions.py file contains the model global functions. Alternatively, the model can be executed by running the main.py script, which does not involve a GUI. The user may directly change the input variables within the file main.py.

An academic article describing the model details is available in Ecological Modelling:

  *Hartzell, S., Bartlett, M.S. and A. Porporato (2018) Unified representation of the C3, C4, and CAM photosynthetic pathways with the       Photo3 model. Ecological Modelling, doi: 10.1016/j.ecolmodel.2018.06.012.*

# Instructions for formatting weather data for the model input
Half hourly data for solar radiation, temperature, and specific humidity is needed to run the Photo3 model. Hourly data may be obtained from a variety of sources and then interpolated to the model timestep of 30 minutes using the script cleanData.py, which will produce a new excel file which can be read into the model. In order to use cleanData.py, the data must first be formatted with the following columns: Year, Month, Day, Hour, Minute, GHI, Temperature, Relative Humidity. Data cannot contain gaps larger than 23 hours.

The final weather data file supplied to the model should have the headings: Temperature, Relative Humidity, GHI, and should match the model timestep of 30 minutes (see example files in the sample_data folder).

# Publications

Miller, G., Hartzell, S., and A. Porporato. Ecohydrology of epiphytes: modeling water balance, CAM photosynthesis, and their climate impacts. In press, Ecohydrology.

Hartzell, S., Bartlett, M.S., Inglese, P., Consoli, S., Yin, J. and A. Porporato (2020) Modeling nonlinear dynamics of CAM productivity and water use for global predictions. Plant, Cell and Environment, 44(1) [10.1111/pce.13918].

Hartzell, S., Bartlett, M.S. and A. Porporato (2018) Unified representation of the C3, C4, and CAM photosynthetic pathways with the Photo3 model. Ecological Modelling, 384 [10.1016/j.ecolmodel.2018.06.012].

# Model requirements
Photo3 was developed in Python 3.6 with the following packages: SciPy, NumPy, Pandas, tkinter, SymPy, Matplotlib. We suggest intalling a Python distribution such as [Anaconda][An] to meet these requirements. 

[An]: https://www.continuum.io/downloads

# License

Copyright (C) 2021, Samantha Hartzell  

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.



