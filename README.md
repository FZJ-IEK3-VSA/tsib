[![Build Status](https://img.shields.io/gitlab/pipeline/l-kotzur/tsib/master.svg)](https://gitlab.com/l-kotzur/tsib/pipelines)
[![Version](https://img.shields.io/pypi/v/tsib.svg)](https://pypi.python.org/pypi/tsib)

<a href="https://www.fz-juelich.de/iek/iek-3/EN/Forschung/_Process-and-System-Analysis/_node.html"><img src="https://www.fz-juelich.de/SharedDocs/Bilder/INM/INM-1/EN/FZj_Logo.jpg?__blob=normal" alt="Forschungszentrum Juelich Logo" width="230px"></a> 

# tsib - Time Series Initialization for Buildings

tsib is a python package that builds up on different databases and models for creating consistent demand and production time series of residential buildings. This could be either occupancy behavior, electricity demand or heat demand time series as well as photovoltaic (PV) and solar thermal production time series.


If you want to use tsib in a published work, please [**cite following publication**](http://juser.fz-juelich.de/record/858675) which applies tsib for the creation of time series for residential buildings in Germany. 


## Features
* flexible configuration of single buildings by different input arguments
* simple building definition based on an archetype building catalogue
* consideration of the occupancy behavior
* derivation of the electric device load or the demand for thermal comfort
* calculation of the heat load based on a thermal building model
* provision of location specific time series for solar irradiation and temperature based on weather data


## Applied databases and models
tsib is a flexible tool which allows the use of different models and databases for the generation of time series for buildings. In Version 0.1.0 the following databases and models are included is tsib:
* [CREST](https://www.lboro.ac.uk/research/crest/demand-model/) demand model for the simulaton of the occupancy behavior
* [5R1C](https://www.sciencedirect.com/science/article/abs/pii/S0306261916314933) thermal building model 
* [pvlib](https://github.com/pvlib/pvlib-python) for solar irradiance calculation and photovoltaic simulation
* [TABULA/EPISCOPE](http://episcope.eu/) archetype building catalogue
* [DWD Testreferenzjahre](https://www.dwd.de/DE/leistungen/testreferenzjahre/testreferenzjahre.html)  for providing weather data


## Installation
Directly install via pip as follows:

	pip install tsib

Alternatively, clone a local copy of the repository to your computer

	git clone https://github.com/FZJ-IEK3-VSA/tsib.git
	
Then install tsib via pip as follow
	
	cd tsib
	pip install . 
	
Or install directly via python as 

	python setup.py install
	
In order to use the 5R1C thermal building model, make sure that you have installed a MILP solver. As default solver coin-cbc is used which can either installed by

	sudo apt-get install coinor-cbc

or for Anaconda under windows as

	conda install -c conda-forge coincbc

. Other solvers can be defined by defining the environment variable $SOLVER. 
	
To get flexible weather data from the Climate Data Store, register [here](https://cds.climate.copernicus.eu/api-how-to) and follow the instructions to get an own key. Make sure that you have agreed on the [license terms](https://cds.climate.copernicus.eu/cdsapp/#!/terms/licence-to-use-copernicus-products).

	
## Examples

This [jupyter notebook](examples/showcase.ipynb) shows the capabilites of tsib to create all relevant time series. 


## License

MIT License

Copyright (C) 2016-2019 Leander Kotzur (FZJ IEK-3), Timo Kannengießer (FZK-IEK-3), Kevin Knosala (FZJ IEK-3), Peter Stenzel (FZJ IEK-3), Peter Markewitz (FZJ IEK-3), Martin Robinius (FZJ IEK-3), Detlef Stolten (FZJ IEK-3)

You should have received a copy of the MIT License along with this program.
If not, see https://opensource.org/licenses/MIT

## About Us 
<a href="http://www.fz-juelich.de/iek/iek-3/EN/Forschung/_Process-and-System-Analysis/_node.html"><img src="https://www.fz-juelich.de/SharedDocs/Bilder/IEK/IEK-3/Abteilungen2015/VSA_DepartmentPicture_2019-02-04_459x244_2480x1317.jpg?__blob=normal" width="400px" alt="Abteilung VSA"></a> 

We are the [Institute of Energy and Climate Research: Techno-Economic Energy Systems Analysis (IEK-3)](https://www.fz-juelich.de/iek/iek-3/EN/Forschung/_Process-and-System-Analysis/_node.html) belonging to the [Forschungszentrum Jülich](https://www.fz-juelich.de/). Our interdisciplinary department's research is focusing on energy-related process and systems analyses. Data searches and system simulations are used to determine energy and mass balances, as well as to evaluate performance, emissions and costs of energy systems. The results are used for performing comparative assessment studies between the various systems. Our current priorities include the development of energy strategies, in accordance with the German Federal Government’s greenhouse gas reduction targets, by designing new infrastructures for sustainable and secure energy supply chains and by conducting cost analysis studies for integrating new technologies into future energy market frameworks.

## Acknowledgement

This work was supported by the Helmholtz Association under the Joint Initiative ["Energy System 2050   A Contribution of the Research Field Energy"](https://www.helmholtz.de/en/research/energy/energy_system_2050/).

<a href="https://www.helmholtz.de/en/"><img src="https://www.helmholtz.de/fileadmin/user_upload/05_aktuelles/Marke_Design/logos/HG_LOGO_S_ENG_RGB.jpg" alt="Helmholtz Logo" width="200px" style="float:right"></a>
