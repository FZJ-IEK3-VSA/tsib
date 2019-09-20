[![Build Status](https://travis-ci.com/FZJ-IEK3-VSA/tsib.svg?branch=master)](https://travis-ci.com/FZJ-IEK3-VSA/tsib) [![Version](https://img.shields.io/pypi/v/tsib.svg)](https://pypi.python.org/pypi/tsib)

<a href="https://www.fz-juelich.de/iek/iek-3/EN/Forschung/_Process-and-System-Analysis/_node.html"><img src="https://www.fz-juelich.de/SharedDocs/Bilder/INM/INM-1/EN/FZj_Logo.html" alt="Forschungszentrum Juelich Logo" width="230px"></a> 

# tsib - Time Series Initialization for Buildings

tsib is a python package that builds up on different databases and models for creating consistent demand time series of residential buildings. This could be either occupancy behavior, electricity demand or heat demand time series.


If you want to use tsam in a published work, please [**kindly cite following publication**](TODO) which applies the time series for Germany. 


## Features
* flexible configuration of single buildings by different input arguments
* simulation of the stochastic occupancy behavior (CREST)
* derivation of the electric device load or the demand for thermal comfort
* prediction of the heat load based on a 5R1C (Schütz)

## Installation
Directly install via pip as follows:

	pip install tsib

Alternatively, clone a local copy of the repository to your computer

	git clone https://github.com/FZJ-IEK3-VSA/tsib.git
	
Then install tsam via pip as follow
	
	cd tsib
	pip install . 
	
Or install directly via python as 

	python setup.py install
	
In order to use the 5R1C thermal building model, make sure that you have installed a MILP solver. As default solver GLPK is used. Nevertheless, in case you have access to a license we recommend commercial solvers (e.g. Gurobi or CPLEX) since they have a better performance.
	
To get flexible weather data from the Climate Data Store, register [here](https://cds.climate.copernicus.eu/api-how-to) and follow the instructions to get an own key. Make sure that you have agreed on the [license terms](https://cds.climate.copernicus.eu/cdsapp/#!/terms/licence-to-use-copernicus-products).

	
## Examples

### Basic workflow

A small example how tsam can be used is decribed as follows
```python
	import pandas as pd
	import tsib
```


### Detailed examples

A [**first example**](/examples/aggregation_example.ipynb) shows the capabilites of tsib as jupyter notebook. 

The example time series are based on a department [publication](https://www.mdpi.com/1996-1073/10/3/361) and the [test reference years of the DWD](https://www.dwd.de/DE/leistungen/testreferenzjahre/testreferenzjahre.html).

## License

MIT License

Copyright (C) 2016-2019 Leander Kotzur (FZJ IEK-3), Kevin Knosala (FZJ IEK-3), Peter Stenzel (FZJ IEK-3), Peter Markewitz (FZJ IEK-3), Martin Robinius (FZJ IEK-3), Detlef Stolten (FZJ IEK-3)

You should have received a copy of the MIT License along with this program.
If not, see https://opensource.org/licenses/MIT

## About Us 
<a href="http://www.fz-juelich.de/iek/iek-3/EN/Forschung/_Process-and-System-Analysis/_node.html"><img src="https://www.fz-juelich.de/SharedDocs/Bilder/IEK/IEK-3/Abteilungen2015/VSA_DepartmentPicture_2019-02-04_459x244_2480x1317.jpg?__blob=normal" width="400px" alt="Abteilung VSA"></a> 

We are the [Techno-Economic Energy Systems Analysis](https://www.fz-juelich.de/iek/iek-3/EN/Forschung/_Process-and-System-Analysis/_node.html) department at the [Institute of Energy and Climate Research: Electrochemical Process Engineering (IEK-3)](https://www.fz-juelich.de/iek/iek-3/EN/Home/home_node.html) belonging to the [Forschungszentrum Jülich](https://www.fz-juelich.de/). Our interdisciplinary department's research is focusing on energy-related process and systems analyses. Data searches and system simulations are used to determine energy and mass balances, as well as to evaluate performance, emissions and costs of energy systems. The results are used for performing comparative assessment studies between the various systems. Our current priorities include the development of energy strategies, in accordance with the German Federal Government’s greenhouse gas reduction targets, by designing new infrastructures for sustainable and secure energy supply chains and by conducting cost analysis studies for integrating new technologies into future energy market frameworks.

## Acknowledgement

This work was supported by the Helmholtz Association under the Joint Initiative ["Energy System 2050   A Contribution of the Research Field Energy"](https://www.helmholtz.de/en/research/energy/energy_system_2050/).

<a href="https://www.helmholtz.de/en/"><img src="https://www.helmholtz.de/fileadmin/user_upload/05_aktuelles/Marke_Design/logos/HG_LOGO_S_ENG_RGB.jpg" alt="Helmholtz Logo" width="200px" style="float:right"></a>

