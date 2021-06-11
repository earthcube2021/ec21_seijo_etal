# click Binder badge to launch notebook on Binder (no need to install environment):
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/gseijo/EC_test/HEAD)
# Notebook: Investigating upper ocean variability during tropical cyclones and seasonal sea ice formation and melting: Argovis APIs exposed to co-locate oceanic and atmospheric datasets
# authors: Giovanni Seijo-Ellis, Donata Giglio, Sarah Purkey, Megan Scanderbeg, and Tyler Tucker
# Last edit: 06/11/2021 9:19 pm MT
# Version 02:
- Large portions of code were consolidated into pre-defined functions.
- All pre-defined functions were compiled into utilities.py and imported as part of the local library imports (Section 2.2). A list and description of the functions is found in Section 2.3 of the notebook.
- Progress bars were added for data download and plotting functions where useful.
# contact: giovanni.seijo@colorado.edu
#  This notebook has been submitted to the EarthCube 2021 Annual Meeting.
# Purpose:
Argovis is a web app and database that allows easy access to Argo profile observations of the global ocean and other earth science datasets using a browser and/or via Application Programming Interfaces (APIs). This notebook serves two main purposes: (i) introducing two new APIs available to access National Hurricane Center tropical cyclone (TC) track data and sea-ice concentration from the Southern Ocean State Estimate (SOSE), and (ii) leverage the capabilities of these APIs with interactive educational activities suitable for courses in oceanography and air-sea interactions. In addition, the notebook serves as a basis for research applications of the two APIs, e.g. to co-locate these datasets with oceanic observations (e.g. profiles of ocean temperature and salinity from Argo) for interdisciplinary research at the interface of different climate system components: the ocean, the atmosphere and the cryosphere.
