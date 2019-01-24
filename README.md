# pyOSOAA

`pyOSOAA` is a python interface for the Ocean Successive Orders with Atmosphere - Advanced (OSOAA) radiative transfer. The OSOAA is a radiative transfer code deveoped in the Sorbonne Université by Pr. Malik Chami group and supported by the CNES. 

The coded is based in the sucessive orders of scattering method and the OSOA code developed by Malik Chami in 2001 that included the computation of the radiane and polarization for the ocean-atmosphere system with a flat surface.

The OSOAA code simulates:

- **Atmospheric and sea profiles**: The atmosphere can be and sea profiles can be defined by the user both for the molecules and aerosol in the atmosphere and the water column, chlorophyll and mineral-like particles in the sea. Detritus and yellow substance absorption can also be modeled.
- **Aerosol models**: Aerosol models inclue WMO, LND, Junge mono-modal, bimodal LND and Shettle and Fenn.
- **Hydrosol models**: For phytoplankton and mineral-like particles including scattering and absorving properties.
- **Sea surface interface**: Both for a flat surface or by a roughsurface using Cox and Munk model.

The `pyOSOAA` interface aims to incorporate the creation of run scripts and parsing of output results for the OSOAA model. It also incorporates helpers to perform common task like calculating the radiance for a certain band instead of a wavelength or runnig the model for multiple wavelengths. 

This code was inspired by [py6S](https://github.com/robintw/Py6S) by Robin Wilson.


## Installation

The installation of the `pyOSOAA` has two parts.

First you need to install the OSOAA software package from https://github.com/CNES/RadiativeTransferCode-OSOAA.

Second you need to download the last version of the `pyOSOAA` from [github](https://github.com/fnemina/pyOSOAA/releases/latest).

Once downloaded decompress it, go to the folder containing the code and run

```
python setup.py install
```

To then check that software installed correctly

```
# Load pyOSOAA module
import pyOSOAA
# Run the test suite
pyOSOAA.test()
```
the following output should appear at the end of the screen
```
OSOAA wrapper script by Francisco Nemiña
Inspired by Py6S wrapper by Robin Wilson
Using OSOAA located at /home/fnemina/Documents/SABIA-Mar/RT/OSOAA_V1.5
Running OSOAA using a set of test parameters
The results are:
Expected result: 0.128266
Actual result: 0.128266
#### Results agree, PyOSOAA is working correctly
```
