[![Build Status](https://travis-ci.org/fnemina/pyOSOAA.svg?branch=master)](https://travis-ci.org/fnemina/pyOSOAA) [![Coverage Status](https://coveralls.io/repos/github/fnemina/pyOSOAA/badge.svg?branch=master)](https://coveralls.io/github/fnemina/pyOSOAA?branch=master)

# pyOSOAA

`pyOSOAA` is a python interface for the Ocean Successive Orders with Atmosphere - Advanced (OSOAA) radiative transfer. The OSOAA is a radiative transfer code developed in the Sorbonne Université by Pr. Malik Chami group and supported by the CNES. 

The coded is based in the successive orders of scattering method and the OSOA code developed by Malik Chami in 2001 that included the computation of the radiance and polarization for the ocean-atmosphere system with a flat surface.

The OSOAA code simulates:

- **Atmospheric and sea profiles**: The atmosphere can be and sea profiles can be defined by the user both for the molecules and aerosol in the atmosphere and the water column, chlorophyll and mineral-like particles in the sea. Detritus and yellow substance absorption can also be modelled.
- **Aerosol models**: Aerosol models include WMO, LND, Junge mono-modal, bimodal LND and Shettle and Fenn.
- **Hydrosol models**: For phytoplankton and mineral-like particles including scattering and absorbing properties.
- **Sea surface interface**: Both for a flat surface or by a rough surface using Cox and Munk model.

The `pyOSOAA` interface aims to incorporate the creation of run scripts and parsing of output results for the OSOAA model. It also incorporates helpers to perform common tasks like calculating the radiance for a certain band instead of a wavelength or running the model for multiple wavelengths. 

This code was inspired by [py6S](https://github.com/robintw/Py6S) by Robin Wilson.

You can find the full `pyOSOAA` manual [here](https://pyosoaa.readthedocs.io/en/latest/).

## Installation

The installation of the `pyOSOAA` has two parts.

First, you need to install the OSOAA software package from https://github.com/CNES/RadiativeTransferCode-OSOAA.

Second, install pyOSOAA. There are two ways to install pyOSOAA.

### Install pyOSOAA from pypi

```bash
pip install pyOSOAA
```

### Install pyOSOAA from source code
Download the last version of the `pyOSOAA` from [github](https://github.com/fnemina/pyOSOAA/releases/latest).

Once downloaded decompress it, go to the folder containing the code and run

```bash
pip install .
```

To then check that software installed correctly

```python
# Load pyOSOAA module
import pyOSOAA
# Run the test suite
pyOSOAA.test()
```
the following output should appear at the end of the screen
```
OSOAA wrapper script by Francisco Nemiña
Inspired by Py6S wrapper by Robin Wilson
Using OSOAA located at /home/.../OSOAA_V1.5
Running OSOAA using a set of test parameters
The results are:
Expected result: 0.128266
Actual result: 0.128266
#### Results agree PyOSOAA is working correctly
```
