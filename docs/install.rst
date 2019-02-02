---
layout: default
title: Installation
navigation_weight: 2
---

# Installation

This document describes how to install the `OSOAA` and `pyOSOAA` software for Ubuntu 18.04 or greater with python 3.6 or greater.


## OSOAA installation

To install `OSOAA` we first need to set up some previous software.

1. We begin by installing the `gfortran`  Fortran compiler and the Korn shell from the terminal by running

```bash
sudo apt update
sudo apt install gfortran k
```

2. We then clone the `OSOAA` [repository](https://github.com/CNES/RadiativeTransferCode-OSOAA)

```bash
git clone https://github.com/CNES/RadiativeTransferCode-OSOAA.git
```

3. Once this is done we add the `OSOAA_ROOT` variable to the bash path by doing and sourcing it

```bash
echo "export OSOAA_ROOT="FULL PATH TO OSOAA FOLDER" >> ~/.bashrc
source ~/.bashrc
```

4. We then create a folder to contain the object files in the OSOAA folder

```bash
mkdir $OSOAA_ROOT/obj
```

5. Finally, we compile the `OSOAA`

```bash
make -f $OSOAA_ROOT/gen/Makefile_OSOAA.gfortran
```

if the installation is correct the following should appear in the screen
```
------------> Link of $OSOAA_ROOT/exe/OSOAA_MAIN.exe
------------>ok
```

## pyOSOAA installation

1. We first download and install [miniconda](https://conda.io/en/latest/miniconda.html)  with python 3.7 by doing
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
and following the onscreen instructions.

2. We then install the necessary python libraries
```bash
conda install numpy scipy matplotlib ipython
```

3. Once python is installed and configured we clone or [download](https://github.com/fnemina/pyOSOAA/releases/latest) the lastest `pyOSOAA` version
```bash
git clone https://github.com/fnemina/pyOSOAA.git
```

4. We then go into the `pyOSOAA` folder and install it by running
```bash
python setup.py install
```

## Installation testing

To test the installation we open a python session and run

```python
import pyOSOAA
```
and then
```python
pyOSOAA.test()
```
the following output should appear on screen

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
