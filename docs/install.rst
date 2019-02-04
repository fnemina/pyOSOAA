------------
Installation
------------

This document describes how to install the ``OSOAA`` and ``pyOSOAA`` software for Ubuntu 18.04 or greater with python 3.6 or greater.

------------------
OSOAA installation
------------------

To install ``OSOAA`` we first need to set up some previous software.

1. We begin by installing the ``gfortran``  Fortran compiler and the Korn shell from the terminal by running ::
  $ sudo apt update
  
  $ sudo apt install gfortran ksh


2. We then clone the ``OSOAA`` `repository<https://github.com/CNES/RadiativeTransferCode-OSOAA>`_::
  $ git clone https://github.com/CNES/RadiativeTransferCode-OSOAA.git


3. Once this is done we add the ``OSOAA_ROOT`` variable to the bash path by doing and sourcing it::
  $ echo "export OSOAA_ROOT="FULL PATH TO OSOAA FOLDER" >> ~/.bashrc
  
  $ source ~/.bashrc


4. We then create a folder to contain the object files in the OSOAA folder::
  $ bash
  
  $ mkdir $OSOAA_ROOT/obj


5. Finally, we compile the ``OSOAA``::
  $ make -f $OSOAA_ROOT/gen/Makefile_OSOAA.gfortran

--------------------
pyOSOAA installation
--------------------

1. We first download and install `miniconda<https://conda.io/en/latest/miniconda.html>`_  with python 3.7 by doing::
  $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
  
  $ bash Miniconda3-latest-Linux-x86_64.sh


and following the onscreen instructions.

2. We then install the necessary python libraries::
  $ conda install numpy scipy matplotlib ipython


3. Once python is installed and configured we clone or `download<https://github.com/fnemina/pyOSOAA/releases/latest>`_ the lastest `pyOSOAA` version::
  $ git clone https://github.com/fnemina/pyOSOAA.git


4. We then go into the `pyOSOAA` folder and install it by running::
  $ python setup.py install

