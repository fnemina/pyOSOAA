import numpy as np
import os


class LND(object):
    """This is a lognormal distribution for the different models."""

    def __init__(self, MRwa, MIwa, SDradius, SDvar, rate):
        """ Init function for the lognormal distribution.
            MRwa        Real part of the refractive index at simulation
                        wavelength
            MIwa        Imaginary part of the refractive index at simulation
                        wavelength
            SDradius    Modal radius of the lognormal distribution
            SDvar       Standar deviation of the lognormal distribution
            rate        Ratio of the mode in the global distribution
            """

        self.MRwa = MRwa
        self.MIwa = MIwa
        self.SDradius = SDradius
        self.SDvar = SDvar
        self.rate = rate


class JD(object):
    """This is a Junge distribution for the different models."""

    def __init__(self, MRwa, MIwa, slope, rmin, rmax, rate):
        """ Init function for the lognormal distribution.
            MRwa        Real part of the refractive index at simulation
                        wavelength
            MIwa        Imaginary part of the refractive index at simulation
                        wavelength
            slope       Slope of Junge's law
            rmin        Minimal radius for Junge's law
            rmax        Maximal radius for Junge's law
            rate        Ratio of the mode in the global distribution
            """

        self.MRwa = MRwa
        self.MIwa = MIwa
        self.slope = slope
        self.rmin = rmin
        self.rmax = rmax
        self.rate = rate
