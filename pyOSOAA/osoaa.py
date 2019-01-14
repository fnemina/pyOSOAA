import numpy as np
import os


class LND(object):
    """This is a lognormal distribution for the different models."""

    def __init__(self, mrwa, miwa, sdradius, sdvar, rate):
        """ Init function for the lognormal distribution.
            mrwa        Real part of the refractive index at simulation
                        wavelength
            miwa        Imaginary part of the refractive index at simulation
                        wavelength
            sdradius    Modal radius of the lognormal distribution
            sdvar       Standar deviation of the lognormal distribution
            rate        Ratio of the mode in the global distribution
            """

        self.mrwa = mrwa
        self.miwa = miwa
        self.sdradius = sdradius
        self.sdvar = sdvar
        self.rate = rate


class JD(object):
    """This is a Junge distribution for the different models."""

    def __init__(self, mrwa, miwa, slope, rmin, rmax, rate):
        """ Init function for the lognormal distribution.
            mrwa        Real part of the refractive index at simulation
                        wavelength
            miwa        Imaginary part of the refractive index at simulation
                        wavelength
            slope       Slope of Junge's law
            rmin        Minimal radius for Junge's law
            rmax        Maximal radius for Junge's law
            rate        Ratio of the mode in the global distribution
            """

        self.mrwa = mrwa
        self.miwa = miwa
        self.slope = slope
        self.rmin = rmin
        self.rmax = rmax
        self.rate = rate
