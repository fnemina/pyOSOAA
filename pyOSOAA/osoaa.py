import numpy as np
import os

class LND(object):
    """ This is a lognormal distribution for the different models. """

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
    """ This is a Junge distribution for the different models. """

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


class SEA(object):
    """ This is the SEA class which defines the interfaces at the bottom of the
        ocean and the interface with the air. """

    def __init__(self, surfalb, bottype, botalb, log, wind, dir, ind, depth):
        """ Init function for the SEA class.
            surfalb     Foam lambertian reflectance for the current wavelength.
            bottype     Type of seabottom for albedo definition.
            botalb      Seabottom albedo for botype = 1.
            log         Log filename for surface computations
            wind        Wind speed at seasurface (m/s)
            dir         Directory for surface files.
            ind         Surface/atmosphere refraction index, air = 1.
            depth       Sea depth in meters.
            """

        self.surfalb = surfalb
        self.bottype = bottype
        self.botalb = botalb
        self.log = log
        self.wind = wind
        self.dir = dir
        self.ind = ind
        self.depth = depth


class LOG(object):
    """ Log files for the simulation.
    osoaa       The main logfile gives information on the different routines
                execution (input/output parameters, warnings, error cases, ...)
    ang         This logfile gives information about the angles computations:
                Gauss angles used for the phase functions and radiance
                calculations, solar zenith angle in the air and sea, user’s
                angles.
    profile     This logfile gives information about the atmospheric and marine
                profiles computation.
    aer         This logfile gives information about the radiative properties
                computation for aerosols: phase function of each elementary
                components, mixed-average phase function, truncation
                calculations, scattering and extinction cross-sections.
    aermie      This logfile gives the value of matrix phase function
                (polarized form) versus the size parameter from Mie
                calculations.
    hyd         As for aerosols, this logfile gives information about the
                radiative properties computation for hydrosols.
    hydmie      As for aerosols, this logfile gives information about the Mie
                computations for hydrosols.
    sea         This logfile gives information about the calculation of
                interface matrices (RAA, RWW, TAW, TWA): Fresnel matrices,
                probability function of the waves orientation, etc.
    sos         This logfile gives information about the successive ord
    """
    osoaa = "log_osoaa.txt"
    ang = "log_ang.txt"
    profile = "log_profile.txt"
    aer = "log_aer.txt"
    aermie = "log_aermie.txt"
    hyd = "log_hyd.txt"
    hydmie = "log_aermie.txt"
    sea = "log_sea.txt"
    sos = "log_sos.txt"


class RESULTS(object):
    """ Result files for the simulation."""
    profileatm = None
    profilesea = None
    aer = None
    pytho = None
    mlp = None
    angrad = None
    angmie = None
    sosbin = None
    vsvza = None
    advup = "resfile_advup.txt"
    advdown = "resfile_advdown.txt"
    vsz = "resfile_vsz.txt"
    viewvza = "resfile_viewvza.txt"


class DIRMIE(object):
    """ Directory for hydrosol MIE files storage (complete path)
    aer     aerosols mie storage
    hid     hidrosols mie storage
    """
    hid = None
    aer = None
