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
        """ Init function for the Junges distribution.
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


class GP(object):
    """ Gaussian profile class definition."""

    def __init__(self, chlbg, deep, width):
        """ Init function for the gaussian profiles
        chlbg       Constant background biomass (mg/m3)
        deep        Maximun deep of the gaussian profile (m)
        width       Peak width of the gaussian profile (m)
        """

        self.chlbg = chlbg
        self.deep = deep
        self.width = width


class PHYTO(object):
    """ Phytoplanckton class"""
    # Profile types
    Homogeneous = 1
    Gaussian = 2
    UserDefined = 3

    def __init__(self, chl=0):
        """ Init function for the Phytoplanckton profiles
        chl     Phytoplanckton concentration
        """

        self.chl = chl
        self.profiltype = self.Homogeneous
        self.jd = None
        self.sm = None
        self.tm = None

    def SetPrimaryMode(self, mrwa, miwa, slope, rmin, rmax, rate):
        """ Sets the primary mode using Junge's law
            mrwa        Real part of the refractive index at simulation
                        wavelength
            miwa        Imaginary part of the refractive index at simulation
                        wavelength
            slope       Slope of Junge's law
            rmin        Minimal radius for Junge's law
            rmax        Maximal radius for Junge's law
            rate        Ratio of the mode in the global distribution
            """

        self.jd = JD(mrwa, miwa, slope, rmin, rmax, rate)

    def SetSecondaryMode(self, mrwa, miwa, sdradius, sdvar, rate):
        """ Sets the secondary mode using lognormal distribution
            mrwa        Real part of the refractive index at simulation
                        wavelength
            miwa        Imaginary part of the refractive index at simulation
                        wavelength
            sdradius    Modal radius of the lognormal distribution
            sdvar       Standar deviation of the lognormal distribution
            rate        Ratio of the mode in the global distribution
            """

        self.sm = LND(mrwa, miwa, sdradius, sdvar, rate)

    def SetTertiaryMode(self, mrwa, miwa, sdradius, sdvar, rate):
        """ Sets the secondary mode using lognormal distribution
            mrwa        Real part of the refractive index at simulation
                        wavelength
            miwa        Imaginary part of the refractive index at simulation
                        wavelength
            sdradius    Modal radius of the lognormal distribution
            sdvar       Standar deviation of the lognormal distribution
            rate        Ratio of the mode in the global distribution
            """

        self.tm = LND(mrwa, miwa, sdradius, sdvar, rate)

    def SetProfilType(self, profiltype, chlbg, deep, width, userfile):
        """ This method sets the profile type for the Phytoplanckton
            distribution. This also configures the parameters for each
            profile type

            profiltype      Profile for the chlorophyl distribution.
                            1 - homogeneous profile
                            2 - Gaussian profile
                            3 - User defined profile

            for the Gaussian profiltype
            chlbg           Constant background biomass (mg/m3)
            deep            Maximun deep of the gaussian profile (m)
            width           Peak width of the gaussian profile (m)

            for the User defined profile
            userfile        Userfile describing the chlorophyll profile
            """

        if profiltype not in [self.Homogeneous,
                              self.Gaussian,
                              self.UserDefined]:
            raise Exception("Invalid profile type.")

        self.profiltype = profiltype

        if profiltype == self.Gaussian:
            # We confifigure the gaussian profile
            self.gp = GP(chlbg, deep, width)

        elif profiltype == self.UserDefined:
            # We configure the user defined profile
            try:
                open(userfile, 'r')
                self.usefile = userfile
            except FileNotFoundError:
                print("File {} not found".format(userfile))


class SED(object):
    """ Sediment class"""

    def __init__(self, csed=0):
        """ Init function for the sediment profiles
        csed     mineral like particles concentration mg/l
        """

        self.csed = csed
        self.jd = None
        self.sm = None
        self.tm = None

    def SetPrimaryMode(self, mrwa, miwa, slope, rmin, rmax, rate):
        """ Sets the primary mode using Junge's law
            mrwa        Real part of the refractive index at simulation
                        wavelength
            miwa        Imaginary part of the refractive index at simulation
                        wavelength
            slope       Slope of Junge's law
            rmin        Minimal radius for Junge's law
            rmax        Maximal radius for Junge's law
            rate        Ratio of the mode in the global distribution
            """

        self.jd = JD(mrwa, miwa, slope, rmin, rmax, rate)

    def SetSecondaryMode(self, mrwa, miwa, sdradius, sdvar, rate):
        """ Sets the secondary mode using lognormal distribution
            mrwa        Real part of the refractive index at simulation
                        wavelength
            miwa        Imaginary part of the refractive index at simulation
                        wavelength
            sdradius    Modal radius of the lognormal distribution
            sdvar       Standar deviation of the lognormal distribution
            rate        Ratio of the mode in the global distribution
            """

        self.sm = LND(mrwa, miwa, sdradius, sdvar, rate)

    def SetTertiaryMode(self, mrwa, miwa, sdradius, sdvar, rate):
        """ Sets the secondary mode using lognormal distribution
            mrwa        Real part of the refractive index at simulation
                        wavelength
            miwa        Imaginary part of the refractive index at simulation
                        wavelength
            sdradius    Modal radius of the lognormal distribution
            sdvar       Standar deviation of the lognormal distribution
            rate        Ratio of the mode in the global distribution
            """

        self.tm = LND(mrwa, miwa, sdradius, sdvar, rate)


class ABS(object):
    """ Absorption class to be used with yellow substance and detritus"""

    def __init__(self, abs440, swa):
        """ Init function for the absorption class.ABS
            abs440      Absorption coefficient 1/m
            swa         spectral variation coefficient
        """

        self.abs440 = abs440

        if abs440 == 0:
            self.swa = None
        elif abs440 > 0:
            self.swa = swa
        else:
            raise Exception("Invalid absorption value.")
