import os
import random
import string
from io import open
from .outputs import OUTPUTS


class SEA(object):
    """ This is the SEA class which defines the interfaces at the bottom of the
        ocean and the interface with the air. """

    def __init__(self, surfalb=0.0, bottype=1, botalb=0.30, wind=7,
                 ind=1.34, depth=15.0):
        """ Init function for the SEA class.
            surfalb     Foam lambertian reflectance for the wavelength of
                        radiance calculation (i.e., albedo of the foam at the
                        sea surface).
            bottype     Type of sea bottom for albedo definition
                        Cases :	1 : Lambertian value (user data SEA.botalb)
                                2 : Light sand	     (tabulated data)
                                3 : Green algua	     (tabulated data)
                                4 : Brown algae	     (tabulated data)
                                5 : Red algae	     (tabulated data)
            botalb      Sea bottom albedo for the wavelength of radiance
                        calculation (lambertian component).
            wind        Wind velocity (m/s).
            ind         Surface / atmosphere refractive index (air = 1).
            depth       Sea depth value  (meters).
                        (if None, euphotic depth will be used from Morel
                        tabulated data with regards to the chlorophyll
                        concentration at sea surface)
            """

        self.surfalb = surfalb
        self.bottype = bottype
        self.botalb = botalb
        self.wind = wind
        self.ind = ind
        self.depth = depth


class LOG(object):
    """ Log files for the simulation.
    osoaa       log filename for ANGLES computations
                (defined without directory tree ==> this file will be located
                in the sub-directory Advanced_outputs of the working folder).
                Only created if the log filename is specified.
    ang         log filename for ANGLES computations (defined without directory
                tree ==> this file will be located in the sub-directory
                Advanced_outputs of the working folder).
                Only created if the log filename is specified.
    profile     log filename for PROFILE computations
                (defined without directory tree ==> this file will be located
                in the sub-directory Advanced_outputs of the working folder).
                Only created if the log filename is specified.
    aer         log filename for AEROSOLS file computations
                (defined without directory tree ==> this file will be located
                in the sub-directory Advanced_outputs of the working folder).
                Only created if the log filename is specified.
    aermie      log filename for MIE aerosols computations
                (defined without directory tree ==> this file will be located
                in the sub-directory Advanced_outputs of the working folder).
                Only created if the log filename is specified.
    hyd         log filename for HYDROSOLS file computations
                (defined without directory tree ==> this file will be located
                in the sub-directory Advanced_outputs of the working folder).
                Only created if the log filename is specified.
    hydmie      log filename for MIE hydrosols computations
                (defined without directory tree ==> this file will be located
                in the sub-directory Advanced_outputs of the working folder).
                Only created if the log filename is specified.
    sea         log filename for SURFACE file computations (defined without
                directory tree ==> this file will be located in the
                sub-directory Advanced_outputs of the working folder).
                Only created if the log filename is specified.
    sos         log filename for SOS computations (defined without directory
                tree ==> this file will be located in the sub-directory
                Advanced_outputs of the working folder).
                Only created if the log filename is specified.
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
    """ Result files for the simulation.
        profileatm  Filename of the result for atmospheric PROFILE computations
                    (defined without directory tree ==> this file will be
                    located in the sub-directory Advanced_outputs of the
                    working folder).
        profilesea  Filename of the result for the marine PROFILE computations
                    (defined without directory tree ==> this file will be
                    located in the sub-directory Advanced_outputs of the
                    working folder).
        aer         Filename of the result OSOAA_AEROSOLS computations
                    (defined without directory tree ==> this file will be
                    located in the sub-directory Advanced_outputs of the
                    working folder).
        phyto       Filename of the result OSOAA_HYDROSOLS computations or
                    phytoplankton particles (defined without directory tree ==>
                    this file will be located in the sub-directory
                    Advanced_outputs of the working folder).
                    NB : Include the result of phase matrix development from
                    user file global phase function in case of an user file is
                    used (-HYD.ExtData).
        mlp         Filename of the result OSOAA_HYDROSOLS computations for
                    Mineral Like particles (defined without directory tree ==>
                    this file will be located in the sub-directory
                    Advanced_outputs of the working folder).
                    NB : Nul coefficients of phase matrix development
                    in case of an user file is used (-HYD.ExtData).
        angrad      Filename of list of angles used to BRDF/BPDF and radiance
                    computations. (defined without directory tree ==> this file
                    will be located in the sub-directory Advanced_outputs of
                    the working folder).
        angmie      Filename of list of angles used to the matrix phase
                    function computations. (defined without directory tree ==>
                    this file will be located in the sub-directory
                    Advanced_outputs of the working folder).
        sosbin      Filename of the binary file resulting from SOS computations
                    (defined without directory tree ==> this file will be
                    located in the sub-directory Advanced_outputs of the
                    working folder).
        vsvza       Filename of the ascii file resulting from SOS computations
                    (defined without directory tree ==> this file will be
                    located in the sub-directory Standard_outputs of the
                    working folder) : ==> Output radiance field versus the
                    viewing zenith angle (for the given relative azimuth angle
                    and given altitude or depth).
        advup       Filename of the ascii file resulting from SOS computations
                    (defined without directory tree ==> this file will be
                    located in the sub-directory Advanced_outputs of the
                    working folder) : ==> Advanced output upward radiance field
                    versus the depth (or altitude) AND versus the viewing
                    zenith angle (for the given relative azimuth angle).
        advdown     Filename of the ascii file resulting from SOS computations
                    (defined without directory tree ==> this file will be
                    located in the sub-directory Advanced_outputs of the
                    working folder) : ==> Advanced output upward radiance field
                    versus the depth (or altitude) AND versus the viewing
                    zenith angle (for the given relative azimuth angle).
        vsz         Filename of the ascii file resulting from SOS computations
                    (defined without directory tree ==> this file will be
                    located in the sub-directory Standard_outputs of the
                    working folder) : ==> Output radiance field versus the
                    depth (or altitude) (for the given relative azimuth angle
                    and given viewing zenith angle).
        """
    profileatm = None
    profilesea = None
    aer = None
    phyto = None
    mlp = None
    angrad = None
    angmie = None
    sosbin = None
    vsvza = None
    advup = "resfile_advup.txt"
    advdown = "resfile_advdown.txt"
    vsz = "resfile_vsz.txt"


class DIRMIE(object):
    """ Directory for hydrosol MIE files storage
    """

    def __init__(self, osoaaroot, hid="/DATABASE/MIE_HYD",
                 aer="/DATABASE/MIE_AER", sea="/DATABASE/SURF_MATR"):
        """Directory for hydrosol MIE files storage (complete path)
            aer         Storage directory for MIE files producted by
                        OSOAA_AEROSOLS computations (complete path).
            hid         Storage directory for MIE files producted by
                        HYDROSOLS_AEROSOLS computations (complete path).
            SEA         Directory for SURFACE files storage (complete path).
            """
        self.hid = osoaaroot+hid
        self.aer = osoaaroot+aer
        self.sea = osoaaroot+sea


class GP(object):
    """ Gaussian profile class definition."""

    def __init__(self, chlbg, deep, width):
        """ Init function for the gaussian profiles
        chlbg       Constant biomass background (mg/m3)
        deep        Maximum deep of the gaussian chlorophyll profile (m)
        width       Peak width of the gaussian chlorophyll profile (m)
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

    class JD(object):
        """ This is a Junge distribution for the different models. """

        def __init__(self, mrwa=1.05, miwa=-0.0, slope=4.0,
                     rmin=0.01, rmax=200, rate=1.0):
            """ Init function for the Junges distribution.
                mrwa        Real part of the refractive index for phytoplankton
                            particles at the simulation wavelength: main mode
                            of particles (Junge distribution)
                miwa        Imaginary part of the refractive index for
                            phytoplankton particles at the simulation
                            wavelength: main mode of particles (Junge
                            distribution)
                slope       Slope of Junge's law for phytoplankton particles
                rmin        Minimal radius of Junge's law for phytoplankton
                            particles (microns)
                rmax        Maximal radius of Junge's law for phytoplankton
                            particles (microns)
                rate        Ratio of the Junge's law in the global distribution
                            for phytoplankton particles  ==> as a proportion of
                            the Junge distribution particles versus the global
                            amount of phytoplankton particles.
                """

            self.mrwa = mrwa
            self.miwa = miwa
            self.slope = slope
            self.rmin = rmin
            self.rmax = rmax
            self.rate = rate

    class LND(object):
        """ This is a lognormal distribution for the different models. """

        def __init__(self, mrwa, miwa, sdradius, sdvar, rate):
            """ Init function for the lognormal distribution.
                mrwa        Real part of the refractive index for phytoplankton
                            particles at the simulation wavelength (LND
                            distribution)
                miwa        Imaginary part of the refractive index for
                            phytoplankton particles at the simulation
                            wavelength (LND distribution)
                sdradius    Modal radius of the LND of phytoplankton particles
                            (mic)
                sdvar       Standard deviation of the LND of phytoplankton
                            particle
                rate        Ratio of the LND mode in the global distribution
                            for phytoplankton particles ==> as a proportion of
                            the number of LND particles versus the global
                            amount of phytoplankton particles.
                """

            self.mrwa = mrwa
            self.miwa = miwa
            self.sdradius = sdradius
            self.sdvar = sdvar
            self.rate = rate

    def __init__(self, chl=0.2):
        """ Init function for the Phytoplanckton profiles
        chl         chlorophyll concentration at sea surface (mg/m3)
        """

        self.chl = chl
        self.profiltype = self.Homogeneous
        self.jd = self.JD()
        self.sm = None
        self.tm = None

    def SetPrimaryMode(self, mrwa=1.05, miwa=-0.0, slope=4.0,
                       rmin=0.01, rmax=200, rate=1.0):
        """ Sets the primary mode using Junge's law
                mrwa        Real part of the refractive index for phytoplankton
                            particles at the simulation wavelength: main mode
                            of particles (Junge distribution)
                miwa        Imaginary part of the refractive index for
                            phytoplankton particles at the simulation
                            wavelength: main mode of particles (Junge
                            distribution)
                slope       Slope of Junge's law for phytoplankton particles
                rmin        Minimal radius of Junge's law for phytoplankton
                            particles (microns)
                rmax        Maximal radius of Junge's law for phytoplankton
                            particles (microns)
                rate        Ratio of the Junge's law in the global distribution
                            for phytoplankton particles  ==> as a proportion of
                            the Junge distribution particles versus the global
                            amount of phytoplankton particles.
            """

        self.jd = self.JD(mrwa, miwa, slope, rmin, rmax, rate)

    def SetSecondaryMode(self, mrwa, miwa, sdradius, sdvar, rate):
        """ Sets the secondary mode using lognormal distribution
                mrwa        Real part of the refractive index for phytoplankton
                            particles at the simulation wavelength (LND
                            distribution)
                miwa        Imaginary part of the refractive index for
                            phytoplankton particles at the simulation
                            wavelength (LND distribution)
                sdradius    Modal radius of the LND of phytoplankton particles
                            (mic)
                sdvar       Standard deviation of the LND of phytoplankton
                            particle
                rate        Ratio of the LND mode in the global distribution
                            for phytoplankton particles ==> as a proportion of
                            the number of LND particles versus the global
                            amount of phytoplankton particles.
            """

        self.sm = self.LND(mrwa, miwa, sdradius, sdvar, rate)

    def SetTertiaryMode(self, mrwa, miwa, sdradius, sdvar, rate):
        """ Sets the secondary mode using lognormal distribution
                mrwa        Real part of the refractive index for phytoplankton
                            particles at the simulation wavelength (LND
                            distribution)
                miwa        Imaginary part of the refractive index for
                            phytoplankton particles at the simulation
                            wavelength (LND distribution)
                sdradius    Modal radius of the LND of phytoplankton particles
                            (mic)
                sdvar       Standard deviation of the LND of phytoplankton
                            particle
                rate        Ratio of the LND mode in the global distribution
                            for phytoplankton particles ==> as a proportion of
                            the number of LND particles versus the global
                            amount of phytoplankton particles.
            """

        self.tm = self.LND(mrwa, miwa, sdradius, sdvar, rate)

    def SetProfilType(self, profiltype, chlbg, deep, width, userfile):
        """ This method sets the profile type for the Phytoplanckton
            distribution. This also configures the parameters for each
            profile type

            profiltype      Profile for the chlorophyl distribution.
                            1 - homogeneous profile
                            2 - Gaussian profile
                            3 - User defined profile

            for the Gaussian profiltype
            chlbg       Constant biomass background (mg/m3)
            deep        Maximum deep of the gaussian chlorophyll profile (m)
            width       Peak width of the gaussian chlorophyll profile (m)

            for the User defined profiltype
            userfile        Name of user phytoplankton profile file
                            (complete access)
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

    class JD(object):
        """ This is a Junge distribution for the different models. """

        def __init__(self, mrwa, miwa, slope, rmin, rmax, rate):
            """ Init function for the Junges distribution.
                mrwa        Real part of the refractive index for mineral-like
                            particles at the simulation wavelength: main mode
                            of particles (Junge distribution)
                miwa        Imaginary part of the refractive index for
                            mineral-like particles at the simulation
                            wavelength: main mode of particles (Junge
                            distribution)
                slope       Slope of Junge's law for mineral-like particles
                rmin        Minimal radius of Junge's law for phytoplankton
                            particles (microns)
                rmax        Maximal radius of Junge's law for mineral-like
                            particles (microns)
                rate        Ratio of the Junge's law in the global distribution
                            for mineral-like particles  ==> as a proportion of
                            the Junge distribution particles versus the global
                            amount of mineral-like particles.
                """

            self.mrwa = mrwa
            self.miwa = miwa
            self.slope = slope
            self.rmin = rmin
            self.rmax = rmax
            self.rate = rate

    class LND(object):
        """ This is a lognormal distribution for the different models. """

        def __init__(self, mrwa, miwa, sdradius, sdvar, rate):
            """ Init function for the lognormal distribution.
                mrwa        Real part of the refractive index for mineral-like
                            particles at the simulation wavelength (LND
                            distribution)
                miwa        Imaginary part of the refractive index for
                            mineral-like particles at the simulation
                            wavelength (LND distribution)
                sdradius    Modal radius of the LND of mineral-like particles
                            (mic)
                sdvar       Standard deviation of the LND of mineral-like
                            particle
                rate        Ratio of the LND mode in the global distribution
                            for mineral-like particles ==> as a proportion of
                            the number of LND particles versus the global
                            amount of mineral-like particles.
                """

            self.mrwa = mrwa
            self.miwa = miwa
            self.sdradius = sdradius
            self.sdvar = sdvar
            self.rate = rate

    def __init__(self, csed=0.0):
        """ Init function for the sediment profiles
        csed    Concentration of mineral-like particles at sea surface
                (mg/liter)
        """

        self.csed = csed
        self.jd = None
        self.sm = None
        self.tm = None

    def SetPrimaryMode(self, mrwa, miwa, slope, rmin, rmax, rate):
        """ Sets the primary mode using Junge's law
                mrwa        Real part of the refractive index for mineral-like
                            particles at the simulation wavelength: main mode
                            of particles (Junge distribution)
                miwa        Imaginary part of the refractive index for
                            mineral-like particles at the simulation
                            wavelength: main mode of particles (Junge
                            distribution)
                slope       Slope of Junge's law for mineral-like particles
                rmin        Minimal radius of Junge's law for phytoplankton
                            particles (microns)
                rmax        Maximal radius of Junge's law for mineral-like
                            particles (microns)
                rate        Ratio of the Junge's law in the global distribution
                            for mineral-like particles  ==> as a proportion of
                            the Junge distribution particles versus the global
                            amount of mineral-like particles.
            """

        self.jd = self.JD(mrwa, miwa, slope, rmin, rmax, rate)

    def SetSecondaryMode(self, mrwa, miwa, sdradius, sdvar, rate):
        """ Sets the secondary mode using lognormal distribution
                mrwa        Real part of the refractive index for mineral-like
                            particles at the simulation wavelength (LND
                            distribution)
                miwa        Imaginary part of the refractive index for
                            mineral-like particles at the simulation
                            wavelength (LND distribution)
                sdradius    Modal radius of the LND of mineral-like particles
                            (mic)
                sdvar       Standard deviation of the LND of mineral-like
                            particle
                rate        Ratio of the LND mode in the global distribution
                            for mineral-like particles ==> as a proportion of
                            the number of LND particles versus the global
                            amount of mineral-like particles.
            """

        self.sm = self.LND(mrwa, miwa, sdradius, sdvar, rate)

    def SetTertiaryMode(self, mrwa, miwa, sdradius, sdvar, rate):
        """ Sets the secondary mode using lognormal distribution
                mrwa        Real part of the refractive index for mineral-like
                            particles at the simulation wavelength (LND
                            distribution)
                miwa        Imaginary part of the refractive index for
                            mineral-like particles at the simulation
                            wavelength (LND distribution)
                sdradius    Modal radius of the LND of mineral-like particles
                            (mic)
                sdvar       Standard deviation of the LND of mineral-like
                            particle
                rate        Ratio of the LND mode in the global distribution
                            for mineral-like particles ==> as a proportion of
                            the number of LND particles versus the global
                            amount of mineral-like particles.
            """

        self.tm = self.LND(mrwa, miwa, sdradius, sdvar, rate)


class YS(object):
    """ Absorption class to be used with yellow substance"""

    def __init__(self, abs440=0.0, swa=None):
        """ Init function for the absorption class.ABS
            abs440      Yellow substance abs. coef. (m-1)  at 440 nm
            swa         Coefficient of spectral variation for yellow substance
                        absorption (m-1)
        """

        self.abs440 = abs440

        if abs440 == 0:
            self.swa = None
        elif abs440 > 0:
            self.swa = swa
        else:
            raise Exception("Invalid absorption value.")


class DET(object):
    """ Absorption class to be used with detritus"""

    def __init__(self, abs440=0.0, swa=None):
        """ Init function for the absorption class.ABS
            abs440      Detritus absorption coef. (m-1)  at 440 nm
            swa         sCoefficient of spectral variation for detritus
                        absorption (m-1)
        """

        self.abs440 = abs440

        if abs440 == 0:
            self.swa = None
        elif abs440 > 0:
            self.swa = swa
        else:
            raise Exception("Invalid absorption value.")


class AP(object):
    """ Atmospheric profile parameters object."""

    def __init__(self, mot=None, pressure=1013.00, hr=8.0, ha=2.0):
        """ Init function for the atmospheric profile
            mot         Molecular optical thickness for the wavelength of
                        radiance simulation
            pressure    Atmospheric pressure at sea level (mbar)
            hr          Molecular heigth scale (km).
            ha          Aerosol heigth scale (km).
            """

        self.mot = mot
        self.pressure = pressure
        self.hr = hr
        self.ha = ha

        def SetPressure(self, pressure=1013.25):
            """ Define the molecular optical thickness with pressure
                pressure    atmospheric pressure at sea level (mbar)"""
            self.mOT = None
            self.pressure = pressure

        def SetMot(self, mot=0.1, hr=8):
            """ Define the molecular optical thickness
                mot         molecular optical thicknes
                hr          height scale for molecules (km)"""
            self.mot = mot
            self.hr = hr
            self.pressure = None


class AEROSOLMODELS(object):
    """ Aerosol models for the AER class."""
    class SF(object):
        """ Shettle and Fenn atmosphere model class."""
        # Shettle and Fenn models
        Tropospheric = 1
        Urban = 2
        Maritime = 3
        Coastal = 4

        def __init__(self, model=3, rh=98):
            """ Init method for the Shettle-Fenn model.
                model       Type of Shettle & Fenn model.
                                1 : Tropospheric S&F model.
                                2 : Urban S&F model.
                                3 : Maritime S&F model.
                                4 : Coastal S&F model.
                rh          Relative humidity (%) for Shettle & Fenn model.
                """
            self.model = model
            self.rh = rh


class AER(object):
    """ This class contains everything related to the aerosol components
        of the atmosphere."""

    def __init__(self, waref=0.550, aotref=0.1, tronca=None, model=2):
        """ Init method for the aerosol componentes class
            waref       Wavelength (microns) for reference aerosol optical
                        thickness.
            aotref      Aerosol optical thickness for the reference wavelength.
                        --> real value, without applied truncation.
            tronca      Option for no aerosol phase function troncation
                        (0 to not apply a troncation). Default is 1.
            model       Type of aerosol model
                            0 : Mono-modal
                            1 : WMO multi-modal
                            2 : Shettle & Fenn bi-modal
                            3 : Log-Normal bi-modal
                            4 : Phase function from an external source
            mm          Mono-modal model parameters
            sf          Shettle and Fenn model parameters
            wmo         WMO model parameters
            lnd         Log-Normal bi-modal model parameters
            external    External phase function
            """

        self.waref = waref
        self.aotref = aotref
        self.tronca = tronca
        self.model = model
        self.sf = AEROSOLMODELS.SF()

    def SetModel(self, model=2):
        """ This methods sets the model for the AER class."""
        self.model = model
        if model is 0:
            self.mm = None
            self.sf = None
            self.wmo = None
            self.lnd = None
            self.external = None
        elif model is 1:
            self.mm = None
            self.sf = None
            self.wmo = None
            self.lnd = None
            self.external = None
        elif model is 2:
            self.mm = None
            self.sf = AEROSOLMODELS.SF()
            self.wmo = None
            self.lnd = None
            self.external = None
        elif model is 3:
            self.mm = None
            self.sf = None
            self.wmo = None
            self.lnd = None
            self.external = None
        elif model is 4:
            self.mm = None
            self.sf = None
            self.wmo = None
            self.lnd = None
            self.external = None


class HYD(object):
    """ This class contains everything related to the hydrosol components
        of the sea."""

    def __init__(self, model=1, extdata=None):
        """ Init method for the aerosol componentes class
            model       Type of hydrosol characterization
                            1 : From size distribution models.
                            2 : Use of external phase functions
            extdata     Filename (complete path) of user's external phase
                        functions and radiative parameters (extinction and
                        scattering coefficients)
            """

        self.model = model
        self.extdata = extdata


class ANG(object):
    """ Angle definitions class."""

    class ANGLES(object):
        """ Angle class to use within object"""

        def __init__(self, nbgauss, userangfile):
            """ Init of the angles class
                nbgaus      Number of gauss angles
                userangfile Filename of the complementary list of user's angles
                """

            self.nbgauss = nbgauss
            self.userangfile = userangfile

    def __init__(self, thetas=30.0, radnb=None, raduser=None,
                 mienb=None, mieuser=None):
        """ Init the angle class
            thetas      Solar zenith angle (degrees)
            radnb       Number of Gauss angles to be used for radiance
                        computations
            raduser     Filename of the complementary list of user's angles to
                        complete the ANG.Rad.NbGauss angles (complete path).
            mienb       Number of Gauss angles to be used for mie computations
            mieuser     Filename of the complementary list of user's angles to
                        complete the ANG.mie.NbGauss angles (complete path).
            """

        self.rad = self.ANGLES(radnb, raduser)
        self.mie = self.ANGLES(mienb, mieuser)
        self.thetas = thetas


class SOS(object):
    """ SOS class definition"""

    def __init__(self, igmax=None):
        """ Init method for the SOS class
            igmax   Maximal order of interaction (scattering & surface
                    reflexion).
        """

        self.igmax = igmax


class VIEW(object):
    """ View class for the osoaa object"""

    def __init__(self, phi=0, level=5, z=-10, vza=0):
        """ This method inits the class for certain view conditions
            phi     Relative azimuth angle (degrees)
            level   Index for the output level definition
                        1 : Top of Atmosphere
                        2 : Sea Bottom
                        3 : Sea Surface 0+
                        4 : Sea Surface 0-
                        5 : User's definition of altitude or depth
                        (user data -OSOAA.View.Z)
            z       Altitude or depth (meters) for which the radiance has to be
                    given versus the viewing zenith angle (for the given
                    relative azimuth angle).
            vza     Viewing zenith angle (degrees) for which the radiance has
                    to be given versus the depth (or altitude) (for the given
                    relative azimuth angle).
            """

        self.phi = phi
        self.level = level
        self.z = z
        self.vza = vza


class OSOAA(object):
    """ This class creates the OSOAA objecto which configures and runs the
        simulation"""

    def __init__(self, wa=0.440, resroot=None):
        """ This method initiates the OSOAA class
            wa          Wavelength of radiance calculation (microns).
            resroot     Working folder for the OSOAA computations (complete
                        path).
        """

        self.wa = wa
        self.root = os.getenv("OSOAA_ROOT")
        if resroot is None:
            rnd = ''.join(random.choice(string.ascii_uppercase
                                        + string.ascii_lowercase
                                        + string.digits) for _ in range(16))
            self.resroot = self.root+"/results/"+rnd
        else:
            self.resroot = resroot

        self.sea = SEA()
        self.log = LOG()
        self.results = RESULTS()
        self.dirmie = DIRMIE(osoaaroot=self.root)
        self.phyto = PHYTO()
        self.sed = SED()
        self.sef = SED()
        self.ys = YS()
        self.det = DET()
        self.ap = AP()
        self.aer = AER()
        self.hyd = HYD()
        self.ang = ANG()
        self.sos = SOS()
        self.view = VIEW()

    def run(self):
        sc = "{}/exe/OSOAA_MAIN.exe \\".format(self.root)
        #   Definition of the working folder :
        #   ----------------------------------
        sc = sc+"\n"+"-OSOAA.ResRoot {} \\".format(self.resroot)
        #
        #   Angles calculation parameters :
        #   --------------------------------
        sc = sc+"\n"+"-ANG.Thetas {} \\".format(self.ang.thetas)
        if self.ang.rad.nbgauss is not None:
            sc = sc+"\n"+"-ANG.Rad.NbGauss {} \\".format(self.ang.rad.nbgauss)
        if self.ang.rad.userangfile is not None:
            sc = sc+"\n"+"-ANG.Rad.UserAngFile {} \\".format(
                self.ang.rad.userangfile)
        if self.results.angrad is not None:
            sc = sc+"\n"+"-ANG.Rad.ResFile {} \\".format(self.results.angrad)
        if self.ang.mie.nbgauss is not None:
            sc = sc+"\n"+"-ANG.Mie.NbGauss {} \\".format(self.ang.mie.nbgauss)
        if self.ang.mie.userangfile is not None:
            sc = sc+"\n"+"-ANG.Mie.UserAngFile {} \\".format(
                self.ang.mie.userangfile)
        if self.results.angmie is not None:
            sc = sc+"\n"+"-ANG.Mie.ResFile {} \\".format(self.results.angmie)
        if self.log.ang is not None:
            sc = sc+"\n"+"-ANG.Log {} \\".format(self.log.ang)
        #
        #   Radiance calculation parameters :
        #   --------------------------------
        if self.log.osoaa is not None:
            sc = sc+"\n"+"-OSOAA.Log {} \\".format(self.log.osoaa)
        sc = sc+"\n"+"-OSOAA.Wa  {} \\".format(self.wa)
        #
        sc = sc+"\n"+"-SEA.SurfAlb {} \\".format(self.sea.surfalb)
        sc = sc+"\n"+"-SEA.BotType {} \\".format(self.sea.bottype)
        if self.sea.bottype is 1:
            sc = sc+"\n"+"-SEA.BotAlb {} \\".format(self.sea.botalb)
        #
        sc = sc+"\n"+"-OSOAA.View.Phi {} \\".format(self.view.phi)
        sc = sc+"\n"+"-OSOAA.View.Level {} \\".format(self.view.level)
        if self.view.level is 5:
            sc = sc+"\n"+"-OSOAA.View.Z {} \\".format(self.view.z)
            sc = sc+"\n"+"-OSOAA.ResFile.Adv.Up {} \\".format(
                self.results.advup)
            sc = sc+"\n"+"-OSOAA.ResFile.Adv.Down {} \\".format(
                self.results.advdown)
            sc = sc+"\n"+"-OSOAA.View.VZA {} \\".format(self.view.vza)
            sc = sc+"\n"+"-OSOAA.ResFile.vsZ  {} \\".format(self.results.vsz)
        if self.results.vsvza is not None:
            sc = sc+"\n"+"-OSOAA.ResFile.vsVZA {} \\".format(
                self.results.vsvza)
        #
        if self.log.sos is not None:
            sc = sc+"\n"+"-SOS.Log {} \\".format(self.log.sos)
        if self.sos.igmax is not None:
            sc = sc+"\n"+"-SOS.IGmax {} \\".format(self.sos.igmax)
        if self.results.sosbin is not None:
            sc = sc+"\n"+"-SOS.ResFile.Bin {} \\".format(self.results.sosbin)
        #
        #   Profile parameters :
        #   -------------------
        if self.log.profile is not None:
            sc = sc+"\n"+"-PROFILE.Log {} \\".format(self.log.profile)
        #
        #     Atmospheric Profile parameters
        if self.results.profileatm is not None:
            sc = sc+"\n"+"-PROFILE_ATM.ResFile {} \\".format(
                self.results.profileatm)
        if self.ap.mot is not None:
            sc = sc+"\n"+"-AP.MOT {} \\".format(self.ap.mot)
        sc = sc+"\n"+"-AP.HR {} \\".format(self.ap.hr)
        if self.ap.pressure is not None:
            sc = sc+"\n"+"-AP.Pressure {} \\".format(self.ap.pressure)
        sc = sc+"\n"+"-AP.HA {} \\".format(self.ap.ha)
        #
        #     Sea Profile parameters
        if self.results.profilesea is not None:
            sc = sc+"\n"+"-PROFILE_SEA.ResFile {} \\".format(
                self.results.profilesea)
        if self.sea.depth is not None:
            sc = sc+"\n"+"-SEA.Depth {} \\".format(self.sea.depth)
        sc = sc+"\n"+"-PHYTO.Chl {} \\".format(self.phyto.chl)
        if self.phyto.chl >= 0:
            sc = sc+"\n"+"-PHYTO.ProfilType {} \\".format(
                self.phyto.profiltype)
        if self.phyto.profiltype == 2:
            sc = sc+"\n"+"-PHYTO.GP.Chlbg {} \\".format(self.phyto.gp.chlbg)
            sc = sc+"\n"+"-PHYTO.GP.Deep {} \\".format(self.phyto.gp.deep)
            sc = sc+"\n"+"-PHYTO.GP.Width {} \\".format(self.phyto.gp.width)
        if self.phyto.profiltype == 3:
            sc = sc+"\n"+"-PHYTO.Userfile {} \\".format(self.phyto.usefile)
        sc = sc+"\n"+"-SED.Csed {} \\".format(self.sed.csed)
        sc = sc+"\n"+"-YS.Abs440 {} \\".format(self.ys.abs440)
        if self.ys.abs440 > 0:
            if self.ys.swa is not None:
                sc = sc+"\n"+"-YS.Swa {} \\".format(self.ys.swa)
        sc = sc+"\n"+"-DET.Abs440 {} \\".format(self.det.abs440)
        if self.det.abs440 > 0:
            if self.det.swa is not None:
                sc = sc+"\n"+"-DET.Swa {} \\".format(self.det.swa)
        #
        #   Aerosols parameters :
        #   ---------------------
        if self.results.aer is not None:
            sc = sc+"\n"+"-AER.ResFile {} \\".format(self.results.aer)
        if self.log.aer is not None:
            sc = sc+"\n"+"-AER.Log {} \\".format(self.log.aer)
        if self.aer.aotref >= 0.0001:
            sc = sc+"\n"+"-AER.DirMie  {} \\".format(self.dirmie.aer)
        if self.log.aermie is not None:
            sc = sc+"\n"+"-AER.MieLog {} \\".format(self.log.aermie)
        sc = sc+"\n"+"-AER.Waref  {} \\".format(self.aer.waref)
        sc = sc+"\n"+"-AER.AOTref {} \\".format(self.aer.aotref)
        if self.aer.tronca is not None:
            sc = sc+"\n"+"-AER.Tronca {} \\".format(self.aer.tronca)
        if self.aer.aotref >= 0.0001:
            sc = sc+"\n"+"-AER.Model {} \\".format(self.aer.model)
        #     Aerosols parameters for mono-modal models :
        if self.aer.model is 0:
            raise ValueError('Model {} not supported '.format(self.aer.model))
            # sc = sc+"\n"+"-AER.MMD.MRwa {} \\".format()
            # sc = sc+"\n"+"-AER.MMD.MIwa {} \\".format()
            # sc = sc+"\n"+"-AER.MMD.MRwaref {} \\".format()
            # sc = sc+"\n"+"-AER.MMD.MIwaref {} \\".format()
            # sc = sc+"\n"+"-AER.MMD.SDtype {} \\".format()
            # sc = sc+"\n"+"-AER.MMD.LNDradius {} \\".format()
            # sc = sc+"\n"+"-AER.MMD.LNDvar {} \\".format()
            # sc = sc+"\n"+"-AER.MMD.JD.slope {} \\".format()
            # sc = sc+"\n"+"-AER.MMD.JD.rmin {} \\".format()
            # sc = sc+"\n"+"-AER.MMD.JD.rmax {} \\".format()
        #     Aerosols parameters for WMO models :
        elif self.aer.model is 1:
            raise ValueError('Model {} not supported '.format(self.aer.model))
            # sc = sc+"\n"+"-AER.WMO.Model {} \\".format()
            # sc = sc+"\n"+"-AER.WMO.DL {} \\".format()
            # sc = sc+"\n"+"-AER.WMO.WS {} \\".format()
            # sc = sc+"\n"+"-AER.WMO.OC {} \\".format()
            # sc = sc+"\n"+"-AER.WMO.SO {} \\".format()
        #     Aerosols parameters for Shettle&Fenn models :
        elif self.aer.model is 2:
            sc = sc+"\n"+"-AER.SF.Model {} \\".format(self.aer.sf.model)
            sc = sc+"\n"+"-AER.SF.RH {} \\".format(self.aer.sf.rh)
        #     Aerosols parameters for LND bi-modal models :
        elif self.aer.model is 3:
            raise ValueError('Model {} not supported '.format(self.aer.model))
            # sc = sc+"\n"+"-AER.BMD.VCdef {} \\".format()
            # sc = sc+"\n"+"-AER.BMD.CoarseVC {} \\".format()
            # sc = sc+"\n"+"-AER.BMD.FineVC {} \\".format()
            # sc = sc+"\n"+"-AER.BMD.RAOT {} \\".format()
            # sc = sc+"\n"+"-AER.BMD.CM.MRwa {} \\".format()
            # sc = sc+"\n"+"-AER.BMD.CM.MIwa {} \\".format()
            # sc = sc+"\n"+"-AER.BMD.CM.MRwaref {} \\".format()
            # sc = sc+"\n"+"-AER.BMD.CM.MIwaref {} \\".format()
            # sc = sc+"\n"+"-AER.BMD.CM.SDradius {} \\".format()
            # sc = sc+"\n"+"-AER.BMD.CM.SDvar {} \\".format()
            # sc = sc+"\n"+"-AER.BMD.FM.MRwa {} \\".format()
            # sc = sc+"\n"+"-AER.BMD.FM.MIwa {} \\".format()
            # sc = sc+"\n"+"-AER.BMD.FM.MRwaref {} \\".format()
            # sc = sc+"\n"+"-AER.BMD.FM.MIwaref {} \\".format()
            # sc = sc+"\n"+"-AER.BMD.FM.SDradius {} \\".format()
            # sc = sc+"\n"+"-AER.BMD.FM.SDvar {} \\".format()
        #    Aerosols parameters for external data (phase functions, scattering
        #    and extinction coefficients) :
        elif self.aer.model is 4:
            raise ValueError('Model {} not supported '.format(self.aer.model))
            # sc = sc+"\n"+"-AER.ExtData {} \\".format()
        #
        #   Hydrosols parameters :
        #   ---------------------
        if self.results.phyto is not None:
            sc = sc+"\n"+"-PHYTO.ResFile {} \\".format(self.results.phyto)
        if self.results.mlp is not None:
            sc = sc+"\n"+"-MLP.ResFile {} \\".format(self.results.mlp)
        if self.log.hyd is not None:
            sc = sc+"\n"+"-HYD.Log {} \\".format(self.log.hyd)
        if self.phyto.chl > 0 or self.sed.csed > 0:
            sc = sc+"\n"+"-HYD.DirMie {} \\".format(self.dirmie.hid)
        if self.log.hydmie is not None:
            sc = sc+"\n"+"-HYD.MieLog {} \\".format(self.log.hydmie)
        if self.phyto.chl > 0 or self.sed.csed > 0:
            sc = sc+"\n"+"-HYD.Model {} \\".format(self.hyd.model)
        #     Phytoplankton model :
        if self.hyd.model is 1:
            #     Junge main mode :
            if self.phyto.jd is not None:
                sc = sc+"\n"+"-PHYTO.JD.slope {} \\".format(
                    self.phyto.jd.slope)
                sc = sc+"\n"+"-PHYTO.JD.rmin {} \\".format(self.phyto.jd.rmin)
                sc = sc+"\n"+"-PHYTO.JD.rmax {} \\".format(self.phyto.jd.rmax)
                sc = sc+"\n"+"-PHYTO.JD.MRwa {} \\".format(self.phyto.jd.mrwa)
                sc = sc+"\n"+"-PHYTO.JD.MIwa {} \\".format(self.phyto.jd.miwa)
                sc = sc+"\n"+"-PHYTO.JD.rate {} \\".format(self.phyto.jd.rate)
            #     Secondary LND mode :
            if self.phyto.sm is not None:
                sc = sc+"\n"+"-PHYTO.LND.SM.SDradius {} \\".format(
                    self.phyto.sm.sdradius)
                sc = sc+"\n"+"-PHYTO.LND.SM.SDvar {} \\".format(
                    self.phyto.sm.sdvar)
                sc = sc+"\n"+"-PHYTO.LND.SM.MRwa {} \\".format(
                    self.phyto.sm.mrwa)
                sc = sc+"\n"+"-PHYTO.LND.SM.MIwa {} \\".format(
                    self.phyto.sm.miwa)
                sc = sc+"\n"+"-PHYTO.LND.SM.rate {} \\".format(
                    self.phyto.sm.rate)
            #     Tertiary LND mode :"
            if self.phyto.tm is not None:
                sc = sc+"\n"+"-PHYTO.LND.TM.SDradius {} \\".format(
                    self.phyto.tm.sdradius)
                sc = sc+"\n"+"-PHYTO.LND.TM.SDvar {} \\".format(
                    self.phyto.tm.sdvar)
                sc = sc+"\n"+"-PHYTO.LND.TM.MRwa {} \\".format(
                    self.phyto.tm.mrwa)
                sc = sc+"\n"+"-PHYTO.LND.TM.MIwa {} \\".format(
                    self.phyto.tm.miwa)
                sc = sc+"\n"+"-PHYTO.LND.TM.rate {} \\".format(
                    self.phyto.tm.rate)
        if self.sed.csed > 0.0:
            #     Mineral-like particles model :
            #     Junge main mode :
            if self.sed.jd is not None:
                sc = sc+"\n"+"-SED.JD.slope {} \\".format(self.sed.jd.slope)
                sc = sc+"\n"+"-SED.JD.rmin {} \\".format(self.sed.jd.rmin)
                sc = sc+"\n"+"-SED.JD.rmax {} \\".format(self.sed.jd.rmax)
                sc = sc+"\n"+"-SED.JD.MRwa {} \\".format(self.sed.jd.mrwa)
                sc = sc+"\n"+"-SED.JD.MIwa {} \\".format(self.sed.jd.miwa)
                sc = sc+"\n"+"-SED.JD.rate {} \\".format(self.sed.jd.rate)
            #     Secondary LND mode :
            if self.sed.sm is not None:
                sc = sc+"\n"+"-SED.LND.SM.SDradius {} \\".format(
                    self.sed.sm.sdradius)
                sc = sc+"\n"+"-SED.LND.SM.SDvar {} \\".format(
                    self.sed.sm.sdvar)
                sc = sc+"\n"+"-SED.LND.SM.MRwa {} \\".format(self.sed.sm.mrwa)
                sc = sc+"\n"+"-SED.LND.SM.MIwa {} \\".format(self.sed.sm.miwa)
                sc = sc+"\n"+"-SED.LND.SM.rate {} \\".format(self.sed.sm.rate)
            #     Tertiary LND mode :
            if self.sed.tm is not None:
                sc = sc+"\n"+"-SED.LND.TM.SDradius {} \\".format(
                    self.sed.tm.sdradius)
                sc = sc+"\n"+"-SED.LND.TM.SDvar {} \\".format(
                    self.sed.tm.sdvar)
                sc = sc+"\n"+"-SED.LND.TM.MRwa {} \\".format(self.sed.tm.mrwa)
                sc = sc+"\n"+"-SED.LND.TM.MIwa {} \\".format(self.sed.tm.miwa)
                sc = sc+"\n"+"-SED.LND.TM.rate {} \\".format(self.sed.tm.rate)
        #     Hydrosols parameters for external data (phase functions,
        #     scattering and extinction coefficients) :
        if self.hyd.model is 2:
            sc = sc+"\n"+"-HYD.ExtData {} \\".format(self.hyd.extdata)
        #
        #   Sea / atmosphere interface parameters :
        #   --------------------------------------
        if self.log.sea is not None:
            sc = sc+"\n"+"-SEA.Log {} \\".format(self.log.sea)
        if self.sea.wind > 0:
            sc = sc+"\n"+"-SEA.Dir {} \\".format(self.dirmie.sea)
            sc = sc+"\n"+"-SEA.Ind {} \\".format(self.sea.ind)
        sc = sc+"\n"+"-SEA.Wind {} ".format(self.sea.wind)
        # Check if directory exists

        if not os.path.exists(self.resroot):
            os.makedirs(self.resroot)
        if not os.path.exists(self.dirmie.aer):
            os.makedirs(self.dirmie.aer)
        if not os.path.exists(self.dirmie.hid):
            os.makedirs(self.dirmie.hid)
        if not os.path.exists(self.dirmie.sea):
            os.makedirs(self.dirmie.sea)

        # We generate the script
        with open(self.resroot+"/script.kzh", 'w') as file:
            file.write(sc)

        # Run script with ksh
        os.system("ksh "+self.resroot+"/script.kzh")

        # read OUTPUTS
        self.outputs = OUTPUTS(self.resroot, self.results)
