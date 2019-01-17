import numpy as np
import os



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
        pytho       Filename of the result OSOAA_HYDROSOLS computations or
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
    pytho = None
    mlp = None
    angrad = None
    angmie = None
    sosbin = None
    vsvza = None
    advup = "resfile_advup.txt"
    advdown = "resfile_advdown.txt"
    vsz = "resfile_vsz.txt"


class DIRMIE(object):
    """ Directory for hydrosol MIE files storage (complete path)
        aer         Storage directory for MIE files producted by OSOAA_AEROSOLS
                    computations (complete path).
        hid         Storage directory for MIE files producted by
                    HYDROSOLS_AEROSOLS computations (complete path).
        SEA         Directory for SURFACE files storage (complete path).
    """
    hid = "/DATABASE/MIE_HYD"
    aer = "/DATABASE/MIE_AER"
    sea = "/DATABASE/SURF_MATR"


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
        self.jd = JD()
        self.sm = None
        self.tm = None

    def SetPrimaryMode(self, self, mrwa=1.05, miwa=-0.0, slope=4.0,
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

        self.jd = JD(mrwa, miwa, slope, rmin, rmax, rate)

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

        self.sm = LND(mrwa, miwa, sdradius, sdvar, rate)

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

        self.jd = JD(mrwa, miwa, slope, rmin, rmax, rate)

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

        self.sm = LND(mrwa, miwa, sdradius, sdvar, rate)

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

        self.tm = LND(mrwa, miwa, sdradius, sdvar, rate)


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

    def __init__(self, waref=0.0, aotref=0.1, tronca=None, model=2):
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
            """

        self.waref = waref
        self.aotref = aotref
        self.tronca = tronca
        self.model = model


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
        if resroot is None:
            self.resroot = os.getenv("OSOAA_ROOT")
        else:
            self.resroot = resroot

        self.sea = SEA()
        self.log = LOG()
        self.results = RESULTS()
        self.dirmie = DIRMIE()
        self.phyto = PHYTO()
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
