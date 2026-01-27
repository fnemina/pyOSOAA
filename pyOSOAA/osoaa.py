# coding=utf-8

import hashlib
import os
import platform
import shlex
import shutil
import subprocess
from tempfile import TemporaryDirectory

from .outputs import OUTPUTS


class SEA(object):
    """This is the SEA class which defines the interfaces at the bottom of the
    ocean and the interface with the air."""

    def __init__(
        self, surfalb=0.0, bottype=1, botalb=0.30, wind=7, ind=1.34, depth=15.0
    ):
        """Init function for the SEA class.
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
    """Log files for the simulation.
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
    hydmie = "log_hydmie.txt"
    sea = "log_sea.txt"
    sos = "log_sos.txt"


class RESULTS(object):
    """Result files for the simulation.
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
    advphi      Filename of the ascii file resulting from SOS computations
                (defined without directory tree ==> this file will be
                located in the sub-directory Advanced_outputs of the
                working folder) : ==> Advanced output upward radiance field
                versus the depth (or altitude) AND versus the viewing
                zenith angle AND versus the relative azimuth angle).
                This parameter only works with a custom OSOAA version from
                https://github.com/fnemina/RadiativeTransferCode-OSOAA
                in the branch advphi.
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
    advup = None
    advdown = None
    advphi = None
    vsz = "resfile_vsz.txt"


class DIRMIE(object):
    """Directory for hydrosol MIE files storage"""

    def __init__(
        self,
        osoaaroot,
        hid=os.path.join("DATABASE", "MIE_HYD"),
        aer=os.path.join("DATABASE", "MIE_AER"),
        sea=os.path.join("DATABASE", "SURF_MATR"),
    ):
        """Directory for hydrosol MIE files storage (complete path)
        aer         Storage directory for MIE files producted by
                    OSOAA_AEROSOLS computations (complete path).
        hid         Storage directory for MIE files producted by
                    HYDROSOLS_AEROSOLS computations (complete path).
        SEA         Directory for SURFACE files storage (complete path).
        """
        self.hyd = os.path.join(osoaaroot, hid)
        self.aer = os.path.join(osoaaroot, aer)
        self.sea = os.path.join(osoaaroot, sea)


class GP(object):
    """Gaussian profile class definition."""

    def __init__(self, chlbg, deep, width):
        """Init function for the gaussian profiles
        chlbg       Constant biomass background (mg/m3)
        deep        Maximum deep of the gaussian chlorophyll profile (m)
        width       Peak width of the gaussian chlorophyll profile (m)
        """

        self.chlbg = chlbg
        self.deep = deep
        self.width = width


class PHYTO(object):
    """Phytoplanckton class"""

    # Profile types
    Homogeneous = 1
    Gaussian = 2
    UserDefined = 3

    class JD(object):
        """This is a Junge distribution for the different models."""

        def __init__(
            self, mrwa=1.05, miwa=-0.0, slope=4.0, rmin=0.01, rmax=200, rate=1.0
        ):
            """Init function for the Junges distribution.
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
        """This is a lognormal distribution for the different models."""

        def __init__(self, mrwa, miwa, sdradius, sdvar, rate):
            """Init function for the lognormal distribution.
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
        """Init function for the Phytoplanckton profiles
        chl         chlorophyll concentration at sea surface (mg/m3)
        """

        self.chl = chl
        self.profiltype = self.Homogeneous
        self.jd = self.JD()
        self.sm = None
        self.tm = None

    def SetPrimaryMode(
        self, mrwa=1.05, miwa=-0.0, slope=4.0, rmin=0.01, rmax=200, rate=1.0
    ):
        """Sets the primary mode using Junge's law
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
        """Sets the secondary mode using lognormal distribution
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
        """Sets the tertiary mode using lognormal distribution
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

    def SetProfilType(
        self, profiltype, chlbg=None, deep=None, width=None, userfile=None
    ):
        """This method sets the profile type for the Phytoplanckton
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

        if profiltype not in [self.Homogeneous, self.Gaussian, self.UserDefined]:
            raise Exception("Invalid profile type.")

        self.profiltype = profiltype

        if profiltype == self.Gaussian:
            # We confifigure the gaussian profile
            self.gp = GP(chlbg, deep, width)

        elif profiltype == self.UserDefined:
            # We configure the user defined profile

            try:
                tmp = open(userfile, "r")
                tmp.close()
                self.usefile = userfile
            except FileNotFoundError:
                print("File {} not found".format(userfile))


class SED(object):
    """Sediment class"""

    class JD(object):
        """This is a Junge distribution for the different models."""

        def __init__(
            self, mrwa=1.15, miwa=0.0, slope=4.0, rmin=0.01, rmax=200, rate=1.0
        ):
            """Init function for the Junges distribution.
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
        """This is a lognormal distribution for the different models."""

        def __init__(self, mrwa, miwa, sdradius, sdvar, rate):
            """Init function for the lognormal distribution.
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
        """Init function for the sediment profiles
        csed    Concentration of mineral-like particles at sea surface
                (mg/liter)
        """

        self.csed = csed
        self.jd = self.JD()
        self.sm = None
        self.tm = None

    def SetPrimaryMode(self, mrwa=1.2, miwa=0, slope=4, rmin=0.01, rmax=200, rate=1):
        """Sets the primary mode using Junge's law
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
        """Sets the secondary mode using lognormal distribution
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
        """Sets the secondary mode using lognormal distribution
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
    """Absorption class to be used with yellow substance"""

    def __init__(self, abs440=0.0, swa=None):
        """Init function for the absorption class.ABS
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
    """Absorption class to be used with detritus"""

    def __init__(self, abs440=0.0, swa=None):
        """Init function for the absorption class.ABS
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
    """Atmospheric profile parameters object."""

    def __init__(self, mot=None, pressure=1013.00, hr=8.0, ha=2.0):
        """Init function for the atmospheric profile
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
        """Define the molecular optical thickness with pressure
        pressure    atmospheric pressure at sea level (mbar)"""
        self.mot = None
        self.pressure = pressure

    def SetMot(self, mot=0.1, hr=8):
        """Define the molecular optical thickness
        mot         molecular optical thicknes
        hr          height scale for molecules (km)"""
        self.mot = mot
        self.hr = hr
        self.pressure = None


class AEROSOLMODELS(object):
    """Aerosol models for the AER class."""

    class SF(object):
        """Shettle and Fenn atmosphere model class."""

        # Shettle and Fenn models
        Tropospheric = 1
        Urban = 2
        Maritime = 3
        Coastal = 4

        def __init__(self, sfmodel, rh):
            """Init method for the Shettle-Fenn model.
            model       Type of Shettle & Fenn model.
                            1 : Tropospheric S&F model.
                            2 : Urban S&F model.
                            3 : Maritime S&F model.
                            4 : Coastal S&F model.
            rh          Relative humidity (%) for Shettle & Fenn model.
            """
            self.model = sfmodel
            self.rh = rh

    class MM(object):
        """Mono-modal size distribution"""

        def __init__(self, sdtype):
            """Init method for the mono-modal size distribution
            sdtype      Type of mono-modal size distribution
                            1 : Log Normal size distribution
                            2 : Junge's law
            lnd         Log normal size distribution
            jd          Junge's law size distribution
            mrwa        Real part of the aerosol refractive index for the
                        wavelength of radiation calculation
            miwa        Imaginary part of the aerosol refractive index for
                        the  wavelength of radiation calculation
            sdradius    Modal radius (um) of the Log-Noprmal size
                        distribution
            sdvar       Standard deviation of the Log-Normal size
                        distribution
            slope       Slope of the Junge's law.
                        Warning: 3 is a singular value.
            rmin        Minimal radius of Junge's law (um)
            rmax        Maximal radius of Junge's law (um)
            mrwaref     Real part of the aerosol refractive index for the
                        reference wavelength of aerosol properties
                        calculation.
            miwaref     Imaginary part of the aerosol refractive index for
                        the reference wavelength of aerosol properties
                        calculation.
            """

            self.sdtype = sdtype
            self.mrwa = None
            self.miwa = None
            self.mrwaref = None
            self.miwaref = None
            if sdtype == 1:
                self.sdradius = None
                self.sdvar = None
            elif sdtype == 2:
                self.slope = None
                self.rmin = None
                self.rmax = None

    class WMO(object):
        """WMO aerosol models."""

        def __init__(self, wmotype, dl=None, ws=None, oc=None, so=None):
            """Init method for the WMO aerosol model
            wmotype     Type of WMO model
                            1 : Continental WMO model
                            2 : Maritime WMO model
                            3 : Urban WMO model
                            4 : WMO model by usef definition
            dl          Volume concentration (between 0 and 1) for
                        "Dust like" components
            ws          Volume concentration (between 0 and 1) for
                        "Water soluble" components
            oc          Volume concentration (between 0 and 1) for
                        "Oceanic" components
            so          Volume concentration (between 0 and 1) for
                        "Soot" components
            """

            self.model = wmotype

            if wmotype == 4:
                self.dl = dl
                self.ws = ws
                self.oc = oc
                self.so = so

    class LNB(object):
        """Log-normal bi-modal aerosol model"""

        def __init__(self, vcdef):
            """Log-normal bi-modal aerosol model init functions
            vcdef       Choide of the mixture description type
                            1 : Use of predefined volume concentrations.
                            2 : Use of the ratio of aerosol optical
                                thickness (coarse mode ATO / total AOT)
            coarsevc    User volume concentration of the LND coarse mode
            finevc      User volume concentration of the LND fine mode
            raot        User value of the ration AOT_coarse/AOT_total for
                        the aerosol reference wavelength

            cmrwa       Real part of the aerosol refractive index for the
                        wavelength of radiation calculation for the coarse
                        mode
            cmiwa       Imaginary part of the aerosol refractive index for
                        the  wavelength of radiation calculationfor the
                        coarse mode
            csdradius   Modal radius (um) of the Log-Normal size
                        distribution for the coarse mode
            csdvar      Standard deviation of the Log-Normal size
                        distribution for the coarse mode
            cmrwaref    Real part of the aerosol refractive index for the
                        reference wavelength of aerosol properties
                        calculation for the coarse mode
            cmiwaref    Imaginary part of the aerosol refractive index for
                        the reference wavelength of aerosol properties
                        calculation for the coarse mode

            fmrwa       Real part of the aerosol refractive index for the
                        wavelength of radiation calculation for the fine
                        mode
            fmiwa       Imaginary part of the aerosol refractive index for
                        the  wavelength of radiation calculationfor the
                        fine mode
            fsdradius   Modal radius (um) of the Log-Noprmal size
                        distribution for the fine mode
            fsdvar      Standard deviation of the Log-Normal size
                        distribution for the fine mode
            fmrwaref    Real part of the aerosol refractive index for the
                        reference wavelength of aerosol properties
                        calculation for the fine mode
            fmiwaref    Imaginary part of the aerosol refractive index for
                        the reference wavelength of aerosol properties
                        calculation for the fine mode
            """

            self.vcdef = vcdef
            if vcdef == 1:
                self.coarsevc = None
                self.finevc = None
            elif vcdef == 2:
                self.raot = None

            self.cmrwa = None
            self.cmiwa = None
            self.cmrwaref = None
            self.cmiwaref = None
            self.csdradius = None
            self.csdvar = None

            self.fmrwa = None
            self.fmiwa = None
            self.fmrwaref = None
            self.fmiwaref = None
            self.fsdradius = None
            self.fsdvar = None


class AER(object):
    """This class contains everything related to the aerosol components
    of the atmosphere."""

    def __init__(self, waref=0.550, aotref=0.1, tronca=None, model=2):
        """Init method for the aerosol componentes class
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
        self.sf = AEROSOLMODELS.SF(sfmodel=3, rh=98)

    def SetModel(
        self,
        model=2,
        sdtype=1,
        sfmodel=3,
        rh=98,
        wmotype=1,
        dl=None,
        ws=None,
        oc=None,
        so=None,
        vcdef=2,
        extdata="",
    ):
        """This methods sets the model for the AER class.

        model       Type of aerosol model
                        0 : Mono-modal
                        1 : WMO multi-modal
                        2 : Shettle & Fenn bi-modal
                        3 : Log-Normal bi-modal
                        4 : Phase function from an external source

        Mono-modal distribution parameters
        ----------------------------------
        mm          Mono-modal model atribute
        sdtype      Type of mono-modal size distribution
                        1 : Log Normal size distribution
                        2 : Junge's law

        WMO model parameters
        -------------------
        wmo         WMO model atribute
        wmotype     Type of WMO model
                        1 : Continental WMO model
                        2 : Maritime WMO model
                        3 : Urban WMO model
                        4 : WMO model by usef definition
        dl          Volume concentration (between 0 and 1) for
                    "Dust like" components
        ws          Volume concentration (between 0 and 1) for
                    "Water soluble" components
        oc          Volume concentration (between 0 and 1) for
                    "Oceanic" components
        so          Volume concentration (between 0 and 1) for
                    "Soot" components

        Shettle and Fenn model parameters
        ---------------------------------
        sf          Shettle and Fenn model atribute
        sfmodel       Type of Shettle & Fenn model.
                        1 : Tropospheric S&F model.
                        2 : Urban S&F model.
                        3 : Maritime S&F model.
                        4 : Coastal S&F model.
        rh          Relative humidity (%) for Shettle & Fenn model.

        Log-Normal bi-modal model parameters
        ------------------------------------
        lnd         Log-Normal bi-modal model atribute
        vcdef       Choide of the mixture description type
                        1 : Use of predefined volume concentrations.
                        2 : Use of the ratio of aerosol optical
                            thickness (coarse mode ATO / total AOT)

        External phase function
        -----------------------
        extdata     Filename (complete path) of user's external phase
                    function data and radiative parameters (extinction and
                    scattering coefficients).
                    The reference aerosol wavelength and the radiance
                    simulation wavelength must be equal
        """
        self.model = model
        if model == 0:
            self.mm = AEROSOLMODELS.MM(sdtype)
            self.wmo = None
            self.sf = None
            self.lnd = None
            self.external = None
        elif model == 1:
            self.mm = None
            self.wmo = AEROSOLMODELS.WMO(wmotype, dl, ws, oc, so)
            self.sf = None
            self.lnd = None
            self.external = None
        elif model == 2:
            self.mm = None
            self.wmo = None
            self.sf = AEROSOLMODELS.SF(sfmodel, rh)
            self.lnd = None
            self.external = None
        elif model == 3:
            self.mm = None
            self.wmo = None
            self.sf = None
            self.lnb = AEROSOLMODELS.LNB(vcdef)
            self.external = None
        elif model == 4:
            self.mm = None
            self.sf = None
            self.wmo = None
            self.lnd = None
            self.extdata = extdata


class HYD(object):
    """This class contains everything related to the hydrosol components
    of the sea."""

    def __init__(self, model=1, extdata=None):
        """Init method for the aerosol componentes class
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
    """Angle definitions class."""

    class ANGLES(object):
        """Angle class to use within object"""

        def __init__(self, nbgauss, userangfile):
            """Init of the angles class
            nbgaus      Number of gauss angles
            userangfile Filename of the complementary list of user's angles
            """

            self.nbgauss = nbgauss
            self.userangfile = userangfile

    def __init__(self, thetas=30.0, radnb=None, raduser=None, mienb=None, mieuser=None):
        """Init the angle class
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
    """SOS class definition"""

    def __init__(self, igmax=None):
        """Init method for the SOS class
        igmax   Maximal order of interaction (scattering & surface
                reflexion).
        """

        self.igmax = igmax


class VIEW(object):
    """View class for the osoaa object"""

    def __init__(self, phi=0, level=5, z=-10, vza=0):
        """This method inits the class for certain view conditions
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
    """This class creates the OSOAA objecto which configures and runs the
    simulation"""

    def __init__(self, wa=0.440, resroot=None, logfile=None, cleanup=False):
        """This method initiates the OSOAA class
        wa          Wavelength of radiance calculation (microns).
        resroot     Working folder for the OSOAA computations (complete
                    path).
        logfile     Logfile name to output the results.
        cleanup     True to delete the directories with the results.
                    False by default.
        """
        self.salt = ''
        # Determine the type of operating system
        self.os = platform.system()

        self.wa = wa
        self.root = os.getenv("OSOAA_ROOT")
        self.resroot = resroot

        self.cleanup = cleanup

        self.sea = SEA()
        self.log = LOG()
        self.results = RESULTS()
        self.dirmie = DIRMIE(osoaaroot=self.root)
        self.phyto = PHYTO()
        self.sed = SED()
        self.ys = YS()
        self.det = DET()
        self.ap = AP()
        self.aer = AER()
        self.hyd = HYD()
        self.ang = ANG()
        self.sos = SOS()
        self.view = VIEW()
        self.logfile = logfile

    def _bake_arguments(self, resroot, fatm_null=False):
        if fatm_null:
            exe = "OSOAA_MAIN_FATM_NULL.exe"
        else:
            exe = "OSOAA_MAIN.exe"

        sc = [os.path.join(self.root, "exe", exe)]
        #
        #   Angles calculation parameters :
        #   --------------------------------
        sc += ["-ANG.Thetas", str(self.ang.thetas)]
        if self.ang.rad.nbgauss is not None:
            sc += ["-ANG.Rad.NbGauss", str(self.ang.rad.nbgauss)]
        if self.ang.rad.userangfile is not None:
            sc += ["-ANG.Rad.UserAngFile", str(self.ang.rad.userangfile)]
        if self.results.angrad is not None:
            sc += ["-ANG.Rad.ResFile", str(self.results.angrad)]
        if self.ang.mie.nbgauss is not None:
            sc += ["-ANG.Mie.NbGauss", str(self.ang.mie.nbgauss)]
        if self.ang.mie.userangfile is not None:
            sc += ["-ANG.Mie.UserAngFile", str(self.ang.mie.userangfile)]
        if self.results.angmie is not None:
            sc += ["-ANG.Mie.ResFile", str(self.results.angmie)]
        if self.log.ang is not None:
            sc += ["-ANG.Log", str(self.log.ang)]
        #
        #   Radiance calculation parameters :
        #   --------------------------------
        if self.log.osoaa is not None:
            sc += ["-OSOAA.Log", str(self.log.osoaa)]
        sc += ["-OSOAA.Wa", str(self.wa)]
        #
        sc += ["-SEA.SurfAlb", str(self.sea.surfalb)]
        sc += ["-SEA.BotType", str(self.sea.bottype)]
        if self.sea.bottype == 1:
            sc += ["-SEA.BotAlb", str(self.sea.botalb)]
        #
        sc += ["-OSOAA.View.Phi", str(self.view.phi)]
        sc += ["-OSOAA.View.Level", str(self.view.level)]
        if self.results.advup is not None:
            sc += ["-OSOAA.ResFile.Adv.Up", str(self.results.advup)]
        if self.results.advdown is not None:
            sc += ["-OSOAA.ResFile.Adv.Down", str(self.results.advdown)]
        if self.results.advphi is not None:
            sc += ["-OSOAA.ResFile.Adv.Phi", str(self.results.advphi)]
        if self.view.level == 5:
            sc += ["-OSOAA.View.Z", str(self.view.z)]
            sc += ["-OSOAA.View.VZA", str(self.view.vza)]
            sc += ["-OSOAA.ResFile.vsZ", str(self.results.vsz)]
        if self.results.vsvza is not None:
            sc += ["-OSOAA.ResFile.vsVZA", str(self.results.vsvza)]
        #
        if self.log.sos is not None:
            sc += ["-SOS.Log", str(self.log.sos)]
        if self.sos.igmax is not None:
            sc += ["-SOS.IGmax", str(self.sos.igmax)]
        if self.results.sosbin is not None:
            sc += ["-SOS.ResFile.Bin", str(self.results.sosbin)]
        #
        #   Profile parameters :
        #   -------------------
        if self.log.profile is not None:
            sc += ["-PROFILE.Log", str(self.log.profile)]
        #
        #     Atmospheric Profile parameters
        if self.results.profileatm is not None:
            sc += ["-PROFILE_ATM.ResFile", str(self.results.profileatm)]
        if self.ap.mot is not None:
            sc += ["-AP.MOT", str(self.ap.mot)]
        sc += ["-AP.HR", str(self.ap.hr)]
        if self.ap.pressure is not None:
            sc += ["-AP.Pressure", str(self.ap.pressure)]
        sc += ["-AP.HA", str(self.ap.ha)]
        #
        #     Sea Profile parameters
        if self.results.profilesea is not None:
            sc += ["-PROFILE_SEA.ResFile", str(self.results.profilesea)]
        if self.sea.depth is not None:
            sc += ["-SEA.Depth", str(self.sea.depth)]
        sc += ["-PHYTO.Chl", str(self.phyto.chl)]
        if self.phyto.chl >= 0:
            sc += ["-PHYTO.ProfilType", str(self.phyto.profiltype)]
        if self.phyto.profiltype == 2:
            sc += ["-PHYTO.GP.Chlbg", str(self.phyto.gp.chlbg)]
            sc += ["-PHYTO.GP.Deep", str(self.phyto.gp.deep)]
            sc += ["-PHYTO.GP.Width", str(self.phyto.gp.width)]
        if self.phyto.profiltype == 3:
            sc += ["-PHYTO.Userfile", str(self.phyto.usefile)]
        sc += ["-SED.Csed", str(self.sed.csed)]
        sc += ["-YS.Abs440", str(self.ys.abs440)]
        if self.ys.abs440 > 0:
            if self.ys.swa is not None:
                sc += ["-YS.Swa", str(self.ys.swa)]
        sc += ["-DET.Abs440", str(self.det.abs440)]
        if self.det.abs440 > 0:
            if self.det.swa is not None:
                sc += ["-DET.Swa", str(self.det.swa)]
        #
        #   Aerosols parameters :
        #   ---------------------
        if self.results.aer is not None:
            sc += ["-AER.ResFile", str(self.results.aer)]
        if self.log.aer is not None:
            sc += ["-AER.Log", str(self.log.aer)]
        if self.aer.aotref >= 0.0:
            # FIXME
            sc += ["-AER.DirMie", str(self.dirmie.aer)]
        if self.log.aermie is not None:
            sc += ["-AER.MieLog", str(self.log.aermie)]
        sc += ["-AER.Waref", str(self.aer.waref)]
        sc += ["-AER.AOTref", str(self.aer.aotref)]
        if self.aer.tronca is not None:
            sc += ["-AER.Tronca", str(self.aer.tronca)]
        if self.aer.aotref > 0.0:
            sc += ["-AER.Model", str(self.aer.model)]
        #     Aerosols parameters for mono-modal models :
        if self.aer.model == 0:
            sc += ["-AER.MMD.MRwa", str(self.aer.mm.mrwa)]
            sc += ["-AER.MMD.MIwa", str(self.aer.mm.miwa)]
            if self.wa is not self.aer.waref:
                sc += ["-AER.MMD.MRwaref", str(self.aer.mm.mrwaref)]
                sc += ["-AER.MMD.MIwaref", str(self.aer.mm.miwaref)]
            sc += ["-AER.MMD.SDtype", str(self.aer.mm.sdtype)]
            if self.aer.mm.sdtype == 1:
                sc += ["-AER.MMD.LNDradius", str(self.aer.mm.sdradius)]
                sc += ["-AER.MMD.LNDvar", str(self.aer.mm.sdvar)]
            elif self.aer.mm.sdtype == 2:
                sc += ["-AER.MMD.JD.slope", str(self.aer.mm.slope)]
                sc += ["-AER.MMD.JD.rmin", str(self.aer.mm.rmin)]
                sc += ["-AER.MMD.JD.rmax", str(self.aer.mm.rmax)]
        #     Aerosols parameters for WMO models :
        elif self.aer.model == 1:
            sc += ["-AER.WMO.Model", str(self.aer.wmo.model)]
            if self.aer.wmo.model == 4:
                sc += ["-AER.WMO.DL", str(self.aer.wmo.dl)]
                sc += ["-AER.WMO.WS", str(self.aer.wmo.ws)]
                sc += ["-AER.WMO.OC", str(self.aer.wmo.oc)]
                sc += ["-AER.WMO.SO", str(self.aer.wmo.so)]
        #     Aerosols parameters for Shettle&Fenn models :
        elif self.aer.model == 2:
            sc += ["-AER.SF.Model", str(self.aer.sf.model)]
            sc += ["-AER.SF.RH", str(self.aer.sf.rh)]
        #     Aerosols parameters for LND bi-modal models :
        elif self.aer.model == 3:
            sc += ["-AER.BMD.VCdef", str(self.aer.lnb.vcdef)]
            if self.aer.lnb.vcdef == 1:
                sc += ["-AER.BMD.CoarseVC", str(self.aer.lnb.coarsevc)]
                sc += ["-AER.BMD.FineVC", str(self.aer.lnb.finevc)]
            elif self.aer.lnb.vcdef == 2:
                sc += ["-AER.BMD.RAOT", str(self.aer.lnb.raot)]
            sc += ["-AER.BMD.CM.MRwa", str(self.aer.lnb.cmrwa)]
            sc += ["-AER.BMD.CM.MIwa", str(self.aer.lnb.cmiwa)]
            sc += ["-AER.BMD.CM.MRwaref", str(self.aer.lnb.cmrwaref)]
            sc += ["-AER.BMD.CM.MIwaref", str(self.aer.lnb.cmiwaref)]
            sc += ["-AER.BMD.CM.SDradius", str(self.aer.lnb.csdradius)]
            sc += ["-AER.BMD.CM.SDvar", str(self.aer.lnb.csdvar)]
            sc += ["-AER.BMD.FM.MRwa", str(self.aer.lnb.fmrwa)]
            sc += ["-AER.BMD.FM.MIwa", str(self.aer.lnb.fmiwa)]
            sc += ["-AER.BMD.FM.MRwaref", str(self.aer.lnb.fmrwaref)]
            sc += ["-AER.BMD.FM.MIwaref", str(self.aer.lnb.fmiwaref)]
            sc += ["-AER.BMD.FM.SDradius", str(self.aer.lnb.fsdradius)]
            sc += ["-AER.BMD.FM.SDvar", str(self.aer.lnb.fsdvar)]
        #    Aerosols parameters for external data (phase functions, scattering
        #    and extinction coefficients) :
        elif self.aer.model == 4:
            sc += ["-AER.ExtData", str(self.aer.extdata)]
        #
        #   Hydrosols parameters :
        #   ---------------------
        if self.results.phyto is not None:
            sc += ["-PHYTO.ResFile", str(self.results.phyto)]
        if self.results.mlp is not None:
            sc += ["-MLP.ResFile", str(self.results.mlp)]
        if self.log.hyd is not None:
            sc += ["-HYD.Log", str(self.log.hyd)]
        if self.phyto.chl > 0 or self.sed.csed > 0:
            sc += ["-HYD.DirMie", str(self.dirmie.hyd)]
        if self.log.hydmie is not None:
            sc += ["-HYD.MieLog", str(self.log.hydmie)]
        if self.phyto.chl > 0 or self.sed.csed > 0:
            sc += ["-HYD.Model", str(self.hyd.model)]
        #     Phytoplankton model :
        if self.hyd.model == 1:
            #     Junge main mode :
            if self.phyto.jd is not None:
                sc += ["-PHYTO.JD.slope", str(self.phyto.jd.slope)]
                sc += ["-PHYTO.JD.rmin", str(self.phyto.jd.rmin)]
                sc += ["-PHYTO.JD.rmax", str(self.phyto.jd.rmax)]
                sc += ["-PHYTO.JD.MRwa", str(self.phyto.jd.mrwa)]
                sc += ["-PHYTO.JD.MIwa", str(self.phyto.jd.miwa)]
                sc += ["-PHYTO.JD.rate", str(self.phyto.jd.rate)]
            #     Secondary LND mode :
            if self.phyto.sm is not None:
                sc += ["-PHYTO.LND.SM.SDradius", str(self.phyto.sm.sdradius)]
                sc += ["-PHYTO.LND.SM.SDvar", str(self.phyto.sm.sdvar)]
                sc += ["-PHYTO.LND.SM.MRwa", str(self.phyto.sm.mrwa)]
                sc += ["-PHYTO.LND.SM.MIwa", str(self.phyto.sm.miwa)]
                sc += ["-PHYTO.LND.SM.rate", str(self.phyto.sm.rate)]
            #     Tertiary LND mode :"
            if self.phyto.tm is not None:
                sc += ["-PHYTO.LND.TM.SDradius", str(self.phyto.tm.sdradius)]
                sc += ["-PHYTO.LND.TM.SDvar", str(self.phyto.tm.sdvar)]
                sc += ["-PHYTO.LND.TM.MRwa", str(self.phyto.tm.mrwa)]
                sc += ["-PHYTO.LND.TM.MIwa", str(self.phyto.tm.miwa)]
                sc += ["-PHYTO.LND.TM.rate", str(self.phyto.tm.rate)]
        if self.sed.csed > 0.0:
            #     Mineral-like particles model :
            #     Junge main mode :
            if self.sed.jd is not None:
                sc += ["-SED.JD.slope", str(self.sed.jd.slope)]
                if self.sed.jd.rmin is not None:
                    sc += ["-SED.JD.rmin", str(self.sed.jd.rmin)]
                if self.sed.jd.rmax is not None:
                    sc += ["-SED.JD.rmax", str(self.sed.jd.rmax)]
                sc += ["-SED.JD.MRwa", str(self.sed.jd.mrwa)]
                sc += ["-SED.JD.MIwa", str(self.sed.jd.miwa)]
                sc += ["-SED.JD.rate", str(self.sed.jd.rate)]
            #     Secondary LND mode :
            if self.sed.sm is not None:
                sc += ["-SED.LND.SM.SDradius", str(self.sed.sm.sdradius)]
                sc += ["-SED.LND.SM.SDvar", str(self.sed.sm.sdvar)]
                sc += ["-SED.LND.SM.MRwa", str(self.sed.sm.mrwa)]
                sc += ["-SED.LND.SM.MIwa", str(self.sed.sm.miwa)]
                sc += ["-SED.LND.SM.rate", str(self.sed.sm.rate)]
            #     Tertiary LND mode :
            if self.sed.tm is not None:
                sc += ["-SED.LND.TM.SDradius", str(self.sed.tm.sdradius)]
                sc += ["-SED.LND.TM.SDvar", str(self.sed.tm.sdvar)]
                sc += ["-SED.LND.TM.MRwa", str(self.sed.tm.mrwa)]
                sc += ["-SED.LND.TM.MIwa", str(self.sed.tm.miwa)]
                sc += ["-SED.LND.TM.rate", str(self.sed.tm.rate)]
        #     Hydrosols parameters for external data (phase functions,
        #     scattering and extinction coefficients) :
        if self.hyd.model == 2:
            sc += ["-HYD.ExtData", str(self.hyd.extdata)]
        #
        #   Sea / atmosphere interface parameters :
        #   --------------------------------------
        if self.log.sea is not None:
            sc += ["-SEA.Log", str(self.log.sea)]
        sc += ["-SEA.Dir", str(self.dirmie.sea)]
        sc += ["-SEA.Ind", str(self.sea.ind)]

        sc += ["-SEA.Wind", str(self.sea.wind)]

        key = hashlib.sha1(shlex.join(sc).encode()).hexdigest()

        sc += ["-OSOAA.ResRoot", str(resroot)]

        return key, sc

    def run(self, root=None, forcerun=False, fatm_null=False):
        """Run OSOAA. If no root directory is given for OSOAA the one
        configured by the system is used.

        If forcerun is set to true the simulation will be run even
        if it exists.

        fatm_null is True, run OSOAA_MAIN_FATM_NULL.exe.
        """

        if root is not None:
            self.root = root

        if self.root is None:
            raise Exception(...)

        # XXX
        if not os.path.exists(self.dirmie.aer):
            os.makedirs(self.dirmie.aer)
        if not os.path.exists(self.dirmie.hyd):
            os.makedirs(self.dirmie.hyd)
        if not os.path.exists(self.dirmie.sea):
            os.makedirs(self.dirmie.sea)

        with TemporaryDirectory() as temp_dir:
            key, sc = self._bake_arguments(resroot=temp_dir, fatm_null=fatm_null)
            # print(sc)

            perm_dir = self.resroot or os.path.join(self.root, "results", key)
            if os.path.exists(perm_dir):
                if not forcerun:
                    self.outputs = OUTPUTS(perm_dir, self.results)
                    return
                else:
                    shutil.rmtree(perm_dir)

            r = subprocess.run(sc, cwd=temp_dir, capture_output=True)

            if r.returncode or r.stderr:
                raise Exception(f'OSOAA finished with code {r.returncode}: {r.stderr} {r.stdout}')

            if self.logfile is not None:
                with open(self.logfile, 'wb') as logfile:
                    logfile.write(r.stdout)

            self.outputs = OUTPUTS(temp_dir, self.results)

            if not self.cleanup:
                shutil.copytree(temp_dir, perm_dir)


def test():
    s = OSOAA()
    s.run()
    print("OSOAA wrapper script by Francisco Nemiña")
    print("Inspired by Py6S wrapper by Robin Wilson")
    print("Using OSOAA located at {}".format(s.root))
    print("Running OSOAA using a set of test parameters")
    print("The results are:")
    print("Expected result: 0.128266")
    print("Actual result: {}".format(s.outputs.vsvza.I[51]))
    if s.outputs.vsvza.I[51] == 0.128266:
        print("#### Results agree PyOSOAA is working correctly")
    if s.outputs.vsvza.I[51] != 0.128266:
        print("#### Results do not agree PyOSOAA is not working correctly")
