import numpy as np
from scipy.io import FortranFile


def ExtractValue(text, reference):
    """ This function extracts the value for a given reference string for the
        text string given and returns it as a float number
        text        Text to search
        reference   Reference string
        """
    for a in text:
        if reference in a:
            res = a.replace(reference, "")

    return float(res)


class VSVZA(object):
    """ This file provides the upwelling radiance field (i.e., the Stokes
        paramaters I,Q,U where I is the radiance) versus the viewing zenith
        angle, for the given relative azimuth angle (-OSOAA.View.Phi) and
        for the given altitude or depth (-OSOAA.View.Z associated to
        -OSOAA.View.Level equals 5).
        This ascii file is composed of a header which describes in detail
        the structure of the file and columns data.
        """

    def __init__(self, resroot, filename="LUM_vsVZA.txt"):
        """ Filename defined by -OSOAA.ResFile.vsVZA, if none is defined
            then LUM_vsVZA.txt is used:
            This file provides the upwelling radiance field (i.e., the
            Stokes paramaters I,Q,U where I is the radiance) versus the
            viewing zenith angle, for the given relative azimuth angle
            (-OSOAA.View.Phi) and for the given altitude or depth
            (-OSOAA.View.Z associated to -OSOAA.View.Level equals 5).
            This ascii file is composed of a header which describes in
            detail the structure of the file and columns data.

            resroot     OSOAA results root directory
            filename    Filename to look for the results.
            fulltext    Full file text.
            vzaless     Simulated relative azimuth (degrees) for VZA < 0
                        (sign convention)
            vzamore     Simulated relative azimuth (degrees) for VZA > 0
                        (sign convention)
            depth       Value of the depth selected for the output (m)
            vza         Viewing Zenith Angle (deg)
            scaang      Scattering angle (deg)
            I           Stokes parameter at output level Z (in sr-1)
                        normalised to the extraterrestrial solar irradiance
                        (PI * L(z) / Esun)
            refl        Reflectance at output level Z (PI * L(z) / Ed(z))
            polrate     Degree of polarization (%)
            lpol        Polarized intensity at output level Z (in sr-1)
                        normalised to the extraterrestrial solar irradiance
                        (PI * Lpol(z) / Esun)
            reflpol     Polarized reflectance at output level Z
                        (PI * Lpol(z) / Ed(z))
            """

        # We open the file with the corresponding encoding and convert it
        # to a text string.
        if filename is None:
            filename = "LUM_vsVZA.txt"
        with open(resroot+"/Standard_outputs/"+filename,
                  encoding="iso-8859-15") as file:
            self.fulltext = file.readlines()
        # Read variables not tabulated
        self.vzaless = ExtractValue(self.fulltext,
                                    "for VZA < 0 (sign convention):")
        self.vzamore = ExtractValue(self.fulltext,
                                    "for VZA > 0 (sign convention):")
        self.depth = ExtractValue(self.fulltext,
                                  "Value of the depth selected for the "
                                  + "output (m) :")

        # Get header length to skip it
        skipheader = [idx for idx, text in enumerate(self.fulltext)
                      if "VZA    SCA_ANG" in text][0]
        self.vza, self.scaang, self.I, self.refl, self.polrate, self.lpol,\
        self.reflpol = np.genfromtxt(resroot+"/Standard_outputs/"+filename,
                                     skip_header=skipheader+1, unpack=True,
                                     encoding="iso-8859-15")


class VSZ(object):
    """ Filename defined by -OSOAA.ResFile.vsZ:
        This optional file provides the radiance field (I,Q,U) versus the
        altitude or depth, for the given relative azimuth angle
        (-OSOAA.View.Phi) and for the given viewing zenith angle
        (-OSOAA.View.VZA).
        """

    def __init__(self, resroot, filename=None):
        """ Filename defined by -OSOAA.ResFile.vsZ:
            This optional file provides the radiance field (I,Q,U) versus
            the altitude or depth, for the given relative azimuth angle
            (-OSOAA.View.Phi) and for the given viewing zenith angle
            (-OSOAA.View.VZA).
            resroot     OSOAA results root directory.
            filename    Filename to look for the results.
            fulltext    Full file text.
            z           Depth in the sea (meters)
            scaang      Scattering angle (deg)
            I           Stokes parameter at level Z (in sr-1) normalised to the
                        extraterrestrial solar irradiance
                        (PI * L(z) / Esun)
            refl        Reflectance at level Z
                        (PI * L(z) / Ed(z))
            polrate     Degree of polarization (%)
            lpol        Polarized intensity at level Z (in sr-1) normalised to
                        the extraterrestrial solar  irradiance
                        (PI * Lpol(z) / Esun)
            reflpol     Polarized reflectance at level Z
                        (PI * Lpol(z) / Ed(z))
            """

        with open(resroot+"/Standard_outputs/"+filename,
                  encoding="iso-8859-15") as file:
            self.fulltext = file.readlines()
        # Read variables not tabulated
        self.simazimuth = ExtractValue(self.fulltext,
                                       "Simulated relative azimuth (degrees) :")
        self.updirecion = ExtractValue(self.fulltext,
                                       "Upward direction VZA (degrees) :")

        # Get header length to skip it
        skipheader = [idx for idx, text in enumerate(self.fulltext)
                      if "       Z     SCA_ANG" in text][0]
        self.z, self.scaang, self.I, self.refl, self.polrate,\
        self.lpol, self.reflpol = np.genfromtxt(resroot+"/Standard_outputs/"+filename,
                                                skip_header=skipheader+1,
                                                unpack=True,
                                                encoding="iso-8859-15")


class BIN(object):
    """ The file defined by -OSOAA.ResFile.Bin is a binary file containing the
        (I,Q,U) Stoke’s parameters in term of Fourier series expansion. This
        file is NOT devoted to be used by the user; it is an intermediary file
        for the computations.
        If not defined the default filename in OSOAA.h is applied:
        CTE_DEFAULT_FICSOS_RES_BIN .
        """

    def __init__(self, resroot, filename="LUM_SF.bin", NT_TOT=107, NBMU=51):
        """ Filename defined by -OSOAA.ResFile.vsZ:
            This optional file provides the radiance field (I,Q,U) versus
            the altitude or depth, for the given relative azimuth angle
            (-OSOAA.View.Phi) and for the given viewing zenith angle
            (-OSOAA.View.VZA).
            resroot     OSOAA results root directory.
            filename    Filename to look for the results.
            I           I Stokes parameters Fourier components. The array
                        is I(LEVEL;WAVE_NUMBER)
            Q           Q Stokes parameters Fourier components. The array
                        is Q(LEVEL;WAVE_NUMBER)
            U           U Stokes parameters Fourier components. The array
                        is U(LEVEL;WAVE_NUMBER)
            """
        # We read the fortran binary file
        with FortranFile(resroot+"/Advanced_outputs/"+filename, "r") as file:
            tmp = file.read_reals()

        # Reshape the array into 3 array for I, Q and U
        tmp = tmp.reshape(3, int(tmp.size/3))
        Itmp = tmp[0, :]
        Qtmp = tmp[1, :]
        Utmp = tmp[2, :]

        # We reshape the array for the number of levels
        Itmp = Itmp.reshape(int(Itmp.size/(NT_TOT+1)), NT_TOT+1)
        Qtmp = Qtmp.reshape(int(Qtmp.size/(NT_TOT+1)), NT_TOT+1)
        Utmp = Utmp.reshape(int(Utmp.size/(NT_TOT+1)), NT_TOT+1)

        # We transpose the result to get and array of ARRAY[LEVEL:ANGLE]
        self.I = Itmp.transpose()
        self.Q = Qtmp.transpose()
        self.U = Utmp.transpose()


class ADVUPDOWN(object):
    """ This file provides the upwelling (respectively downwelling) Stoke’s
        parameters (I,Q,U) versus the viewing zenith angle and versus the
        complete profile (depth or altitude), for the given relative azimuth
        angle (-OSOAA.View.Phi).
        For each level of the profile, the radiance field (I,Q,U) is provided
        for each viewing zenith angle. The scattering angle, the polarization
        angle and the degree of polarization (also called polarization rate),
        as well as the polarized intensity are also provided.
        """

    def __init__(self, resroot, filename):
        """ Filename defined by -OSOAA.ResFile.Adv.Up/Down
            This file provides the upwelling (respectively downwelling) Soke’s
            parameters (I,Q,U) versus the viewing zenith angle and versus the
            complete profile (depth or altitude), for the given relative
            azimuth angle (-OSOAA.View.Phi).

            For each level of the profile, the radiance field (I,Q,U) is
            provided for each viewing zenith angle. The scattering angle, the
            polarization angle and the degree of polarization (also called
            polarization rate), as well as the polarized intensity are also
            provided.

            resroot     OSOAA results root directory.
            filename    Filename to look for the results.
            fulltext    Full file text.
            vzaless     Simulated relative azimuth (degrees) for VZA < 0
                        (sign convention)
            vzamore     Simulated relative azimuth (degrees) for VZA > 0
                        (sign convention)
            level       Level I of the profile
                            TOA is level 0
                            0+ (over surface) is level           26
                            0- (under surface) is level          27
            z           Altitude or depth (meters)
            vza         Viewing Zenith Angle (deg)
            scaang      Scattering angle (deg)
            I, Q, U     Stokes parameters at level Z (in sr-1)
                        normalised to the extraterrestrial solar irradiance
                        (PI * L(z) / Esun)
            polang      Polarization angle (deg)
            polrate     Degree of polarization (%)
            lpol        Polarized intensity at level Z (in sr-1)
                        normalised to the extraterrestrial solar irradiance
                        (PI * Lpol(z) / Esun)

            """

        # We open the file with the corresponding encoding and convert it
        # to a text string.
        with open(resroot+"/Advanced_outputs/"+filename,
                  encoding="iso-8859-15") as file:
            self.fulltext = file.readlines()
        # Read variables not tabulated
        self.vzaless = ExtractValue(self.fulltext,
                                    "for VZA < 0 (sign convention):")
        self.vzamore = ExtractValue(self.fulltext,
                                    "for VZA > 0 (sign convention):")

        # Get header length to skip it
        skipheader = [idx for idx, text in enumerate(self.fulltext)
                      if "  LEVEL    Z         VZA" in text][0]
        self.level, self.z, self.vza, self.scaang, self.I, self.Q, self.U,\
        self.polang, self.polrate, self.lpol = np.genfromtxt(resroot+"/Advanced_outputs/"+filename,
                                                             skip_header=skipheader+1, unpack=True,
                                                             encoding="iso-8859-15")


class PROFILE_SEA(object):
    """ This file contains the sea optical thickness vertical profile
        calculated by OSOAA.
        """

    def __init__(self, resroot, filename="PROFILE_SEA.txt"):
        """ This file contains the sea optical thickness vertical profile
            calculated by OSOAA.

            resroot     OSOAA results root directory.
            filename    Filename to look for the results.
            fulltext    Full file text.
            level       The level number (or index) p
            z           The depth of level p (in meters)
            tau         The extinction optical thickness for the level p,
            mixmol      The scattering mixing rate of pure sea water for
                        the layer p
            mixphy      The scattering mixing rate of phytoplankton for the
                        layer p
            mixpc       The scattering mixing rate of Mineral-Like
                        particles for the layer p
            Note that the radiative parameters define the real profile
            without any adjustment to the phase function truncation of the
            components.
            """

        # Get header length to skip it
        skipheader = 10
        self.level, self.z, self.tau, self.mixmol, self.mixphy,\
        self.mixpc = np.genfromtxt(resroot+"/Advanced_outputs/"+filename,
                                   skip_header=skipheader+1, unpack=True,
                                   encoding="iso-8859-15")


class PROFILE_ATM(object):
    """ This file contains the atmospheric optical thickness vertical profile
        calculated by OSOAA.
        """

    def __init__(self, resroot, filename="PROFILE_ATM.txt"):
        """ This file contains the atmospheric optical thickness vertical
        profile calculated by OSOAA.

            resroot     OSOAA results root directory.
            filename    Filename to look for the results.
            fulltext    Full file text.
            level       Level I of the atmospheric profile
            z           Altitude of the level I (in km)
            tau         Total extinction opt. thickness at level I
            mixaer      Mixing rate aer. for scattering in layer I
            mixray      Mixing rate mol. for scattering in layer I

            The radiative parameters define the real profile without any
            adjustment to the aerosol phase function truncation.
            """

        # Get header length to skip it
        skipheader = 9
        self.level, self.z, self.tau, self.mixaer,\
        self.mixray = np.genfromtxt(resroot+"/Advanced_outputs/"+filename,
                                   skip_header=skipheader+1, unpack=True,
                                   encoding="iso-8859-15")


class PM(object):
    """ Files containing the radiative properties of Aerosols, chlorophyll or
        Mineral-Like particles
        """

    def __init__(self, resroot, filename):
        """ First thirteen header lines provide comments and formatted data on
            the extinction and scattering cross-sections (in μm 2 ), the
            asymmetry factor, the volume of the equivalent mean particle
            (in μm 3 ), the real part of the mean refractive index, the phase
            function truncation coefficient and the single scattering albedo
            (adjusted to the truncation).

            The following lines contain the phase matrix coefficients of the
            development of the Legendre Polynomials of the phase matrix for
            each order k ranging from k= 0 up to the maximum order of
            computations (OS_NB).

            resroot     OSOAA results root directory.
            filename    Filename to look for the results.
            fulltext    Full file text.
            extcs       Extinction cross section (mic^2)
            scacs       Scattering cross section (mic^2)
            asymm       Asymmetry factor (no truncation)
            mpd         Mean particules altitude/depth (m)
            vol         Volume of the mean particule (mic^3)
            rindex      Mean refractive index (real part)
            trunca      Truncation coefficient
            singlesca   Single scattering albedo (truncation)

            alpha       the coefficient alpha Related to the polarized phase
                        functions
            beta        the coefficient beta Related to the polarized phase
                        functions
            gamma       the coefficient gamma Related to the polarized phase
                        functions
            xi          the coefficient xi Related to the polarized phase
                        functions

            These coefficients are adjusted to a phase function truncation if
            applied.
            """
        # We open the file with the corresponding encoding and convert it
        # to a text string.
        with open(resroot+"/Advanced_outputs/"+filename,
                  encoding="iso-8859-15") as file:
            self.fulltext = file.readlines()
        # Read variables not tabulated
        self.extcs = ExtractValue(self.fulltext,
                                    "EXTINCTION CROSS SECTION (mic^2)     :")
        self.scacs = ExtractValue(self.fulltext,
                                    "SCATTERING CROSS SECTION (mic^2)     :")
        self.asymm = ExtractValue(self.fulltext,
                                    "ASYMMETRY FACTOR (no truncation)     :")
        self.mpd = ExtractValue(self.fulltext,
                                    "MEAN PARTICULES ALTITUDE/DEPTH (m)   :")
        self.vol = ExtractValue(self.fulltext,
                                    "VOLUME OF A MEAN PARTICULE (mic^3)   :")
        self.rindex = ExtractValue(self.fulltext,
                                    "MEAN REFRACTIVE INDEX (real part)    :")
        self.trunca = ExtractValue(self.fulltext,
                                    "TRUNCATION COEFFICIENT               :")
        self.singlesca = ExtractValue(self.fulltext,
                                    "SINGLE SCATTERING ALBEDO (truncation): ")

        # Get header length to skip it
        skipheader = [idx for idx, text in enumerate(self.fulltext)
                      if "ALPHA(K)        BETA11(K)" in text][0]
        self.alpha, self.beta, self.gamma, self.xi = np.genfromtxt(resroot+"/Advanced_outputs/"+filename,
                                                                   skip_header=skipheader+1, unpack=True,
                                                                   encoding="iso-8859-15")


class FLUX(object):
    """ 'Flux.txt' file containing the profile of downward and upward fluxes
        from TOA to the sea bottom, normalised to the solar irradiance at TOA.
        """

    def __init__(self, resroot, filename="Flux.txt"):
        """ First five header lines provide comments.
            The following lines contain the values of fluxes.
            For each line, this file provides:
            levels      The level number (or index) k
            z           The altitude/depth of level k (in m): positive value in
                        the atmosphere, negative value in the sea
            Eddir       The direct downward flux
            Eddif       The diffuse downward flux
            Ed          The total downward flux (Ed)
            Eudir       The direct upward flux
            Eudif       The diffuse upward flux
            Eu          The total upward flux (Eu)
            ratio       The ratio Eu/Ed: total upward flux / total downward
                        flux
            The fluxes are normalised to the solar irradiance at TOA (no unit).
            """
        # Get header length to skip it
        skipheader = 4
        self.level, self.z, self.Eddir, self.Eddif, self.Ed, self.Eudir,\
        self.Eudif, self.Eu, self.ratio = np.genfromtxt(resroot+"/Advanced_outputs/"+filename,
                                                        skip_header=skipheader+1, unpack=True,
                                                        encoding="iso-8859-15")


class OUTPUTS(object):
    """ This class contains the standard and advanced outputs generated by the
        OSOAA software"""

    def __init__(self, resroot, filenames):
        """ This methods inits the output class with all the avalaible outputs.
            resroot     Results root Directory
            filenames   Object with all filenames
            """

        self.vsvza = VSVZA(resroot, filenames.vsvza)

        if filenames.vsz is not None:
            self.vsz = VSZ(resroot, filenames.vsz)

        if filenames.advup is not None:
            self.advup = ADVUPDOWN(resroot, filenames.advup)

        if filenames.advdown is not None:
            self.advdown = ADVUPDOWN(resroot, filenames.advdown)

        if filenames.profilesea is None:
            self.profilesea = PROFILE_SEA(resroot)
        else:
            self.profilesea = PROFILE_SEA(resroot, filenames.profilesea)

        if filenames.profileatm is None:
            self.profileatm = PROFILE_ATM(resroot)
        else:
            self.profileatm = PROFILE_ATM(resroot, filenames.profileatm)

        if filenames.aer is None:
            self.aer = PM(resroot, filename="PM_AER.txt")
        else:
            self.aer = PM(resroot, filenames.aer)

        if filenames.phyto is None:
            self.phyto = PM(resroot, filename="PM_PHYTO.txt")
        else:
            self.phyto = PM(resroot, filenames.phyto)

        if filenames.mlp is None:
            self.mlp = PM(resroot, filename="PM_MLP.txt")
        else:
            self.mlp = PM(resroot, filenames.mlp)

        self.flux = FLUX(resroot)

        self.bin = BIN(resroot)
