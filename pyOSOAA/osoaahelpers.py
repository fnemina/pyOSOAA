# coding=utf-8

import numpy as np
from scipy.interpolate import interp1d


def ConfigureOcean(s, ocean_type="black"):
    """ This function defines custom ocean to use with the OSOAA. The
        The functions returns a reconfigured pyOSOAA object.
        Parameters
        ----------
        s               The pyOSOAA object for which we want to set the ocean
        ocean_type      The ocean type
                        black   A black ocean where water leaving radiance is
                                zero
        """

    if ocean_type not in ["black"]:
        raise(ValueError("Wrong ocean type."))

    if ocean_type is "black":
        # Sea bottom configuration
        s.sea.depth = 0.05
        s.sea.bottype = 1
        s.sea.botalb = 0

        # Sea particles configuration
        s.phyto.chl = 0
        s.sed.csed = 0
        s.det.abs440 = 1e5
        s.ys.abs440 = 1e5

        return s

def RunWavelengths(s, wavelengths=[0.550], angle=0, output="I", tau=False):
    """ This method run the simulation for a given pyOSOAA object for a set of
        wavelengths and angles and returns the output from the file vsVZA

        Parameters
        ----------
        s           The pyOSOAA object for which we want to run the simulation.
        wavelengths An interable with the wavelengths in micrometers for which
                    to run the simulations
        angles      The view angle in degrees for which to run the simulation.
        output      The name of the output we want to compute as an string
                    I           Stokes parameter at output level Z (in sr-1)
                                normalised to the extraterrestrial solar
                                irradiance
                                (PI * L(z) / Esun)
                    refl        Reflectance at output level Z
                    (PI * L(z) / Ed(z))
                    polrate     Degree of polarization (%)
                    lpol        Polarized intensity at output level Z (in sr-1)
                                normalised to the extraterrestrial solar
                                irradiance
                                (PI * Lpol(z) / Esun)
                    reflpol     Polarized reflectance at output level Z
                                (PI * Lpol(z) / Ed(z))
                    tau         True to compute optical thickness
        """

    if output not in ["I", "refl", "polrate", "lpol", "reflpol"]:
        raise(ValueError("Wrong output variable."))

    values = np.array([])
    tauv = np.array([])

    if type(angle) is int or type(angle) is float:
        angle = np.zeros(np.size(wavelengths))+angle

    for idx, wl in np.ndenumerate(wavelengths):
        # We set the wavelength and run the simulation
        s.wa = wl
        s.run()
        # Convert the output to a directory
        results = vars(s.outputs.vsvza)
        # We interpolte the values and add it to a numpy array
        f = interp1d(results['vza'], results[output])
        values = np.append(values, f(angle[idx[0]]))
        tauv = np.append(tauv, s.outputs.profileatm.tau[-1])

    if tau is True:
        return values, tauv

    return values


def RunAngles(s, wavelength=0.550, thetav=0, thetas=40, phi=90,
              output="I"):
    """ This method run the simulation for a given pyOSOAA object for a set of
        angles and a wavelength and returns the output from the file vsVZA

        Parameters
        ----------
        s           The pyOSOAA object for which we want to run the simulation.
        wavelengths An interable with the wavelengths in micrometers for which
                    to run the simulations
        thetav      The view angle in degrees for which to run the simulation.
        thetas      The sun angle in degrees for which to run the simulation.
        phi         The relative acimuth angle in degrees for which to run the
                    simulation.
        output      The name of the output we want to compute as an string
                    I           Stokes parameter at output level Z (in sr-1)
                                normalised to the extraterrestrial solar
                                irradiance
                                (PI * L(z) / Esun)
                    refl        Reflectance at output level Z
                    (PI * L(z) / Ed(z))
                    polrate     Degree of polarization (%)
                    lpol        Polarized intensity at output level Z (in sr-1)
                                normalised to the extraterrestrial solar
                                irradiance
                                (PI * Lpol(z) / Esun)
                    reflpol     Polarized reflectance at output level Z
                                (PI * Lpol(z) / Ed(z))
        """

    if output not in ["I", "refl", "polrate", "lpol", "reflpol"]:
        raise(ValueError("Wrong output variable."))

    # We check the type for the angles
    if type(thetav) is int or type(thetav) is float:
        thetav = np.array([thetav])
    if type(thetas) is int or type(thetas) is float:
        thetas = np.array([thetas])
    if type(phi) is int or type(phi) is float:
        phi = np.array([phi])

    values = -np.ones((thetav.size, thetas.size, phi.size))

    # We set the wavelength and run the simulation
    s.wa = wavelength

    for idxthetas, valthetas in np.ndenumerate(thetas):
        for idxphi, valphi in np.ndenumerate(phi):
            s.view.phi = valphi
            s.ang.thetas = valthetas
            s.run()
            # Convert the output to a directory
            results = vars(s.outputs.vsvza)
            # We interpolte the values and add it to a numpy array
            f = interp1d(results['vza'], results[output])
            # We reshape the array to make the asignement
            tempval = f(thetav).reshape(thetav.size, 1)
            values[:, idxthetas, idxphi] = tempval[:]

    return values
