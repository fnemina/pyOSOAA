import numpy as np
from scipy.interpolate import interp1d


def RunWavelengths(s, wavelengths=[0.550], angle=0, output="I"):
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
        """

    if output not in ["I", "refl", "polrate", "lpol", "reflpol"]:
        raise("Not a valid output parameter.")
        return -1

    values = np.array([])

    if angle.type is int or angle.type is float:
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

    return values
