import numpy as np
import matplotlib.pyplot as plt
import pyOSOAA
from pyOSOAA.osoaahelpers import ConfigureOcean


def getLglint():
    s = pyOSOAA.OSOAA()

    # simulate wavelength [um]
    s.wa = 0.56

    # relative azimuth angle
    s.view.phi = 0

    # aerosol optical depth
    s.aer.aotref = 0.2

    # solar zenith angle
    s.ang.thetas = 30

    # wind speed [m s^-1]
    s.sea.wind = 5

    # over the surface 0+
    s.view.level = 3

    s = ConfigureOcean(s)
    s.run()

    Lglint = s.outputs.vsvza.I

    return s.outputs.vsvza.vza, Lglint


def getLt():
    s = pyOSOAA.OSOAA()

    # simulate wavelength [um]
    s.wa = 0.56

    # relative azimuth angle
    s.view.phi = 0

    # aerosol optical depth
    s.aer.aotref = 0.2

    # solar zenith angle
    s.ang.thetas = 30

    # wind speed [m s^-1]
    s.sea.wind = 5

    # view at TOA
    s.view.level = 1

    # seabed depth [m]
    s.sea.depth = 50

    # seabed composition [Sand]
    s.sea.bottype = 2

    # Chl-a concentration [mg m^-3]
    s.phyto.chl = 1.8

    # refractive index for phytoplankton-like particles
    s.phyto.jd.mrwa = 1.05

    # SPM concentration [g m^-3]
    s.sed.csed = 9.0

    # refractive index for SPM particles
    s.sed.jd.mrwa = 1.15

    # a_CDOM(443nm) [m^-1]
    s.ys.abs440 = 0.07

    # spectral slope of a_CDOM
    s.ys.swa = 0.0176

    s.run()

    return s.outputs.vsvza.vza, s.outputs.vsvza.I


def getLsur():
    s = pyOSOAA.OSOAA()

    # simulate wavelength [um]
    s.wa = 0.56

    # relative azimuth angle
    s.view.phi = 0

    # aerosol optical depth
    s.aer.aotref = 0.2

    # solar zenith angle
    s.ang.thetas = 30

    # wind speed [m s^-1]
    s.sea.wind = 5

    # over the surface 0+
    s.view.level = 3

    # seabed depth [m]
    s.sea.depth = 50

    # seabed composition [Sand]
    s.sea.bottype = 2

    # Chl-a concentration [mg m^-3]
    s.phyto.chl = 1.8

    # refractive index for phytoplankton-like particles
    s.phyto.jd.mrwa = 1.05

    # SPM concentration [g m^-3]
    s.sed.csed = 9.0

    # refractive index for SPM particles
    s.sed.jd.mrwa = 1.15

    # a_CDOM(443nm) [m^-1]
    s.ys.abs440 = 0.07

    # spectral slope of a_CDOM
    s.ys.swa = 0.0176

    s.run()
    return s.outputs.vsvza.vza, s.outputs.vsvza.I


if __name__ == "__main__":
    vza, Lglint = getLglint()
    vza, Lsur = getLsur()
    Lw = Lsur - Lglint
    vza, Lt = getLt()
    Es = 1767.558   # W m^-2
    fig, ax = plt.subplots()
    ax.plot(vza, Lw *
            Es/np.pi, color='g', label=r'$L_{W}$')
    ax.plot(vza, Lglint*Es/np.pi,
            color='b', label=r'$L_{GLINT}$')
    ax.plot(vza, Lt*Es/np.pi,
            color='r', label=r'$L_{TOA}$')
    ax.set_xlabel(r'$\theta_{v}[^{\circ}]$')
    ax.set_ylabel(r'Radiance L $[W\cdot m^{-2}\cdot um^{-1}\cdot sr^{-1}]$')
    ax.set_xlim(-75, 75)
    plt.legend()
    plt.savefig('figure1a.png', dpi=300, bbox_inches='tight')
    plt.close()

    fig, ax = plt.subplots()
    ax.plot(vza, Lw/Lt*100, color='g')

    ax.set_xlabel(r'$\theta_{v}[^{\circ}]$')
    ax.set_ylabel(r'$L_{w}/L_{TOA}$[%]')
    ax.set_xlim(-75, 75)
    ax.set_ylim(0, 100)
    plt.savefig('figure2a.png', dpi=300, bbox_inches='tight')
    plt.close()
