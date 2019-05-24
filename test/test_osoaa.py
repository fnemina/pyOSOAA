# coding=utf-8

import unittest
import os
import pyOSOAA
import numpy as np


class TestOUTPUTSClasses(unittest.TestCase):
    def testExtractValue(self):
        from pyOSOAA.outputs import ExtractValue
        string = ["Zero 0", "One 1"]
        self.assertEqual(ExtractValue(string, "Zero"), 0)
        self.assertEqual(ExtractValue(string, "One"), 1)

    def testVSVZA(self):
        resroot = os.path.dirname(os.path.abspath(__file__))+"/OSOAA_RESULTS"
        vsvza = pyOSOAA.outputs.VSVZA(resroot, "RESLUM_vsVZA.txt")
        self.assertEqual(vsvza.vzaless, 180.0)
        self.assertEqual(vsvza.vzamore, 0.0)
        self.assertEqual(vsvza.vza[0], -89.07)
        self.assertEqual(vsvza.vza.size, 102)
        self.assertEqual(vsvza.scaang[0], 112.84)
        self.assertEqual(vsvza.scaang.size, 102)
        self.assertEqual(vsvza.I[0], 0.122260E-04)
        self.assertEqual(vsvza.I.size, 102)
        self.assertEqual(vsvza.refl[0], 0.139485E-02)
        self.assertEqual(vsvza.refl.size, 102)
        self.assertEqual(vsvza.polrate[0], 66.45)
        self.assertEqual(vsvza.polrate.size, 102)
        self.assertEqual(vsvza.lpol[0], 0.812479E-05)
        self.assertEqual(vsvza.lpol.size, 102)
        self.assertEqual(vsvza.reflpol[0], 0.926945E-03)
        self.assertEqual(vsvza.reflpol.size, 102)

    def testVSZ(self):
        resroot = os.path.dirname(os.path.abspath(__file__))+"/OSOAA_RESULTS"
        vsz = pyOSOAA.outputs.VSZ(resroot, "RESLUM_vsZ.txt")
        self.assertEqual(vsz.simazimuth, 0.0)
        self.assertEqual(vsz.updirecion, 0.0)
        self.assertEqual(vsz.z[0], -0)
        self.assertEqual(vsz.z.size, 81)
        self.assertEqual(vsz.scaang[0], 158.09)
        self.assertEqual(vsz.scaang.size, 81)
        self.assertEqual(vsz.I[0], 0.621761E-03)
        self.assertEqual(vsz.I.size, 81)
        self.assertEqual(vsz.refl[0], 0.760278E-03)
        self.assertEqual(vsz.refl.size, 81)
        self.assertEqual(vsz.polrate[0], 6.45)
        self.assertEqual(vsz.polrate.size, 81)
        self.assertEqual(vsz.lpol[0], 0.401063E-04)
        self.assertEqual(vsz.lpol.size, 81)
        self.assertEqual(vsz.reflpol[0], 0.490412E-04)
        self.assertEqual(vsz.reflpol.size, 81)

    def testBIN(self):
        # This code check the BIN reading method.
        # To test it the code LUM_SF.f95 is used to generate a matrix and which
        # Is then read by this code.
        # By the way FOTRAN saves it matrices there is a missmatch for indices
        # We solve that by computing k-NBMU here and k in the fortran code.
        # Comparison is done for an error less than 0.001 that should be enough

        resroot = os.path.dirname(os.path.abspath(__file__))+"/OSOAA_RESULTS"
        NB_TOT = 107
        NBMU = 51
        I3 = np.zeros((NB_TOT+1, 2*NBMU+1))
        Q3 = np.zeros((NB_TOT+1, 2*NBMU+1))
        U3 = np.zeros((NB_TOT+1, 2*NBMU+1))

        for j in range(0, NB_TOT+1):
            for k in range(0, 2*NBMU+1):
                I3[j, k] = 0+k-NBMU+j/108.
                Q3[j, k] = 1000+k-NBMU+j/108.
                U3[j, k] = 2000+k-NBMU+j/108.
        bin = pyOSOAA.outputs.BIN(resroot)

        self.assertTrue((bin.I-Q3 < 1e-3).all())
        self.assertTrue((bin.Q-Q3 < 1e-3).all())
        self.assertTrue((bin.U-U3 < 1e-3).all())

    def testADVUPDOWN(self):
        resroot = os.path.dirname(os.path.abspath(__file__))+"/OSOAA_RESULTS"

        advup = pyOSOAA.outputs.ADVUPDOWN(resroot, "RESLUM_Advanced_UP.txt")
        self.assertEqual(advup.vzaless, 180.0)
        self.assertEqual(advup.vzamore, 0.0)
        self.assertEqual(advup.level[0], 0)
        self.assertEqual(advup.level.size, 11016)
        self.assertEqual(advup.z[0], 300000.0)
        self.assertEqual(advup.z.size, 11016)
        self.assertEqual(advup.vza[0], -89.07)
        self.assertEqual(advup.vza.size, 11016)
        self.assertEqual(advup.scaang[0], 120.93)
        self.assertEqual(advup.scaang.size, 11016)
        self.assertEqual(advup.I[0], 0.232266E+00)
        self.assertEqual(advup.I.size, 11016)
        self.assertEqual(advup.Q[0], -0.119399E+00)
        self.assertEqual(advup.Q.size, 11016)
        self.assertEqual(advup.U[0], 0.165581E-16)
        self.assertEqual(advup.U.size, 11016)
        self.assertEqual(advup.polang[0], 90)
        self.assertEqual(advup.polang.size, 11016)
        self.assertEqual(advup.polrate[0], 51.41)
        self.assertEqual(advup.polrate.size, 11016)
        self.assertEqual(advup.lpol[0], 0.119399E+00)
        self.assertEqual(advup.lpol.size, 11016)

        resroot = os.path.dirname(os.path.abspath(__file__))+"/OSOAA_RESULTS"

        advdown = pyOSOAA.outputs.ADVUPDOWN(resroot, "RESLUM_Advanced_DOWN.txt")
        self.assertEqual(advdown.vzaless, 180.0)
        self.assertEqual(advdown.vzamore, 0.0)
        self.assertEqual(advdown.level[0], 0)
        self.assertEqual(advdown.level.size, 11016)
        self.assertEqual(advdown.z[0], 300000.0)
        self.assertEqual(advdown.z.size, 11016)
        self.assertEqual(advdown.vza[0], -89.07)
        self.assertEqual(advdown.vza.size, 11016)
        self.assertEqual(advdown.scaang[0], 119.07)
        self.assertEqual(advdown.scaang.size, 11016)
        self.assertEqual(advdown.I[0], 0.0)
        self.assertEqual(advdown.I.size, 11016)
        self.assertEqual(advdown.Q[0], 0.0)
        self.assertEqual(advdown.Q.size, 11016)
        self.assertEqual(advdown.U[0], 0.0)
        self.assertEqual(advdown.U.size, 11016)
        self.assertEqual(advdown.polang[0], -999.0)
        self.assertEqual(advdown.polang.size, 11016)
        self.assertEqual(advdown.polrate[0], -999.0)
        self.assertEqual(advdown.polrate.size, 11016)
        self.assertEqual(advdown.lpol[0], 0.0)
        self.assertEqual(advdown.lpol.size, 11016)

    def testPROFILE_SEA(self):
        resroot = os.path.dirname(os.path.abspath(__file__))+"/OSOAA_RESULTS"
        profile = pyOSOAA.outputs.PROFILE_SEA(resroot)
        self.assertEqual(profile.level[0], 0)
        self.assertEqual(profile.level.size, 81)
        self.assertEqual(profile.z[0], 0.0)
        self.assertEqual(profile.z.size, 81)
        self.assertEqual(profile.tau[0], 0.0)
        self.assertEqual(profile.tau.size, 81)
        self.assertEqual(profile.mixmol[0], 0.00175)
        self.assertEqual(profile.mixmol.size, 81)
        self.assertEqual(profile.mixphy[0], 0.18193)
        self.assertEqual(profile.mixphy.size, 81)
        self.assertEqual(profile.mixpc[0], -0.0)
        self.assertEqual(profile.mixpc.size, 81)

    def testPROFILE_ATM(self):
        resroot = os.path.dirname(os.path.abspath(__file__))+"/OSOAA_RESULTS"
        profile = pyOSOAA.outputs.PROFILE_ATM(resroot)
        self.assertEqual(profile.level[0], 0)
        self.assertEqual(profile.level.size, 27)
        self.assertEqual(profile.z[0], 300.0)
        self.assertEqual(profile.z.size, 27)
        self.assertEqual(profile.tau[0], 0.0)
        self.assertEqual(profile.tau.size, 27)
        self.assertEqual(profile.mixaer[0], 0.0)
        self.assertEqual(profile.mixaer.size, 27)
        self.assertEqual(profile.mixray[0], 1.0)
        self.assertEqual(profile.mixray.size, 27)

    def testPM(self):
        resroot = os.path.dirname(os.path.abspath(__file__))+"/OSOAA_RESULTS"
        pm = pyOSOAA.outputs.PM(resroot, "PM_AER.txt")
        self.assertEqual(pm.extcs, 0.17015E+00)
        self.assertEqual(pm.scacs, 0.16985E+00)
        self.assertEqual(pm.asymm, 0.80875E+00)
        self.assertEqual(pm.mpd, -999.000)
        self.assertEqual(pm.vol, 0.41155E+00)
        self.assertEqual(pm.rindex, 1.3677)
        self.assertEqual(pm.trunca, 0.63643)
        self.assertEqual(pm.singlesca, 0.99741)
        self.assertEqual(pm.alpha[0], 0.0)
        self.assertEqual(pm.alpha.size, 81)
        self.assertEqual(pm.beta[0], 0.1E+01)
        self.assertEqual(pm.beta.size, 81)
        self.assertEqual(pm.gamma[0], 0.0)
        self.assertEqual(pm.gamma.size, 81)
        self.assertEqual(pm.xi[0], 0.0)
        self.assertEqual(pm.xi.size, 81)

    def testFLUX(self):
        resroot = os.path.dirname(os.path.abspath(__file__))+"/OSOAA_RESULTS"
        pm = pyOSOAA.outputs.FLUX(resroot, "Flux.txt")
        self.assertEqual(pm.level[0], 0)
        self.assertEqual(pm.level.size, 108)
        self.assertEqual(pm.z[0], 300000.0)
        self.assertEqual(pm.z.size, 108)
        self.assertEqual(pm.Eddir[0], 0.866025)
        self.assertEqual(pm.Eddir.size, 108)
        self.assertEqual(pm.Eddif[0], 0.0)
        self.assertEqual(pm.Eddif.size, 108)
        self.assertEqual(pm.Ed[0], 0.866025E+00)
        self.assertEqual(pm.Ed.size, 108)
        self.assertEqual(pm.Eudir[0], 0.151412E-01)
        self.assertEqual(pm.Eudir.size, 108)
        self.assertEqual(pm.Eudif[0], 0.344081E-01)
        self.assertEqual(pm.Eudif.size, 108)
        self.assertEqual(pm.Eu[0], 0.495493E-01)
        self.assertEqual(pm.Eu.size, 108)
        self.assertEqual(pm.ratio[0], 0.572147E-01)
        self.assertEqual(pm.ratio.size, 108)


class TestOSOAAClasses(unittest.TestCase):

    def testSEA(self):
        self.assertEqual(pyOSOAA.SEA().surfalb, 0.0)
        self.assertEqual(pyOSOAA.SEA().bottype, 1)
        self.assertEqual(pyOSOAA.SEA().botalb, 0.30)
        self.assertEqual(pyOSOAA.SEA().wind, 7)
        self.assertEqual(pyOSOAA.SEA().depth, 15.0)

    def testLOG(self):
        self.assertEqual(pyOSOAA.LOG().osoaa, "log_osoaa.txt")
        self.assertEqual(pyOSOAA.LOG().ang, "log_ang.txt")
        self.assertEqual(pyOSOAA.LOG().profile, "log_profile.txt")
        self.assertEqual(pyOSOAA.LOG().aer, "log_aer.txt")
        self.assertEqual(pyOSOAA.LOG().aermie, "log_aermie.txt")
        self.assertEqual(pyOSOAA.LOG().hyd, "log_hyd.txt")
        self.assertEqual(pyOSOAA.LOG().hydmie, "log_hydmie.txt")
        self.assertEqual(pyOSOAA.LOG().sea, "log_sea.txt")
        self.assertEqual(pyOSOAA.LOG().sos, "log_sos.txt")

    def testRESULTS(self):
        # We test default values
        self.assertIsNone(pyOSOAA.RESULTS().profileatm)
        self.assertIsNone(pyOSOAA.RESULTS().profilesea)
        self.assertIsNone(pyOSOAA.RESULTS().aer)
        self.assertIsNone(pyOSOAA.RESULTS().phyto)
        self.assertIsNone(pyOSOAA.RESULTS().mlp)
        self.assertIsNone(pyOSOAA.RESULTS().angrad)
        self.assertIsNone(pyOSOAA.RESULTS().angmie)
        self.assertIsNone(pyOSOAA.RESULTS().sosbin)
        self.assertIsNone(pyOSOAA.RESULTS().vsvza)
        self.assertEqual(pyOSOAA.RESULTS().advup, "resfile_advup.txt")
        self.assertEqual(pyOSOAA.RESULTS().advdown, "resfile_advdown.txt")
        self.assertEqual(pyOSOAA.RESULTS().vsz, "resfile_vsz.txt")
        # We test different Values
        s = pyOSOAA.OSOAA()
        s.results.vsz = "test_vsz.txt"
        s.results.vsvza = "test_vsvza.txt"

        s.results.profileatm = "test_profileatm.txt"
        s.results.profilesea = "test_profilesea.txt"
        s.results.aer = "test_aer.txt"
        s.results.phyto = "test_phyto.txt"
        s.results.mlp = "test_mlp.txt"
        s.results.angrad = "test_angrad.txt"
        s.results.angmie = "test_angmie.txt"
        s.results.sosbin = "test_sosbin"
        s.results.advup = "test_advup.txt"
        s.results.advdown = "test_advdown.txt"
        s.run()
        self.assertTrue(os.path.exists(s.resroot+"/Standard_outputs/"
                                       + "test_vsz.txt"))
        self.assertTrue(os.path.exists(s.resroot+"/Standard_outputs/"
                                       + "test_vsvza.txt"))

        self.assertTrue(os.path.exists(s.resroot+"/Advanced_outputs/"
                                       + "test_profileatm.txt"))
        self.assertTrue(os.path.exists(s.resroot+"/Advanced_outputs/"
                                       + "test_profilesea.txt"))
        self.assertTrue(os.path.exists(s.resroot+"/Advanced_outputs/"
                                       + "test_aer.txt"))
        self.assertTrue(os.path.exists(s.resroot+"/Advanced_outputs/"
                                       + "test_phyto.txt"))
        self.assertTrue(os.path.exists(s.resroot+"/Advanced_outputs/"
                                       + "test_mlp.txt"))
        self.assertTrue(os.path.exists(s.resroot+"/Advanced_outputs/"
                                       + "test_angrad.txt"))
        self.assertTrue(os.path.exists(s.resroot+"/Advanced_outputs/"
                                       + "test_sosbin"))
        self.assertTrue(os.path.exists(s.resroot+"/Advanced_outputs/"
                                       + "test_advup.txt"))
        self.assertTrue(os.path.exists(s.resroot+"/Advanced_outputs/"
                                       + "test_advdown.txt"))

    def testDIRMIE(self):
        self.assertEqual(pyOSOAA.DIRMIE("").hyd, "/DATABASE/MIE_HYD")
        self.assertEqual(pyOSOAA.DIRMIE("").aer, "/DATABASE/MIE_AER")
        self.assertEqual(pyOSOAA.DIRMIE("").sea, "/DATABASE/SURF_MATR")

    def testGP(self):
        self.assertEqual(pyOSOAA.GP(1, 5, 2).chlbg, 1)
        self.assertEqual(pyOSOAA.GP(1, 5, 2).deep, 5)
        self.assertEqual(pyOSOAA.GP(1, 5, 2).width, 2)

    def testPHYTO(self):
        phyto = pyOSOAA.PHYTO()

        self.assertEqual(phyto.chl, 0.2)
        self.assertEqual(phyto.profiltype, phyto.Homogeneous)
        self.assertEqual(phyto.jd.mrwa, 1.05)
        self.assertEqual(phyto.jd.miwa, -0.0)
        self.assertEqual(phyto.jd.slope, 4.0)
        self.assertEqual(phyto.jd.rmin, 0.01)
        self.assertEqual(phyto.jd.rmax, 200)
        self.assertEqual(phyto.jd.rate, 1.0)
        self.assertIsNone(phyto.sm)
        self.assertIsNone(phyto.tm)

        # Test the method to set the primary mode
        phyto.SetPrimaryMode(1.06, 0.01, 4.2, 0.02, 210, 0.8)
        self.assertEqual(phyto.jd.mrwa, 1.06)
        self.assertEqual(phyto.jd.miwa, 0.01)
        self.assertEqual(phyto.jd.slope, 4.2)
        self.assertEqual(phyto.jd.rmin, 0.02)
        self.assertEqual(phyto.jd.rmax, 210)
        self.assertEqual(phyto.jd.rate, 0.8)
        self.assertIsNone(phyto.sm)
        self.assertIsNone(phyto.tm)

        # Test the method to set the secundary mode
        phyto.SetSecondaryMode(1.05, -0.0, 1, 0.2, 0.1)
        self.assertEqual(phyto.sm.mrwa, 1.05)
        self.assertEqual(phyto.sm.miwa, -0.0)
        self.assertEqual(phyto.sm.sdradius, 1)
        self.assertEqual(phyto.sm.sdvar, 0.2)
        self.assertEqual(phyto.sm.rate, 0.1)
        self.assertIsNone(phyto.tm)

        # Test the method to set the tertiary mode
        phyto.SetTertiaryMode(1.05, -0.0, 1, 0.2, 0.1)
        self.assertEqual(phyto.tm.mrwa, 1.05)
        self.assertEqual(phyto.tm.miwa, -0.0)
        self.assertEqual(phyto.tm.sdradius, 1)
        self.assertEqual(phyto.tm.sdvar, 0.2)
        self.assertEqual(phyto.tm.rate, 0.1)

        # Check profile type

        phyto.SetProfilType(phyto.Homogeneous)
        self.assertEqual(phyto.profiltype, phyto.Homogeneous)

        phyto.SetProfilType(phyto.Gaussian, chlbg=0.1, deep=1, width=2)
        self.assertEqual(phyto.profiltype, phyto.Gaussian)
        self.assertEqual(phyto.gp.chlbg, 0.1)
        self.assertEqual(phyto.gp.deep, 1)
        self.assertEqual(phyto.gp.width, 2)

        phyto.SetProfilType(phyto.UserDefined, userfile="test/emptyfile.txt")
        self.assertEqual(phyto.usefile, "test/emptyfile.txt")

        phyto.SetProfilType(phyto.UserDefined, userfile="fail")
        self.assertRaises(FileNotFoundError)

        # We check the fail check
        with self.assertRaises(Exception) as context:
            phyto.SetProfilType("Bad")
        self.assertTrue("Invalid profile type." in str(context.exception))

    def testSED(self):
        sed = pyOSOAA.SED()

        self.assertEqual(sed.csed, 0.0)
        self.assertIsNone(sed.jd)
        self.assertIsNone(sed.sm)
        self.assertIsNone(sed.tm)

        # Test the method to set the primary mode
        sed.SetPrimaryMode(1.06, 0.01, 4.2, 0.02, 210, 0.8)
        self.assertEqual(sed.jd.mrwa, 1.06)
        self.assertEqual(sed.jd.miwa, 0.01)
        self.assertEqual(sed.jd.slope, 4.2)
        self.assertEqual(sed.jd.rmin, 0.02)
        self.assertEqual(sed.jd.rmax, 210)
        self.assertEqual(sed.jd.rate, 0.8)
        self.assertIsNone(sed.sm)
        self.assertIsNone(sed.tm)

        # Test the method to set the secundary mode
        sed.SetSecondaryMode(1.05, -0.0, 1, 0.2, 0.1)
        self.assertEqual(sed.sm.mrwa, 1.05)
        self.assertEqual(sed.sm.miwa, -0.0)
        self.assertEqual(sed.sm.sdradius, 1)
        self.assertEqual(sed.sm.sdvar, 0.2)
        self.assertEqual(sed.sm.rate, 0.1)
        self.assertIsNone(sed.tm)

        # Test the method to set the tertiary mode
        sed.SetTertiaryMode(1.05, -0.0, 1, 0.2, 0.1)
        self.assertEqual(sed.tm.mrwa, 1.05)
        self.assertEqual(sed.tm.miwa, -0.0)
        self.assertEqual(sed.tm.sdradius, 1)
        self.assertEqual(sed.tm.sdvar, 0.2)
        self.assertEqual(sed.tm.rate, 0.1)

    def testYS(self):
        self.assertEqual(pyOSOAA.YS().abs440, 0.0)
        self.assertIsNone(pyOSOAA.YS().swa)
        self.assertEqual(pyOSOAA.YS(1.0, 1.0).abs440, 1.0)
        self.assertEqual(pyOSOAA.YS(1.0, 2.0).swa, 2.0)
        # We check the fail check
        with self.assertRaises(Exception) as context:
            pyOSOAA.YS(-1)
        self.assertTrue("Invalid absorption value." in str(context.exception))

    def testDET(self):
        self.assertEqual(pyOSOAA.DET().abs440, 0.0)
        self.assertIsNone(pyOSOAA.DET().swa)
        self.assertEqual(pyOSOAA.DET(1.0, 1.0).abs440, 1.0)
        self.assertEqual(pyOSOAA.DET(1.0, 2.0).swa, 2.0)
        # We check the fail check
        with self.assertRaises(Exception) as context:
            pyOSOAA.DET(-1)
        self.assertTrue("Invalid absorption value." in str(context.exception))

    def testAP(self):
        ap = pyOSOAA.AP()

        self.assertEqual(ap.pressure, 1013.00)
        self.assertEqual(ap.hr, 8.0)
        self.assertEqual(ap.ha, 2.0)
        self.assertIsNone(ap.mot)

        ap.SetMot(mot=0.1, hr=6.0)
        self.assertEqual(ap.mot, 0.1)
        self.assertEqual(ap.hr, 6.0)
        self.assertIsNone(ap.pressure)

        ap.SetPressure(1013.25)
        self.assertEqual(ap.pressure, 1013.25)
        self.assertIsNone(ap.mot)

    def testAER(self):
        aer = pyOSOAA.AER()
        self.assertEqual(aer.waref, 0.55)
        self.assertEqual(aer.aotref, 0.1)
        self.assertIsNone(aer.tronca)
        self.assertEqual(aer.model, 2)
        self.assertEqual(aer.sf.model, 3)
        self.assertEqual(aer.sf.rh, 98)

        aer.SetModel(model=0, sdtype=1)
        self.assertIsNone(aer.wmo)
        self.assertIsNone(aer.sf)
        self.assertIsNone(aer.lnd)
        self.assertIsNone(aer.external)
        self.assertEqual(aer.model, 0)
        self.assertEqual(aer.mm.sdtype, 1)
        self.assertListEqual(list(vars(aer.mm).keys()), ["sdtype",
                                                         "mrwa", "miwa",
                                                         "mrwaref", "miwaref",
                                                         "sdradius", "sdvar"])

        aer.SetModel(model=0, sdtype=2)
        self.assertIsNone(aer.wmo)
        self.assertIsNone(aer.sf)
        self.assertIsNone(aer.lnd)
        self.assertIsNone(aer.external)
        self.assertEqual(aer.model, 0)
        self.assertEqual(aer.mm.sdtype, 2)
        self.assertListEqual(list(vars(aer.mm).keys()), ["sdtype", "mrwa",
                                                         "miwa", "mrwaref",
                                                         "miwaref", "slope",
                                                         "rmin", "rmax"])

        aer.SetModel(model=1, wmotype=2)
        self.assertIsNone(aer.mm)
        self.assertIsNone(aer.sf)
        self.assertIsNone(aer.lnd)
        self.assertIsNone(aer.external)
        self.assertEqual(aer.model, 1)
        self.assertEqual(aer.wmo.model, 2)

        aer.SetModel(model=1, wmotype=4, dl=0.4, ws=0.3, oc=0.2, so=0.1)
        self.assertIsNone(aer.mm)
        self.assertIsNone(aer.sf)
        self.assertIsNone(aer.lnd)
        self.assertIsNone(aer.external)
        self.assertEqual(aer.model, 1)
        self.assertEqual(aer.wmo.model, 4)
        self.assertEqual(aer.wmo.dl, 0.4)
        self.assertEqual(aer.wmo.ws, 0.3)
        self.assertEqual(aer.wmo.oc, 0.2)
        self.assertEqual(aer.wmo.so, 0.1)

        aer.SetModel(model=2, sfmodel=2, rh=50)
        self.assertEqual(aer.model, 2)
        self.assertEqual(aer.sf.model, 2)
        self.assertEqual(aer.sf.rh, 50)
        self.assertIsNone(aer.mm)
        self.assertIsNone(aer.wmo)
        self.assertIsNone(aer.lnd)
        self.assertIsNone(aer.external)

        aer.SetModel(model=3, vcdef=1)
        self.assertEqual(aer.model, 3)
        self.assertEqual(aer.lnb.vcdef, 1)
        self.assertIsNone(aer.mm)
        self.assertIsNone(aer.wmo)
        self.assertIsNone(aer.sf)
        self.assertIsNone(aer.external)
        self.assertListEqual(list(vars(aer.lnb).keys()), ["vcdef",
                                                          "coarsevc", "finevc",
                                                          "cmrwa", "cmiwa",
                                                          "cmrwaref",
                                                          "cmiwaref",
                                                          "csdradius",
                                                          "csdvar",
                                                          "fmrwa", "fmiwa",
                                                          "fmrwaref",
                                                          "fmiwaref",
                                                          "fsdradius",
                                                          "fsdvar"])

        aer.SetModel(model=3, vcdef=2)
        self.assertEqual(aer.model, 3)
        self.assertEqual(aer.lnb.vcdef, 2)
        self.assertIsNone(aer.mm)
        self.assertIsNone(aer.wmo)
        self.assertIsNone(aer.sf)
        self.assertIsNone(aer.external)
        self.assertListEqual(list(vars(aer.lnb).keys()), ["vcdef",
                                                          "raot",
                                                          "cmrwa", "cmiwa",
                                                          "cmrwaref",
                                                          "cmiwaref",
                                                          "csdradius",
                                                          "csdvar",
                                                          "fmrwa", "fmiwa",
                                                          "fmrwaref",
                                                          "fmiwaref",
                                                          "fsdradius",
                                                          "fsdvar"])

        aer.SetModel(model=4, extdata="test.txt")
        self.assertEqual(aer.model, 4)
        self.assertEqual(aer.extdata, "test.txt")
        self.assertIsNone(aer.mm)
        self.assertIsNone(aer.wmo)
        self.assertIsNone(aer.sf)
        self.assertIsNone(aer.lnd)

    def tetHYD(self):
        hyd = pyOSOAA.HYD()
        self.assertEqual(hyd.model, 1)
        self.assertIsNone(hyd.extdata)

    def testANG(self):
        ang = pyOSOAA.ANG()

        self.assertIsNone(ang.rad.nbgauss)
        self.assertIsNone(ang.rad.userangfile)

        self.assertIsNone(ang.mie.nbgauss)
        self.assertIsNone(ang.mie.userangfile)

        self.assertEqual(ang.thetas, 30)

    def testSOS(self):
        self.assertIsNone(pyOSOAA.SOS().igmax)

    def testVIEW(self):
        view = pyOSOAA.VIEW()

        self.assertEqual(view.phi, 0)
        self.assertEqual(view.level, 5)
        self.assertEqual(view.z, -10)
        self.assertEqual(view.vza, 0)

    def testOSOAA(self):
        s = pyOSOAA.OSOAA(resroot="test")
        self.assertEqual(s.resroot, "test")
        self.assertEqual(s.wa, 0.440)
        s = pyOSOAA.OSOAA()
        self.assertIsNotNone(s.resroot)
        self.assertIsNotNone(s.root)
        self.assertListEqual(list(vars(s).keys()), ["wa", "root",
                                                    "resroot", "sea",
                                                    "log", "results",
                                                    "dirmie", "phyto",
                                                    "sed", "ys",
                                                    "det", "ap",
                                                    "aer", "hyd",
                                                    "ang", "sos",
                                                    "view"])

        s.run()
        self.assertTrue(os.path.exists(s.resroot))
        self.assertTrue(os.path.exists(s.dirmie.aer))
        self.assertTrue(os.path.exists(s.dirmie.hyd))
        self.assertTrue(os.path.exists(s.dirmie.sea))
        self.assertTrue(os.path.isfile(s.resroot+"/script.kzh"))
        self.assertTrue(os.path.isfile(s.resroot
                                       + "/Standard_outputs/"
                                       + "LUM_vsVZA.txt"))
        self.assertTrue(s.outputs.vsvza.I[51], 0.128266)


class TestOSOAAHelpers(unittest.TestCase):

    def testConfigureOceanBlack(self):
        s = pyOSOAA.OSOAA()
        wl = [0.44, 0.55, 0.66]
        view = 0
        expected = [0.0, 0.0, 0.0]
        s.view.level = 4
        s = pyOSOAA.osoaahelpers.ConfigureOcean(s, "black")
        result = pyOSOAA.osoaahelpers.RunWavelengths(s, wl, view)
        self.assertListEqual(list(result), expected)

    def testConfigureOceanBlackErrorCatch(self):
        s = pyOSOAA.OSOAA()
        with self.assertRaises(Exception) as context:
            pyOSOAA.osoaahelpers.ConfigureOcean(s, "wrong")
        self.assertTrue("Wrong ocean type." in str(context.exception))

    def testRunWavelenghtsSameAngle(self):
        s = pyOSOAA.OSOAA()
        wl = [0.44, 0.55, 0.66]
        view = 0
        expected = [0.128266, 0.639862E-01, 0.399832E-04]
        result = pyOSOAA.osoaahelpers.RunWavelengths(s, wl, view)
        self.assertListEqual(list(result), expected)

    def testRunWavelenghtsDifferentAngle(self):
        s = pyOSOAA.OSOAA()
        wl = [0.44, 0.55, 0.66]
        view = [0, 1.43, -1.43]
        expected = [0.128266, 0.639702E-01, 0.399984E-04]
        result = pyOSOAA.osoaahelpers.RunWavelengths(s, wl, view)
        self.assertListEqual(list(result), expected)

    def testRunWavelenghtsOtherOutput(self):
        s = pyOSOAA.OSOAA()
        wl = [0.44, 0.55, 0.66]
        view = 0
        expected = [0.229523, 0.158827, 0.456162E-02]
        result = pyOSOAA.osoaahelpers.RunWavelengths(s, wl, view, "refl")
        self.assertListEqual(list(result), expected)

    def testRunWavelenghtsErrorCatch(self):
        s = pyOSOAA.OSOAA()
        with self.assertRaises(Exception) as context:
            pyOSOAA.osoaahelpers.RunWavelengths(s, 0.5, 0, "wrong")
        self.assertTrue("Wrong output variable." in str(context.exception))



if __name__ == '__main__':
    unittest.main()
