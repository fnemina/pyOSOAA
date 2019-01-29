import unittest
import pyOSOAA


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

        phyto.SetProfilType(phyto.UserDefined, userfile="test/profile.txt")
        self.assertEqual(phyto.usefile, "test/profile.txt")

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


if __name__ == '__main__':
    unittest.main()
