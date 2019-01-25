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


if __name__ == '__main__':
    unittest.main()
