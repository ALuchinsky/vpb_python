import unittest

import TDAvectorizer
from TDAvectorizer import computeNL, computeVAB, computeVPB, computePL, computePS, \
    computeECC, computePES, computePI, computeFDA,\
        DiagToPD, createEllipse, vect
import ripser
import numpy as np

def lists_are_equal(nums1, nums2, atol=1e-6):
    return np.allclose(nums1, nums2, atol=atol)

class Testing_Functions(unittest.TestCase):

    def setUp(self) -> None:
        self.X = np.loadtxt("../R/unitCircle.csv", skiprows=1, delimiter=",")
        self.D = ripser.ripser(self.X, thresh=2)["dgms"]
        self.D[0][-1, 1] = 2
        self.scaleSeq = np.linspace(0, 2, 11)


    def test_PL_0(self):
        python = computePL(self.D, 0, self.scaleSeq)
        R = np.loadtxt("../R/pl_0.csv", skiprows=1)
        self.assertTrue( lists_are_equal(python, R))
    
    def test_PL_1(self):
        python = computePL(self.D, 1, self.scaleSeq)
        R = np.loadtxt("../R/pl_1.csv", skiprows=1)
        self.assertTrue( lists_are_equal(python, R))

    def test_PS_0(self):
        python = computePS(self.D, 0, self.scaleSeq)
        R = np.loadtxt("../R/ps_0.csv", skiprows=1)
        self.assertTrue( lists_are_equal(python, R))

    def test_PS_1(self):
        python = computePS(self.D, 1, self.scaleSeq)
        R = np.loadtxt("../R/ps_1.csv", skiprows=1)
        self.assertTrue( lists_are_equal(python, R))

    def test_NL_0(self):
        python = computeNL(self.D, 0, self.scaleSeq)
        R = np.loadtxt("../R/nl_0.csv", skiprows=1)
        self.assertTrue( lists_are_equal(python, R))

    def test_NL_1(self):
        python = computeNL(self.D, 1, self.scaleSeq)
        R = np.loadtxt("../R/nl_1.csv", skiprows=1)
        self.assertTrue( lists_are_equal(python, R))

    def test_VAB_0(self):
        python = computeVAB(self.D, 0, self.scaleSeq)
        R = np.loadtxt("../R/vab_0.csv", skiprows=1)
        self.assertTrue( lists_are_equal(python, R))

    def test_VAB_1(self):
        python = computeVAB(self.D, 1, self.scaleSeq)
        R = np.loadtxt("../R/vab_1.csv", skiprows=1)
        self.assertTrue( lists_are_equal(python, R))

    def test_ECC_0(self):
        python = computeECC(self.D, 0, self.scaleSeq)
        R = np.loadtxt("../R/ecc_0.csv", skiprows=1)
        self.assertTrue( lists_are_equal(python, R))

    def test_ECC_1(self):
        python = computeECC(self.D, 1, self.scaleSeq)
        R = np.loadtxt("../R/ecc_1.csv", skiprows=1)
        self.assertTrue( lists_are_equal(python, R))

    def test_PES_0(self):
        python = computePES(self.D, 0, self.scaleSeq)
        R = np.loadtxt("../R/pes_0.csv", skiprows=1)
        self.assertTrue( lists_are_equal(python, R))

    def test_PES_1(self):
        python = computePES(self.D, 1, self.scaleSeq)
        R = np.loadtxt("../R/pes_1.csv", skiprows=1)
        self.assertTrue( lists_are_equal(python, R))

    def test_FDA_0(self):
        maxD = max(self.D[0][:,1])
        python = computeFDA(self.D, maxD, 0, 10)
        R = np.loadtxt("../R/fda_0.csv", skiprows=1)
        self.assertTrue( lists_are_equal(python, R))

    def test_FDA_1(self):
        maxD = max(self.D[1][:,1])
        python = computeFDA(self.D, maxD, 1, 10)
        R = np.loadtxt("../R/fda_1.csv", skiprows=1)
        self.assertTrue( lists_are_equal(python, R))

    def test_VPB_0(self):
        ySeqH0 = np.quantile(self.D[0][:,1] - self.D[0][:,0] , np.arange(0, 1.1, 0.2))
        vpb0 = computeVPB(self.D, homDim=0, xSeq = [], ySeq = ySeqH0)
        R = np.loadtxt("../R/vpb_0.csv", skiprows=1)
        self.assertTrue( lists_are_equal(vpb0, R))

    def test_VPB_1(self):
        xSeqH1 = np.quantile(self.D[1][:,0], np.arange(0, 1.1, 0.2))
        ySeqH1 = np.quantile(self.D[1][:,1]- self.D[1][:,0], np.arange(0, 1.1, 0.2))
        vpb1 = computeVPB(self.D, homDim = 1, xSeq=xSeqH1, ySeq=ySeqH1)
        vpb1 = np.transpose(vpb1).reshape( (25,))
        R = np.loadtxt("../R/vpb_1.csv", skiprows=1)
        self.assertTrue( lists_are_equal(vpb1, R))

    def test_PI_0(self):
        resB, resP = 5, 5
        minPH0, maxPH0 = np.min(self.D[0][:,1]), np.max(self.D[0][:,1])
        ySeqH0 = np.linspace(minPH0, maxPH0, resP+1)
        xSeqH0 = np.zeros( resB+1)
        sigma = 0.5*(maxPH0-minPH0)/resP
        pi0 = computePI(self.D, homDim = 0, xSeq = xSeqH0, ySeq = ySeqH0, sigma = sigma)
        pi0_R = np.loadtxt("../R/pi_0.csv", skiprows=1)
        self.assertTrue( lists_are_equal(pi0, pi0_R))

    # def test_PI_1(self):
        resB, resP = 5, 5
        minBH1, maxBH1 = np.min(self.D[1][:,0]), np.max(self.D[1][:,0])
        xSeqH1 = np.linspace(minBH1, maxBH1, resB+1)
        minPH1, maxPH1 = np.min(self.D[1][:,1] - self.D[1][:,0]), np.max(self.D[1][:,1] - self.D[1][:,0])
        ySeqH1 = np.linspace(minPH1, maxPH1, resP+1)
        sigma = 0.5*(maxPH1-minPH1)/resP
        pi1 = computePI(self.D, homDim = 1, xSeq = xSeqH1, ySeq = ySeqH1, sigma = sigma)
        pi1_R = np.loadtxt("../R/pi_1.csv", skiprows=1)
        self.assertTrue( lists_are_equal(pi1, pi1_R))


class Testing_Class(unittest.TestCase):
    def setUp(self):
        np.random.seed(42)
        self.clouds = []
        self.ratList = np.random.uniform(-0.5, 0.5, 10)
        for ratio in self.ratList:
            self.clouds = self.clouds + [createEllipse(a=1-ratio, b=1, eps=0.1)]
        self.vect = TDAvectorizer.TDAvectorizer()
        self.vect.fit(self.clouds)

    def test_N(self):
        self.assertEqual( len(self.vect.diags), 10)



if __name__ == '__main__':
    unittest.main()