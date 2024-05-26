import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import ripser
import persim

import TDAvec
from TDAvec import pmin, pmax, DiagToPD, \
    computeVPB, computePL, computePS, computeNL, computeVAB, computeECC, computePES, computePI

########### Download point cloud and create diagrams
X = np.loadtxt("../R/unitCircle.csv", skiprows=1, delimiter=",")
D = ripser.ripser(X, thresh=2)["dgms"]
D[0][-1, 1] = 2
PD = DiagToPD(D)

######### computeVPB
print("computeVPB:")
ySeqH0 = np.quantile(PD[0][:,1], np.arange(0, 1.1, 0.2))
xSeqH1 = np.quantile(PD[1][:,0], np.arange(0, 1.1, 0.2))
ySeqH1 = np.quantile(PD[1][:,1], np.arange(0, 1.1, 0.2))

vpb0 = computeVPB(PD, homDim=0, xSeq=[], ySeq=ySeqH0)
vpb0_R = np.loadtxt("../R/vpb_0.csv", skiprows=1)
cond = np.allclose(vpb0, vpb0_R) 
print(f"\t dim={0}: {cond}")
assert cond
vpb1 = computeVPB(PD, homDim = 1, xSeq=xSeqH1, ySeq=ySeqH1)
vpb1 = np.transpose(vpb1).reshape( (25,))
vpb1_R = np.loadtxt("../R/vpb_1.csv", skiprows=1)
cond = np.allclose(vpb1, vpb1_R, atol=1e-7)
print(f"\t dim={1}: {cond}")
assert cond

##### computePI
print("computePI:")
resB, resP = 5, 5
minPH0, maxPH0 = np.min(PD[0][:,1]), np.max(PD[0][:,1])
ySeqH0 = np.linspace(minPH0, maxPH0, resP+1)
xSeqH0 = np.zeros( resB+1)
minBH1, maxBH1 = np.min(PD[1][:,0]), np.max(PD[1][:,0])
xSeqH1 = np.linspace(minBH1, maxBH1, resB+1)
minPH1, maxPH1 = np.min(PD[1][:,1]), np.max(PD[1][:,1])
ySeqH1 = np.linspace(minPH1, maxPH1, resP+1)
sigma = 0.5*(maxPH0-minPH0)/resP
pi0 = computePI(PD, homDim = 0, xSeq = xSeqH0, ySeq = ySeqH0, sigma = sigma)
pi0_R = np.loadtxt("../R/pi_0.csv", skiprows=1)
cond = np.allclose(pi0, pi0_R)
print(f"\t dim={0}: {cond}")
assert cond
sigma = 0.5*(maxPH1-minPH1)/resP
pi1 = computePI(PD, homDim = 1, xSeq = xSeqH1, ySeq = ySeqH1, sigma = sigma)
pi1_R = np.loadtxt("../R/pi_1.csv", skiprows=1)
cond = np.allclose(pi1, pi1_R)
print(f"\t dim={1}: {cond}")
assert cond

## Others

scaleSeq = np.linspace(0, 2, 11)
# universal function for comparison
def compareResults(func, R_prefix, maxDim = 1, D=D, scaleSeq = scaleSeq, atol = 1e-7):
    print(f"Comparing compute{R_prefix.upper()}:")
    # for each of the dims
    for d in range(maxDim+1):
        # calc python result
        pyth = func(D, d, scaleSeq)
        # read R result
        R = np.loadtxt(f"../R/{R_prefix}_{d}.csv", skiprows=1)
        # report if they are equal
        cond = np.allclose(pyth, R, atol=atol)
        print("\t dim=", d, ":", cond)
        assert cond
        
testDict = {"pl":computePL, "ps":computePS, "nl":computeNL, "vab":computeVAB, "ecc":computeECC, "pes": computePES}
for p in testDict.keys():
    compareResults(testDict[p], p)

