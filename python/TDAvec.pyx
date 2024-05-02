import numpy as np

def DiagToPD(D):
    PD = [ np.array([D[dim][:,0], D[dim][:,1] - D[dim][:,0]]) for dim in range(len(D))]
    return PD

def computeVPB_dim0(x, y, ySeq, lam):
    dy = np.diff(ySeq)
    vpb = np.zeros( len(dy))
    for i in range(len(dy)):
        c = ySeq[i]
        d = ySeq[i+1]
        for j in range( len(y)):
            if c - lam[j] < y[j] and y[j] < d + lam[j]:
                y_cd = y[j]
                lam_cd = lam[j]
                yMin = max(c, y_cd - lam_cd)
                yMax = min(d, y_cd + lam_cd)
                vpb[i] += 0.5*(yMax**2 - yMin**2)/dy[i]
    return vpb

def pmax(num, vec):
    return np.array([max(num, vec[i_]) for i_ in range(vec.size)])
def pmin(num, vec):
    return np.array([min(num, vec[i_]) for i_ in range(vec.size)])

def computeVPB_dim1(x, y, xSeq, ySeq, lam):
    dx = np.diff(xSeq)
    dy = np.diff(ySeq)
    vpb = np.zeros( (dx.size, dy.size) )
    for i in range(dx.size):
        a, b = xSeq[i], xSeq[i+1]
        for j in range(dy.size):
            c, d = ySeq[j], ySeq[j+1]
            xCond = (x+lam >= a) & (x-lam <= b)
            yCond = (y+lam >= c) & (y-lam <= d)
            inds = np.where(xCond & yCond)[0]
            if len(inds)>0:
                xInd, yInd, lamInd = x.take(inds), y.take(inds), lam.take(inds)
                xMin, xMax = pmax(a, xInd - lamInd), pmin(b, xInd + lamInd)
                yMin, yMax = pmax(c, yInd - lamInd), pmin(d, yInd + lamInd)
                add = 0.5*np.sum( (xMax-xMin)*(yMax-yMin)*(xMax+xMin+yMax+yMin))/dx[i]/dy[j]
                vpb[i, j] += add
    return vpb

def computeVPB(PD, homDim, xSeq, ySeq, tau=0.3):
    x = PD[homDim][0]
    y = PD[homDim][1]
    lam = tau * y
    if homDim == 0:
        return computeVPB_dim0(x, y, ySeq, lam)
    else:
        return computeVPB_dim1(x, y, xSeq, ySeq, lam)


