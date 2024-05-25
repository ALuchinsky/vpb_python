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

def computePL(D, homDim, scaleSeq, k=1):
    birth, death = D[homDim][:,0], D[homDim][:,1]
    Lambda = [
        np.sort(pmax(0, np.apply_along_axis(min, 0, np.array([s-birth, death-s]))))[-k]
        for s in scaleSeq]
    return np.array(Lambda)

def computePS(D, homDim, scaleSeq, p=1):
    x, y = D[homDim][:,0], D[homDim][:,1]
    pp = (y-x)**p
    w = pp/np.sum(pp)

    phi = []
    for k in range(len(scaleSeq)-1):
        alpha1 = pmax(scaleSeq[k], x)
        alpha2 = pmax(scaleSeq[k], (x+y)/2)
        beta1 = pmin(scaleSeq[k+1], (x+y)/2)
        beta2 = pmin(scaleSeq[k+1], y)
        b1 = pmax(0,beta1-alpha1)*((beta1+alpha1)/2-x)
        b2 = pmax(0,beta2-alpha2)*(y-(beta2+alpha2)/2)
        phi.append( np.sum(w*(b1+b2))/(scaleSeq[k+1]-scaleSeq[k]))
    return np.array(phi)

def computeNL(D, homDim, scaleSeq):
    x, y = D[homDim][:,0], D[homDim][:,1]
    lL = (y-x)/sum(y-x)
    nl = []
    for k in range(len(scaleSeq)-1):
        b = pmin(scaleSeq[k+1],y)-pmax(scaleSeq[k],x)
        nl.append( np.sum(lL*pmax(0,b))/(scaleSeq[k+1]-scaleSeq[k]))
    return np.array(nl)

def computeVAB(D, homDim, scaleSeq):
    x, y = D[homDim][:,0], D[homDim][:,1]
    vab = []
    for k in range( len(scaleSeq)-1):
        b = pmin(scaleSeq[k+1],y)-pmax(scaleSeq[k],x)
        vab.append( sum(pmax(0,b))/(scaleSeq[k+1]-scaleSeq[k]))
    return np.array(vab)
