import ripser
import numpy as np

def pmax(num, vec):
    """
    Compute the element-wise maximum of a scalar value and a NumPy array.

    Parameters:
        num (float): The scalar value.
        vec (numpy.ndarray): The input array.

    Returns:
        numpy.ndarray: The resulting array with the element-wise maximum.
    """
    return np.array([max(num, vec[i_]) for i_ in range(vec.size)])

def pmin(num, vec):
    """
    Compute the element-wise minimum of a scalar value and a NumPy array.

    Parameters:
        num (float): The scalar value.
        vec (numpy.ndarray): The input array.

    Returns:
        numpy.ndarray: The resulting array with the element-wise minimum.
    """
    return np.array([min(num, vec[i_]) for i_ in range(vec.size)])


class tdadiag_vect:
    data = None
    diag = None
    pd = None
    threshold = 2
    max_scale = None
    scale = None
    inf = 2

    def getData(self):
        return self.data
    
    def getThreshold(self):
        return self.threshold
    
    def setThreshold(self, threshold_):
        self.threshold = threshold_

    def setInf(self, inf):
        self.inf= inf
        return

    def setData(self, data_):
        self.data = data_
        self.diag = None
        self.max_scale = None

    def setScale(self, scale = None, nGrid = 11, quantiles = False):
        if scale is None:
            if quantiles:
                pd = self.getDiag()
                scale = np.quantile(pd[1][:,1], np.linspace(0, 1, nGrid))
            else:
                scale = np.linspace(0, self.max_scale, nGrid)
        self.scale = scale
        return self.scale

    def getScale(self):
        return self.scale

    def calcDiag(self, threshold = None, inf = None):
        if threshold is None:
            threshold = self.threshold
        if inf is None:
            inf = self.inf
        self.diag = ripser.ripser(self.data, thresh=threshold)["dgms"]
        if inf is not None:
            self.diag[0][-1, 1] = inf
        else:
            self.diag[0] = self.diag[0][:-1,:]
        self.max_scale = np.max(self.diag[0][:,1])
        return self.diag
    
    def getDiag(self, threshold = None, inf = None):
        if threshold is None:
            threshold = self.threshold
        if inf is None:
            inf = self.inf
        if self.diag is None:
            self.calcDiag(threshold=threshold, inf=inf)
        return self.diag

    def getPD(self, threshold = None, inf = None):
        diag = self.getDiag(threshold=threshold, inf=inf)
        return [
            diag[0], 
            np.vstack([diag[1][:,0], diag[1][:,1]-diag[1][:,0]]).transpose()
            ]
    
    def __init__(self, data = None, threshold = 2, inf = 2) -> None:
        self.setData(data)
        self.setThreshold(threshold)
        self.setInf(inf)

    def computePS(self, homDim = 1, scaleSeq = None, nGrid = 11, p=1, quantiles = False):
        """
        Compute the Persistence Silhouette vectorization for a given homological dimension, scale sequence, and power.

        Parameters:
            D (numpy.ndarray): Persistence diagram (array of birth-death arrays for each dimension).
            homDim (int): The homological dimension along which the PS is computed.
            scaleSeq (numpy.ndarray): The sequence of scale values.
            p (int, optional): The power to raise the difference between y and x. Defaults to 1.

        Returns:
            numpy.ndarray: The persistence spectrum vector.
        """
        scaleSeq = self.setScale(scale = scaleSeq, nGrid = nGrid, quantiles = quantiles)
        D = self.getDiag()
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
    
    def computeNL(self, homDim = 1, scaleSeq = None, nGrid = 11, p=1, quantiles = False):
        """
        Compute theNormalized Life Curve vectorization for a given homological dimension, scale sequence, and power.

        Parameters:
            D (numpy.ndarray): Persistence diagram (array of birth-death arrays for each dimension).
            homDim (int): The homological dimension along which the NL is computed.
            scaleSeq (numpy.ndarray): The sequence of scale values.

        Returns:
            numpy.ndarray: The nonlinear landscape vector.
        """
        scaleSeq = self.setScale(scale = scaleSeq, nGrid = nGrid, quantiles = quantiles)
        D = self.getDiag()
        x, y = D[homDim][:,0], D[homDim][:,1]
        lL = (y-x)/sum(y-x)
        nl = []
        for k in range(len(scaleSeq)-1):
            b = pmin(scaleSeq[k+1],y)-pmax(scaleSeq[k],x)
            nl.append( np.sum(lL*pmax(0,b))/(scaleSeq[k+1]-scaleSeq[k]))
        return np.array(nl)

    def computeVAB(self, homDim = 1, scaleSeq = None, nGrid = 11, quantiles = False):
        """
        Compute the Vector Summary of the Betti Curve    (VAB) vectorization for a given homological dimension, scale sequence, and power.

        Parameters:
            D (numpy.ndarray): Persistence diagram (array of birth-death arrays for each dimension).
            homDim (int): The homological dimension along which the VAB is computed.
            scaleSeq (numpy.ndarray): The sequence of scale values.

        Returns:
            numpy.ndarray: The VAB vector.
        """
        scaleSeq = self.setScale(scale = scaleSeq, nGrid = nGrid, quantiles = quantiles)
        D = self.getDiag()
        x, y = D[homDim][:,0], D[homDim][:,1]
        vab = []
        for k in range( len(scaleSeq)-1):
            b = pmin(scaleSeq[k+1],y)-pmax(scaleSeq[k],x)
            vab.append( sum(pmax(0,b))/(scaleSeq[k+1]-scaleSeq[k]))
        return np.array(vab)
    
    def computeECC(self, homDim = 1, scaleSeq = None, nGrid = 11, quantiles = False):
        """
        Compute the Euler Characteristic Curve (ECC) vectorization for a given homological dimension, maximum homological dimension, and scale sequence.

        Parameters:
            D (numpy.ndarray): Persistence diagram (array of birth-death arrays for each dimension).
            maxhomDim (int): The maximum homological dimension.
            scaleSeq (numpy.ndarray): The sequence of scale values.

        Returns:
            numpy.ndarray: The ECC vector.
        """
        scaleSeq = self.setScale(scale = scaleSeq, nGrid = nGrid, quantiles = quantiles)
        D = self.getDiag()
        ecc = np.zeros( len(scaleSeq)-1)
        for d in range(homDim+1):
            ecc = ecc + (-1)**d * self.computeVAB(d, scaleSeq, nGrid)
        return ecc
    
    def computeVPB_dim0(self, x, y, ySeq, lam):
        """
        Compute the VPB values for dimension 0.

        Parameters:
            x (numpy.ndarray): The x values.
            y (numpy.ndarray): The y values.
            ySeq (numpy.ndarray): The sequence of y values.
            lam (numpy.ndarray): The lambda values.

        Returns:
            numpy.ndarray: The computed VPB values.
        """
        n = ySeq.shape[0] - 1
        nPoints = y.shape[0]
        vpb = np.zeros(n, dtype=np.float64)
        for i in range(n):
            c = ySeq[i]
            d = ySeq[i+1]
            for j in range(nPoints):
                if c - lam[j] < y[j] < d + lam[j]:
                    y_cd = y[j]
                    lam_cd = lam[j]
                    yMin = max(c, y_cd - lam_cd)
                    yMax = min(d, y_cd + lam_cd)
                    vpb[i] += 0.5 * (yMax**2 - yMin**2) / (ySeq[i+1] - ySeq[i])
        return vpb


    def computeVPB_dim1(self, x, y, xSeq, ySeq, lam):
        """
        Compute the Vector Persistence Block (VPB) vectorization for a given set of points in dimension 1.

        Parameters:
            x (numpy.ndarray): The x-coordinates of the points.
            y (numpy.ndarray): The y-coordinates of the points.
            xSeq (numpy.ndarray): The x-coordinates of the grid points.
            ySeq (numpy.ndarray): The y-coordinates of the grid points.
            lam (numpy.ndarray): The lambda values.

        Returns:
            numpy.ndarray: The VPB matrix.
        """
        n = xSeq.shape[0] - 1
        m = ySeq.shape[0] - 1
        dx = np.diff(xSeq)
        dy = np.diff(ySeq)
        vpb = np.zeros((n, m))

        for i in range(n):
            a, b = xSeq[i], xSeq[i+1]
            for j in range(m):
                c, d = ySeq[j], ySeq[j+1]
                # Using bitwise operations for logical conditions
                xCond = (x + lam >= a) & (x - lam <= b)
                yCond = (y + lam >= c) & (y - lam <= d)
                inds = np.where(xCond & yCond)[0]
                num_inds = inds.shape[0]

                if num_inds > 0:
                    xInd = x.take(inds)
                    yInd = y.take(inds)
                    lamInd = lam.take(inds)
                    xMin = np.maximum(a, xInd - lamInd)
                    xMax = np.minimum(b, xInd + lamInd)
                    yMin = np.maximum(c, yInd - lamInd)
                    yMax = np.minimum(d, yInd + lamInd)
                    add = 0.5 * np.sum((xMax - xMin) * (yMax - yMin) * (xMax + xMin + yMax + yMin)) / dx[i] / dy[j]
                    vpb[i, j] += add
        return np.asarray(vpb)

    def computeVPB(self, homDim, xSeq=None, ySeq=None, tau=0.3, flatten = False):
        """
        Compute the VPB vectorization using the given parameters.

        Parameters:
            PD (list): Persistence Diagram (list of birth-persistence arrays for each dimension).
            homDim (int): The dimension along which the homogeneity is computed.
            xSeq (numpy.ndarray): The x-coordinates of the grid points.
            ySeq (numpy.ndarray): The y-coordinates of the grid points.
            tau (float, optional): The tau value. Defaults to 0.3.

        Returns:
            numpy.ndarray: The VPB matrix.
        """
        PD = self.getPD()
        x = PD[homDim][:,0]
        y = PD[homDim][:,1]
        lam = tau * y
        xSeq_ = self.getScale() if xSeq is None else xSeq
        ySeq_ = self.getScale() if ySeq is None else ySeq
        if homDim == 0:
            return self.computeVPB_dim0(x, y, ySeq_, lam)
        else:
            out = self.computeVPB_dim1(x, y, xSeq_, ySeq_, lam)
            if flatten:
                out = np.transpose(out).reshape( out.shape[0]*out.shape[1])
            return out


