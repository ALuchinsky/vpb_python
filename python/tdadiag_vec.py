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

    def getData(self):
        return self.data
    
    def getThreshold(self):
        return self.threshold
    
    def setThreshold(self, threshold_):
        self.threshold = threshold_

    def setData(self, data_):
        self.data = data_
        self.diag = None
        self.pd = None
        self.max_scale = None

    def setScale(self, scale_ = None, nGrid = 11):
        if scale_ is None:
            scale_ = np.linspace(0, self.max_scale, nGrid)
        self.scale = scale_
        return self.scale

    def getScale(self):
        return self.scale

    def calcDiag(self, threshold_ = None, inf = None):
        if threshold_ is None:
            threshold_ = self.threshold
        self.diag = ripser.ripser(self.data, thresh=threshold_)["dgms"]
        if inf is not None:
            self.diag[0][-1, 1] = inf
        else:
            self.diag[0] = self.diag[0][:-1,:]
        self.max_scale = np.max(self.diag[0][:,1])


        return self.diag
    
    def getDiag(self):
        if self.diag is None:
            self.calcDiag()
        return self.diag

    def __init__(self, data_ = None, threshold_ = 2) -> None:
        self.setData(data_)
        self.setThreshold(threshold_)

    def computePS(self, homDim = 1, scaleSeq_ = None, nGrid = 11, p=1):
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
        scaleSeq = self.setScale(scale_ = scaleSeq_, nGrid = nGrid)
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
    
    def computeNL(self, homDim = 1, scaleSeq_ = None, nGrid = 11, p=1):
        """
        Compute theNormalized Life Curve vectorization for a given homological dimension, scale sequence, and power.

        Parameters:
            D (numpy.ndarray): Persistence diagram (array of birth-death arrays for each dimension).
            homDim (int): The homological dimension along which the NL is computed.
            scaleSeq (numpy.ndarray): The sequence of scale values.

        Returns:
            numpy.ndarray: The nonlinear landscape vector.
        """
        scaleSeq = self.setScale(scale_ = scaleSeq_, nGrid = nGrid)
        D = self.getDiag()
        x, y = D[homDim][:,0], D[homDim][:,1]
        lL = (y-x)/sum(y-x)
        nl = []
        for k in range(len(scaleSeq)-1):
            b = pmin(scaleSeq[k+1],y)-pmax(scaleSeq[k],x)
            nl.append( np.sum(lL*pmax(0,b))/(scaleSeq[k+1]-scaleSeq[k]))
        return np.array(nl)

    def computeVAB(self, homDim = 1, scaleSeq_ = None, nGrid = 11):
        """
        Compute the Vector Summary of the Betti Curve    (VAB) vectorization for a given homological dimension, scale sequence, and power.

        Parameters:
            D (numpy.ndarray): Persistence diagram (array of birth-death arrays for each dimension).
            homDim (int): The homological dimension along which the VAB is computed.
            scaleSeq (numpy.ndarray): The sequence of scale values.

        Returns:
            numpy.ndarray: The VAB vector.
        """
        scaleSeq = self.setScale(scale_ = scaleSeq_, nGrid = nGrid)
        D = self.getDiag()
        x, y = D[homDim][:,0], D[homDim][:,1]
        vab = []
        for k in range( len(scaleSeq)-1):
            b = pmin(scaleSeq[k+1],y)-pmax(scaleSeq[k],x)
            vab.append( sum(pmax(0,b))/(scaleSeq[k+1]-scaleSeq[k]))
        return np.array(vab)
    
    def computeVAB(self, homDim = 1, scaleSeq_ = None, nGrid = 11):
        """
        Compute the Vector Summary of the Betti Curve    (VAB) vectorization for a given homological dimension, scale sequence, and power.

        Parameters:
            D (numpy.ndarray): Persistence diagram (array of birth-death arrays for each dimension).
            homDim (int): The homological dimension along which the VAB is computed.
            scaleSeq (numpy.ndarray): The sequence of scale values.

        Returns:
            numpy.ndarray: The VAB vector.
        """
        scaleSeq = self.setScale(scale_ = scaleSeq_, nGrid = nGrid)
        D = self.getDiag()
        x, y = D[homDim][:,0], D[homDim][:,1]
        vab = []
        for k in range( len(scaleSeq)-1):
            b = pmin(scaleSeq[k+1],y)-pmax(scaleSeq[k],x)
            vab.append( sum(pmax(0,b))/(scaleSeq[k+1]-scaleSeq[k]))
        return np.array(vab)


    def computeECC(self, homDim = 1, scaleSeq_ = None, nGrid = 11):
        """
        Compute the Euler Characteristic Curve (ECC) vectorization for a given homological dimension, maximum homological dimension, and scale sequence.

        Parameters:
            D (numpy.ndarray): Persistence diagram (array of birth-death arrays for each dimension).
            maxhomDim (int): The maximum homological dimension.
            scaleSeq (numpy.ndarray): The sequence of scale values.

        Returns:
            numpy.ndarray: The ECC vector.
        """
        scaleSeq = self.setScale(scale_ = scaleSeq_, nGrid = nGrid)
        D = self.getDiag()
        ecc = np.zeros( len(scaleSeq)-1)
        for d in range(homDim+1):
            ecc = ecc + (-1)**d * self.computeVAB(d, scaleSeq, nGrid)
        return ecc
