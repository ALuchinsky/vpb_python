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

    def calcDiag(self, threshold_ = None, inf = None):
        if threshold_ is None:
            threshold_ = self.threshold
        self.diag = ripser.ripser(self.data, thresh=threshold_)["dgms"]
        if inf is not None:
            self.diag[0][-1, 1] = inf
        else:
            self.diag[0] = self.diag[0][:-1,:]


        return self.diag
    
    def getDiag(self):
        if self.diag is None:
            self.calcDiag()
        return self.diag

    def __init__(self, data_ = None, threshold_ = 2) -> None:
        self.setData(data_)
        self.setThreshold(threshold_)

    def computePS(self, homDim, scaleSeq, p=1):
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

        