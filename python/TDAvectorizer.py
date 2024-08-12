import numpy as np
import ripser

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


class TDAvectorizer:

    def __init__(self, params = {"output": "vpb", "threshold": 2, "inf": 2, "maxDim": 1,
                                 "scale": None, "nGrid": 11, "quantiles": False}):
        self.diags = []
        self.params = params
        return 
    def getParams(self, varName = None):
        if varName is not None:
            return self.params[varName]
        return self.params
    def setParams(self, params):
        for k in params.keys():
            self.params[k] = params[k]
        return
    
    def createEllipse(self, n = 100, a = 1, b = 1, eps = 0.1):
        phi = np.random.uniform(0,2*np.pi,n)
        r = np.random.uniform(1-eps,1+eps,n)
        x = a * r * np.cos(phi)
        y = b * r * np.sin(phi)
        return np.vstack([x,y]).T
    
    def __getitem__(self, index):
        self.diags = self.diags[index]
    
    def fit(self, data):
        self.diags = []
        for d in data:
            diag_ = ripser.ripser(d, thresh=self.params["threshold"])["dgms"]
            if self.params["inf"] is not None:
                diag_[0][-1, 1] = self.params["inf"]
            else:
                diag_[0] = diag_[0][:-1,:]
            self.diags.append(diag_)
        if self.params["scale"] is None:
            limits1 = self.findLimits(1)
            self.params["scale"] = np.linspace(limits1["minB"], limits1["maxB"], self.params["nGrid"])
        return
    
    def findLimits(self, homDim = 0):
        if self.diags is not None:
            births = np.concatenate([d[homDim][:,0] for d in self.diags])
            deaths = np.concatenate([d[homDim][:,1] for d in self.diags])
            return {
                "minB": births.min(), "maxB": deaths.max(), "minD": deaths.min(), "maxD": births.min()
            }

    def computeVAB(self, D, homDim, scaleSeq):
        """
        Compute the Vector Summary of the Betti Curve    (VAB) vectorization for a given homological dimension, scale sequence, and power.

        Parameters:
            D (numpy.ndarray): Persistence diagram (array of birth-death arrays for each dimension).
            homDim (int): The homological dimension along which the VAB is computed.
            scaleSeq (numpy.ndarray): The sequence of scale values.

        Returns:
            numpy.ndarray: The VAB vector.
        """
        x, y = D[homDim][:,0], D[homDim][:,1]
        vab = []
        for k in range( len(scaleSeq)-1):
            b = pmin(scaleSeq[k+1],y)-pmax(scaleSeq[k],x)
            vab.append( sum(pmax(0,b))/(scaleSeq[k+1]-scaleSeq[k]))
        return np.array(vab)
    
    def transform(self, homDim = 1, scaleSeq = None, output = None):
        if output is None:
            output = self.params["output"].lower()
        else:
            output = output.lower()
        if scaleSeq is None:
            scaleSeq = self.params["scale"]
        if output == "vab":
            print("VAB")
            return np.array([self.computeVAB(d, homDim = homDim, scaleSeq = scaleSeq) for d in self.diags])
    


    

    

