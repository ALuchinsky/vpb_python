import numpy as np
import ripser

class TDAvectorizer:

    def __init__(self, params = {"output": "vpb", "threshold": 2, "inf": 2, "maxDim": 1}):
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
        return
    
    def findLimits(self, homDim = 0):
        if self.diags is not None:
            births = np.concatenate([d[homDim][:,0] for d in self.diags])
            deaths = np.concatenate([d[homDim][:,1] for d in self.diags])
            return {
                "minB": births.min(), "maxB": deaths.max(), "minD": deaths.min(), "maxD": births.min()
            }


    

    

