import numpy as np 
import matplotlib.pyplot as plt  
from scipy import optimize


class Fitter:
    """System for fitting to multipeak data
    Used to fit data with peaks discribed by lorentz, guassian or voigt profiles
    
    Attributes
    
    fitype:   List[strings N+v] can be "Gauss" "Lorentz" or "Voigt" corisponding to the fit of each peak
                                requires more entries for voigt due to its multiple characteristic width
    centroid: List[double  N] the position of the each of the centriod values for the data to be fit 
    FWHM:     List[double  N+v] the full width half maximum of the data
    y0:       List[double  N] the background level for the data to be fit
    roi:      List[double (2,1)] the limits for the fitting regon of interest 
    xs:       np.array[double M] the x space on which to plot the 
    """
    def __init__(self):
        self.fittype  = ['Gauss']
        self.centroid = [0]
        self.FWHM     = [10] 
        self.y0       = [0]
        self.roi      = [0,50]
        self.x  = np.linspace(0,50,1000)
        
    def getargs(self):
        """grabs the arguments for use in fitters internal and external"""
        args=[]
        j =0
        for i,type in enumerate(self.fittype):
            if type=='Voigt':
                args.append(self.centroid[i])
                args.append(self.FWHM[j])
                j+= 1
                args.append(self.FWHM[j])
                args.append(self.y0[i])
            if (type =='Gauss')^(type =="Lorentz"):
                args.append(self.centroid[i])
                args.append(self.FWHM[j])
                args.append(self.FWHM[j])
                args.append(self.y0[i])
            j+= 1
        return args

                
    def voigt():
        
    def guass():
    
    def lorentz():      

    def fitfunct(self,x,args):
        results = np.zeros(len(x))
        for i,type in enumerate(self.fittype):
            if (type=='Voigt'):
            
            if (type =='Gauss'):
                
            if (type =="Lorentz"):
        
        return results
                
        
    def fitto(self,data):
        