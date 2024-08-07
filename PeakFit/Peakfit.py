import numpy as np 
import matplotlib.pyplot as plt  
from scipy import optimize as o
from scipy import special as s 

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
        self.y0       = 0
        self.A        = [1]
        self.roi      = [0,50]
        
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
                args.append(self.A[i])
            if (type =='Gauss')^(type =="Lorentz"):
                args.append(self.centroid[i])
                args.append(self.FWHM[j])
                args.append(self.A[i])
            j+= 1
        args.append(self.y0)
        return args
    
    def setargs(self,args):
        j=0
        k=0
        self.y0= args[-1]
        #print('arguments in set args',args)
        for i,type in enumerate(self.fittype):
            if type=='Voigt':
                self.centroid[i]= args[k]
                self.FWHM[j]    = args[k+1]
                j+=1
                self.FWHM[j]    = args[k+2]
                self.A[i]       = args[k+3]
                k+= 4 
            if (type =='Gauss')^(type =="Lorentz"):
                self.centroid[i] =  args[k]
                self.FWHM[j]     = args[k+1]
                self.A[i]       = args[k+2]
                k+= 3
            j+=1
            
    def voigt(x,A,FWHM1,FWHM2,centroid):
        sig = FWHM1/(2*np.sqrt(2*np.log(2)))
        gamma = FWHM2/2
        return A*s.voigt_profile(x-centroid,sig,gamma)/s.voigt_profile(0,sig,gamma)
        
    def gauss(x,A,FWHM,centroid):
        sig = FWHM/(2*np.sqrt(2*np.log(2)))
        return A*np.exp(-((x-centroid)**2/(2*sig**2)))
    
    def lorentz(x,A,FWHM,centroid):  
        gamma = FWHM/2    
        return A*1/((1)+((x-centroid)/gamma)**2)
    def pltfunct(self,x):
        j = 0
        results = np.zeros_like(x)
        for i,type in enumerate(self.fittype):
            print(type)
            if (type=='Voigt'):
                results += Fitter.voigt(x,self.A[i],self.FWHM[j],self.FWHM[j+1],self.centroid[i]) 
                j+=1
            if (type =='Gauss'):
                results += Fitter.gauss(x,self.A[i],self.FWHM[j],self.centroid[i])
            if (type =="Lorentz"):
                results += Fitter.lorentz(x,self.A[i],self.FWHM[j],self.centroid[i]) 
            j+=1
        return results + self.y0
    def fitfunct(self,x,*args):
        results = np.zeros(len(x))
        j= 0
        self.setargs(args)
        for i,type in enumerate(self.fittype):
            if (type=='Voigt'):
                results += Fitter.voigt(x,self.A[i],self.FWHM[j],self.FWHM[j+1],self.centroid[i]) 
                j+=1
            if (type =='Gauss'):
                results += Fitter.gauss(x,self.A[i],self.FWHM[j],self.centroid[i])
            if (type =="Lorentz"):
                results += Fitter.lorentz(x,self.A[i],self.FWHM[j],self.centroid[i]) 
            j+=1
        return results+self.y0
        
    def fit1d(self,x,y):
        pinit = self.getargs()
        print('fitparameters before fitting',pinit)
        poutp = o.curve_fit(self.fitfunct,x,y,p0=pinit)
        print('fit parameters after fit',poutp[0])
        self.setargs(poutp[0])
        
                    
        
