
import scipy.constants as c
import numpy as np 
import matplotlib.pyplot as plt 
def fitfunctt(t,Ex1,Ex2,EB,ol21,ol61,oli21,oli61,gm21,gm61,gmi21,gmi61,K,Ad,tau=0e-12):
    u1  = np.asarray([1.0,0.0])#polarization state pulse 1
    u2  = np.asarray([1.0,0.0])#polarization state pulse 2
    u3  = np.asarray([1.0,0.0])#polarization state pulse 3

    N1  = 1.0e15 #Exciton density pulse 1
    N2  = 1.0e15 #Exciton density pulse 2
    N3  = 1.0e15 #Exciton density pulse 3

    N1eh = 1.0e15 #electron hole pair density pulse 1
    N2eh = 1.0e15 #electron hole pair density pulse 2
    N3eh = 1.0e15 #electron hole pair density pulse 3

    N1i = 1.0e15#Impurity bound Exciton density pulse 1 
    N2i = 1.0e15#Impurity bound Exciton density pulse 2
    N3i = 1.0e15#Impurity bound Exciton density pulse 3



    sigx   = 0.00016 #exciton exciton scattering Parameter
    sigeh  = 0.00050 #exciton electron/hole scattering Parameter 
    sigxi  = 0.00016 #exciton impurity scattering Parameter
    sigi   = 0.00000 #impurity impurity scattering Parameter
    sigieh = 0.00050 #impurity electron/hole scattering Parameter
    sigv = np.asarray([sigx,sigxi,sigi])
    sigvchk = np.asarray(['sigx','sigxi','sigi'])

    A =   N1   +N2   +N3         # total exciton density
    Ai =  N1i  +N2i  +N3i         # total impurity bound exciton density
    Aeh = N1eh +N2eh +N3eh # total electron hole density

    
    T1 = 100.0e-12 #exciton Lifetime
    Ti1= 100.0e-12 #impurity bound exciton Lifetime
    

    E0 =c.e*(1.5095)/c.hbar 
    E21 =  E0 
    E61 =  E0-Ex1 
    Ei21 = E0+EB 
    Ei61 = E0-Ex2+EB 
    O21 =  E21 -1.0j*gm21  #Base damped oscillator frequency  
    O61 =  E61 -1.0j*gm61  #Base damped oscillator frequency  
    Oi21 = Ei21 -1.0j*gmi21 #Base damped oscillator frequency
    Oi61 = Ei61 -1.0j*gmi61 #Base damped oscillator frequency

    Ot21 =  O21 -1.0j*(sigx*A+ sigxi*Ai+sigeh*Aeh) #Base Eid oscillator frequency
    Ot61 =  O61 -1.0j*(sigx*A+ sigxi*Ai+sigeh*Aeh)#Base Eid oscillator frequency
    Oti21 = Oi21-1.0j*(sigi*Ai+sigxi*A+ sigieh*Aeh)#Base Eid oscillator frequency
    Oti61 = Oi61-1.0j*(sigi*Ai+sigxi*A+ sigieh*Aeh)#Base Eid oscillator frequency

    O  = np.asarray([ O21, O61, Oi21, Oi61]) 
    Ot = np.asarray([Ot21,Ot61,Oti21,Oti61]) 

    tau13 = 4.0e-12 
    mu21 = ol21*np.asarray([1.0,0.0])   #dipole moment of hh exciton state 
    mu61 = ol61*np.asarray([1.0,0.0])   #dipole moment of lh exciton state 
    mui21 = oli21*np.asarray([1.0,0.0]) #dipole moment of hh bound exciton state 
    mui61 = oli61*np.asarray([1.0,0.0]) #dipole moment of lh bound exciton state 


    mu = np.vstack((mu21,mu61,mui21,mui61))
    print(mu[0])
    imp = 1 
    nw  = 1
    PSF = 1
    EID = 1
    def vecdim(a,b): 
        '''
            Element wise multiplier for vector time variables
            takes 
            a =  n length base vector element 
            b =  m length time axis over which to multiply 
            returns 
            nxm vector field
        
        '''
        v =  []
        for bi in b:
            v.append(a*bi)
        return np.asarray(v)
    def aterm(t,tau,i,j,k):
        '''psf term for optical bloch equation solution 
        takes 
        t = the time delay for time resolving reference pulse
        tau = delay between first and second pulses 
        i = index of primary term for the 
        j = index of secondary term for the 
        k = scattering term for free exciton(0) or impurity(1)
        return 
        vector field 
        '''
        #print('PSF Term',i,j,k)
        if(k==0):
            wf = np.heaviside(tau,0.5)*np.exp(1.0j*(np.conj(O[i])+1.0j*(sigx*N1+sigxi*N1i+sigeh*N1eh))*tau)
            lf = np.exp(-(tau13-tau)/T1)
        else:
            wf =  np.heaviside(tau,0.5)*np.exp(1.0j*(np.conj(O[i])+1.0j*(sigi*N1i+sigxi*N1+sigieh*N1eh))*tau)
            lf = np.exp(-(tau13-tau)/Ti1)
        tf = np.exp(1.0j*Ot[j]*(tau13-t))
        cp = np.conj(mu[j])*mu[j].dot(u3)*mu[i].dot(u2)*np.conj(mu[i]).dot(np.conj(u1))
        return vecdim(wf*lf*tf,cp)
    def bterm(t,tau,i,j,k):
        """
        complex conjugate psf term for optical bloch equation solution 
        takes 
        t = the time delay for time resolving reference pulse
        tau = delay between first and second pulses 
        i = index of primary term for the 
        j = index of secondary term for the 
        k = scattering term for free exciton(0) or impurity(1)
        return 
        vector field
        """
        if(k==0):    
            wf = np.heaviside(-tau,0.5)*np.exp(1.0j*(O[i]-1.0j*(sigx*N2+sigxi*N2i+sigeh*N2eh))*tau)
            lf = np.exp(-(tau13)/T1)
        else:
            wf = np.heaviside(-tau,0.5)*np.exp(1.0j*(O[i]-1.0j*(sigi*N2i+sigxi*N2+sigieh*N2eh))*tau)
            lf = np.exp(-(tau13)/Ti1)
        tf = np.exp(1.0j*Ot[j]*(tau13-t))
        cp = np.conj(mu[j])*mu[j].dot(u3)*np.conj(mu[i]).dot(np.conj(u1))*mu[i].dot(u2)
        return vecdim(wf*lf*tf,cp)
    def cterm(t,tau,i,j,k):
        """
        EID term for optical bloch equation solution 
        takes 
        t = the time delay for time resolving reference pulse
        tau = delay between first and second pulses 
        i = index of primary term for the 
        j = index of secondary term for the 
        k = scattering term for free exciton(0) or impurity(1)
        return 
        vector field
        """
        #print('EID term',i,j,k)
        if(k==0):
            wf = np.heaviside(tau,0.5)*np.exp(1.0j*(np.conj(O[i])+1.0j*(sigx*N1+sigxi*N1i+sigeh*N1eh))*tau)
            lf = np.exp(-(tau13-tau)/T1)
            ei = T1*(1-np.exp((tau13-t)/T1))
        else:
            wf = np.heaviside(tau,0.5)*np.exp(1.0j*(np.conj(O[i])+1.0j*(sigi*N1i+sigxi*N1+sigieh*N1eh))*tau)
            lf = np.exp(-(tau13-tau)/Ti1)
            ei = Ti1*(1-np.exp((tau13-t)/Ti1))
        tf = np.exp(1.0j*Ot[j]*(tau13-t))
        cp = np.conj(mu[j])*mu[j].dot(u3)*mu[i].dot(u2)*np.conj(mu[i]).dot(np.conj(u1))
        return vecdim(wf*ei*lf*tf,cp)
    def dterm(t,tau,i,j,k):
        """
        complex conjugate EID term for optical bloch equation solution 
        takes 
        t = the time delay for time resolving reference pulse
        tau = delay between first and second pulses 
        i = index of primary term for the 
        j = index of secondary term for the 
        k = scattering term for free exciton(0) or impurity(1) type term
        return 
        vector field
        """
        if(k==0):
            wf =np.heaviside(-tau,0.5)*np.exp(1.0j*(O[i]-1.0j*(sigx*N2+sigxi*N2i+sigeh*N2eh))*tau)
            lf = np.exp(-(tau13)/T1)
            ei =T1*(1-np.exp((tau13-t)/T1))
        else: 
            wf =np.heaviside(-tau,0.5)*np.exp(1.0j*(O[i]-1.0j*(sigi*N2+sigxi*N2i+sigieh*N2eh))*tau)
            lf = np.exp(-(tau13)/Ti1)
            ei =Ti1*(1-np.exp((tau13-t)/Ti1))
        tf = np.exp(1.0j*Ot[j]*(tau13-t))
        cp = np.conj(mu[j])*mu[j].dot(u3)*np.conj(mu[i]).dot(np.conj(u1))*mu[i].dot(u2)
        return vecdim(wf*ei*lf*tf,cp)
    def Pnw(t,tau):
        """
        Polarization of the general exciton states of the Ga As Coreshell Cap nanowires 
        including contributions from scattering on eh pairs, bound acceptor states, and 
        exciton-exciton scattering
        takes: 
        t = the time delay between t1 and t2 
        tau = the time delay between t3 and tref
        returns 
        nanowire exciton polarization 2xNtres 
        """
        sig=0
        for j in range(2):
            sig = sig + PSF*(aterm(t,tau,j,j,0)+bterm(t,tau,j,j,0))
            for i in range(4):
                if(i<2):
                    sig = sig+PSF*(aterm(t,tau,i,j,0)+bterm(t,tau,i,j,0))
                #sig =sig+EID*A*sigv[int(np.floor(i/2))]*(cterm(t,tau,i,j,int(np.floor(i/2)))+dterm(t,tau,i,j,int(np.floor(i/2))))
        return 1*K*-1.0j*A/(c.hbar**3)*np.heaviside(t-tau13,.5)*sig
    def Pimp(t,tau):
        """
        Polarization of the general exciton states of the Ga As Coreshell Cap nanowires 
        including contributions from scattering on eh pairs, bound acceptor states, and 
        exciton-exciton scattering
        takes: 
        t = the time delay between t1 and t2 
        tau = the time delay between t3 and tref
        returns 
        impurity stat  polarization 2xNtres 
        """
        sig=0
        for j in range(2):
            sig = sig + 2*PSF*(aterm(t,tau,j+2,j+2,1)+bterm(t,tau,j+2,j+2,1))
            #for i in range(4):
                #if(i>1):
                #    sig = sig+PSF*(aterm(t,tau,i,j+2,1)+bterm(t,tau,i,j+2,1))
                #sig =sig+EID*Ai*sigv[int(np.floor(i/2))]*(cterm(t,tau,i,j+2,int(np.floor(i/2)))+dterm(t,tau,i,j+2,int(np.floor(i/2))))
        return 1*K*(-1.0j)*Ai/(c.hbar**3)*np.heaviside(t-tau13,.5)*sig
    def P(t,tau):
        return nw*Pnw(t,tau)+imp*Pimp(t,tau)
    def FWM(t,tau):
        return np.abs(P(t,tau))[0]
    def dirac(t,tau):
        a = 0+(np.abs(t-tau13)<.03e-12)&(np.abs(tau)<.03e-12)
        return a
    return FWM(t,tau)+Ad*dirac(t,tau)
def fitfuncttau(tau,Ex1,Ex2,EB,ol21,ol61,oli21,oli61,gm21,gm61,gmi21,gmi61,K,t=7.1e-12):
    u1  = np.asarray([1.0,0.0])#polarization state pulse 1
    u2  = np.asarray([1.0,0.0])#polarization state pulse 2
    u3  = np.asarray([1.0,0.0])#polarization state pulse 3

    N1  = 1.0e15 #Exciton density pulse 1
    N2  = 1.0e15 #Exciton density pulse 2
    N3  = 1.0e15 #Exciton density pulse 3

    N1eh = 1.0e15 #electron hole pair density pulse 1
    N2eh = 1.0e15 #electron hole pair density pulse 2
    N3eh = 1.0e15 #electron hole pair density pulse 3

    N1i = 1.0e15#Impurity bound Exciton density pulse 1 
    N2i = 1.0e15#Impurity bound Exciton density pulse 2
    N3i = 1.0e15#Impurity bound Exciton density pulse 3



    sigx   = 0.00016 #exciton exciton scattering Parameter
    sigeh  = 0.00050 #exciton electron/hole scattering Parameter 
    sigxi  = 0.00016 #exciton impurity scattering Parameter
    sigi   = 0.00000 #impurity impurity scattering Parameter
    sigieh = 0.00050 #impurity electron/hole scattering Parameter
    sigv = np.asarray([sigx,sigxi,sigi])
    sigvchk = np.asarray(['sigx','sigxi','sigi'])

    A =   N1   +N2   +N3         # total exciton density
    Ai =  N1i  +N2i  +N3i         # total impurity bound exciton density
    Aeh = N1eh +N2eh +N3eh # total electron hole density

    
    T1 = 100.0e-12 #exciton Lifetime
    Ti1= 100.0e-12 #impurity bound exciton Lifetime
    

    E0 =c.e*(1.5065)/c.hbar 
    E21 =  E0 
    E61 =  E0+Ex1 
    Ei21 = E0+EB 
    Ei61 = E0+Ex2+EB 
    O21 =  E21 -1.0j*gm21  #Base damped oscillator frequency  
    O61 =  E61 -1.0j*gm61  #Base damped oscillator frequency  
    Oi21 = Ei21 -1.0j*gmi21 #Base damped oscillator frequency
    Oi61 = Ei61 -1.0j*gmi61 #Base damped oscillator frequency

    Ot21 =  O21 -1.0j*(sigx*A+ sigxi*Ai+sigeh*Aeh) #Base Eid oscillator frequency
    Ot61 =  O61 -1.0j*(sigx*A+ sigxi*Ai+sigeh*Aeh)#Base Eid oscillator frequency
    Oti21 = Oi21-1.0j*(sigi*Ai+sigxi*A+ sigieh*Aeh)#Base Eid oscillator frequency
    Oti61 = Oi61-1.0j*(sigi*Ai+sigxi*A+ sigieh*Aeh)#Base Eid oscillator frequency

    O  = np.asarray([ O21, O61, Oi21, Oi61]) 
    Ot = np.asarray([Ot21,Ot61,Oti21,Oti61]) 

    tau13 = 4.0e-12 
    mu21 = ol21*np.asarray([1.0,0.0])   #dipole moment of hh exciton state 
    mu61 = ol61*np.asarray([1.0,0.0])   #dipole moment of lh exciton state 
    mui21 = oli21*np.asarray([1.0,0.0]) #dipole moment of hh bound exciton state 
    mui61 = oli61*np.asarray([1.0,0.0]) #dipole moment of lh bound exciton state 


    mu = np.vstack((mu21,mu61,mui21,mui61))
    print(mu[0])
    imp = 1 
    nw  = 1
    PSF = 1
    EID = 1
    def vecdim(a,b): 
        '''
            Element wise multiplier for vector time variables
            takes 
            a =  n length base vector element 
            b =  m length time axis over which to multiply 
            returns 
            nxm vector field
        
        '''
        v =  []
        for bi in b:
            v.append(a*bi)
        return np.asarray(v)
    def aterm(t,tau,i,j,k):
        '''psf term for optical bloch equation solution 
        takes 
        t = the time delay for time resolving reference pulse
        tau = delay between first and second pulses 
        i = index of primary term for the 
        j = index of secondary term for the 
        k = scattering term for free exciton(0) or impurity(1)
        return 
        vector field 
        '''
        if(k==0):
            wf = np.heaviside(tau,1.0)*np.exp(1.0j*(np.conj(O[i])+1.0j*(sigx*N1+sigxi*N1i+sigeh*N1eh))*tau)
            lf = np.exp(-(tau13-tau)/T1)
        else:
            wf =  np.heaviside(tau,1.0)*np.exp(1.0j*(np.conj(O[i])+1.0j*(sigi*N1i+sigxi*N1+sigieh*N1eh))*tau)
            lf = np.exp(-(tau13-tau)/Ti1)
        tf = np.exp(1.0j*Ot[j]*(tau13-t))
        cp = np.conj(mu[j])*mu[j].dot(u3)*mu[i].dot(u2)*np.conj(mu[i]).dot(np.conj(u1))
        return vecdim(wf*lf*tf,cp)
    def bterm(t,tau,i,j,k):
        """
        complex conjugate psf term for optical bloch equation solution 
        takes 
        t = the time delay for time resolving reference pulse
        tau = delay between first and second pulses 
        i = index of primary term for the 
        j = index of secondary term for the 
        k = scattering term for free exciton(0) or impurity(1)
        return 
        vector field
        """
        if(k==0):    
            wf = np.heaviside(-tau,1.0)*np.exp(1.0j*(O[i]-1.0j*(sigx*N2+sigxi*N2i+sigeh*N2eh))*tau)
            lf = np.exp(-(tau13)/T1)
        else:
            wf = np.heaviside(-tau,1.0)*np.exp(1.0j*(O[i]-1.0j*(sigi*N2i+sigxi*N2+sigieh*N2eh))*tau)
            lf = np.exp(-(tau13)/Ti1)
        tf = np.exp(1.0j*Ot[j]*(tau13-t))
        cp = np.conj(mu[j])*mu[j].dot(u3)*np.conj(mu[i]).dot(np.conj(u1))*mu[i].dot(u2)
        return vecdim(wf*lf*tf,cp)
    def cterm(t,tau,i,j,k):
        """
        EID term for optical bloch equation solution 
        takes 
        t = the time delay for time resolving reference pulse
        tau = delay between first and second pulses 
        i = index of primary term for the 
        j = index of secondary term for the 
        k = scattering term for free exciton(0) or impurity(1)
        return 
        vector field
        """
        if(k==0):
            wf = np.heaviside(tau,1.0)*np.exp(1.0j*(np.conj(O[i])+1.0j*(sigx*N1+sigxi*N1i+sigeh*N1eh))*tau)
            lf = np.exp(-(tau13-tau)/T1)
            ei = T1*(1-np.exp((tau13-t)/T1))
        else:
            wf = np.heaviside(tau,1.0)*np.exp(1.0j*(np.conj(O[i])+1.0j*(sigi*N1i+sigxi*N1+sigieh*N1eh))*tau)
            lf = np.exp(-(tau13-tau)/Ti1)
            ei = Ti1*(1-np.exp((tau13-t)/Ti1))
        tf = np.exp(1.0j*Ot[j]*(tau13-t))
        cp = np.conj(mu[j])*mu[j].dot(u3)*mu[i].dot(u2)*np.conj(mu[i]).dot(np.conj(u1))
        return vecdim(wf*ei*lf*tf,cp)
    def dterm(t,tau,i,j,k):
        """
        complex conjugate EID term for optical bloch equation solution 
        takes 
        t = the time delay for time resolving reference pulse
        tau = delay between first and second pulses 
        i = index of primary term for the 
        j = index of secondary term for the 
        k = scattering term for free exciton(0) or impurity(1) type term
        return 
        vector field
        """
        if(k==0):
            wf =np.heaviside(-tau,1.0)*np.exp(1.0j*(O[i]-1.0j*(sigx*N2+sigxi*N2i+sigeh*N2eh))*tau)
            lf = np.exp(-(tau13)/T1)
            ei =T1*(1-np.exp((tau13-t)/T1))
        else: 
            wf =np.heaviside(-tau,1.0)*np.exp(1.0j*(O[i]-1.0j*(sigi*N2+sigxi*N2i+sigieh*N2eh))*tau)
            lf = np.exp(-(tau13)/Ti1)
            ei =Ti1*(1-np.exp((tau13-t)/Ti1))
        tf = np.exp(1.0j*Ot[j]*(tau13-t))
        cp = np.conj(mu[j])*mu[j].dot(u3)*np.conj(mu[i]).dot(np.conj(u1))*mu[i].dot(u2)
        return vecdim(wf*ei*lf*tf,cp)
    def Pnw(t,tau):
        """
        Polarization of the general exciton states of the Ga As Coreshell Cap nanowires 
        including contributions from scattering on eh pairs, bound acceptor states, and 
        exciton-exciton scattering
        takes: 
        t = the time delay between t1 and t2 
        tau = the time delay between t3 and tref
        returns 
        nanowire exciton polarization 2xNtres 
        """
        sig=0
        for j in range(2):
            sig = sig + PSF*(aterm(t,tau,j,j,0)+bterm(t,tau,j,j,0))
            for i in range(4):
                if(i<2):
                    sig = sig+PSF*(aterm(t,tau,i,j,0)+bterm(t,tau,i,j,0))
                sig =sig+EID*A*sigv[int(np.floor(i/2))]*(cterm(t,tau,i,j,int(np.floor(i/2)))+dterm(t,tau,i,j,int(np.floor(i/2))))
        return 1*K*-1.0j*A/(c.hbar**3)*np.heaviside(t-tau13,0.5)*sig
    def Pimp(t,tau):
        """
        Polarization of the general exciton states of the Ga As Coreshell Cap nanowires 
        including contributions from scattering on eh pairs, bound acceptor states, and 
        exciton-exciton scattering
        takes: 
        t = the time delay between t1 and t2 
        tau = the time delay between t3 and tref
        returns 
        impurity stat  polarization 2xNtres 
        """
        sig=0
        for j in range(2):
            sig = sig + PSF*(aterm(t,tau,j+2,j+2,0)+bterm(t,tau,j+2,j+2,0))
            for i in range(4):
                if(i>1):
                    sig = sig+PSF*(aterm(t,tau,i,j+2,0)+bterm(t,tau,i,j+2,0))
                sig =sig+EID*Ai*sigv[int(np.floor(i/2))]*(cterm(t,tau,i,j+2,int(np.floor(i/2)))+dterm(t,tau,i,j+2,int(np.floor(i/2))))
        return 1*K*(-1.0j)*Ai/(c.hbar**3)*np.heaviside(t-tau13,0.5)*sig
    def P(t,tau):
        return nw*Pnw(t,tau)+imp*Pimp(t,tau)
    def FWM(t,tau):
        return np.abs(P(t,tau))[0]
    return FWM(t,tau)
