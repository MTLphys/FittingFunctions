
import scipy.constants as c
import numpy as np 
import matplotlib.pyplot as plt 
exwl = 815 #exciting wavelength

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

gm21  = 5.0e12#hh Exciton dephasing rate /s
gm61  = 5.0e12#lh Exciton dephasing rate /s
gmi21 = 5.0e12#Impurity bound hh Exciton dephasing rate /s
gmi21 = 5.0e12#Impurity bound lh Exciton dephasing rate /s

sigx   = 0.00016 #exciton exciton scattering Parameter
sigeh  = 0.00050 #exciton electron/hole scattering Parameter 
sigxi  = 0.00016 #exciton impurity scattering Parameter
sigi   = 0.00000 #impurity impurity scattering Parameter
sigieh = 0.00050 #impurity electron/hole scattering Parameter

A = N1 + N2 + N3         # total exciton density
Ai = N1i+N2i+N3i         # total impurity bound exciton density
Aeh = N1eh + N2eh + N3eh # total electron hole density

K = 1e-117     #
T1 = 100.0e-12 #exciton Lifetime
Ti1= 100.0e-12 #impurity bound exciton Lifetime
O21 = c.e*1.5139/c.hbar -1.0j*gm21 
O61 = c.e*(1.5139+.0025)/c.hbar -1.0j*gm61
Oi21 = c.e*(1.5139-0.023)/c.hbar -1.0j*gmi21
Oi61 = c.e*(1.5139-0.023+.0025)/c.hbar -1.0j*gmi21

Ot21 = O21-1.0j*(sigx*A+sigxi*Ai+sigeh*Aeh)
Ot61 = O61-1.0j*(sigx*A+sigxi*Ai+sigeh*Aeh)
Oti21 = Oi21-1.0j*(sigi*A+sigxi*Ai+sigieh*Aeh)
Oti61 = Oi61-1.0j*(sigi*A+sigxi*Ai+sigieh*Aeh)

O  = np.asarray([ O21, O61, Oi21, Oi61])
Ot = np.asarray([Ot21,Ot61,Oti21,Oti61])

tau13 = 7.0e-12 

ol21 = 1.0#overlap of hh exciton state
ol61 = 0.8#overlap of lh exciton state
oli21 = 0.5 #overlap of hh bound exciton state
oli61 = 0.4 #overlap of lh bound exciton state


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
    v =  []
    for bi in b:
        v.append(a*bi)
    return np.asarray(v)
def aterm(t,tau,i,j,k):
    if(k==0):
        wf =np.heaviside(tau,1.0)*np.exp(1.0j*(np.conj(O[i])+1.0j*(sigx*N1+sigxi*N1i+sigeh*N1eh))*tau)
        lf = np.exp(-(tau13-tau)/T1)
    else:
        wf =np.heaviside(tau,1.0)*np.exp(1.0j*(np.conj(O[i])+1.0j*(sigi*N1i+sigxi*N1+sigieh*N1eh))*tau)
        lf = np.exp(-(tau13-tau)/Ti1)
    tf = np.exp(1.0j*Ot[j]*(tau13-t))
    cp = np.conj(mu[j])*mu[j].dot(u3)*mu[i].dot(u2)*np.conj(mu[i]).dot(np.conj(u1))
    return vecdim(wf*lf*tf,cp)
def bterm(t,tau,i,j,k):
    if(k==0):    
        wf =np.heaviside(tau,1.0)*np.exp(1.0j*(O[i]-1.0j*(sigx*N2+sigxi*N2i+sigeh*N2eh))*tau)
        lf = np.exp(-(tau13)/T1)
    if(k==1):
        wf =np.heaviside(tau,1.0)*np.exp(1.0j*(O[i]-1.0j*(sigi*N2i+sigxi*N2+sigieh*N2eh))*tau)
        lf = np.exp(-(tau13)/Ti1)
    tf = np.exp(1.0j*Ot[j]*(tau13-t))
    cp = np.conj(mu[j])*mu[j].dot(u3)*np.conj(mu[i]).dot(np.conj(u1))*mu[i].dot(u2)
    return vecdim(wf*lf*tf,cp)
def cterm(t,tau,i,j,k):
    if(k==0):
        wf =np.heaviside(tau,1.0)*np.exp(1.0j*(np.conj(O[i])+1.0j*(sigx*N1+sigxi*N1i+sigeh*N1eh))*tau)
        lf = np.exp(-(tau13-tau)/T1)
        ei =T1*(1-np.exp((tau13-t)/T1))
    else:
        wf =np.heaviside(tau,1.0)*np.exp(1.0j*(np.conj(O[i])+1.0j*(sigx*N1+sigxi*N1i+sigeh*N1eh))*tau)
        lf = np.exp(-(tau13-tau)/T1)
        ei =T1*(1-np.exp((tau13-t)/T1))
    tf = np.exp(1.0j*Ot[j]*(tau13-t))
    cp = np.conj(mu[j])*mu[j].dot(u3)*mu[i].dot(u2)*np.conj(mu[i]).dot(np.conj(u1))
    return vecdim(wf*ei*lf*tf,cp)
def dterm(t,tau,i,j,k):
    if(k==0):
        wf =np.heaviside(tau,1.0)*np.exp(1.0j*(O[i]-1.0j*(sigx*N2+sigxi*N2i+sigeh*N2eh))*tau)
        lf = np.exp(-(tau13)/T1)
        ei =T1*(1-np.exp((tau13-t)/T1))
    else: 
        wf =np.heaviside(tau,1.0)*np.exp(1.0j*(O[i]-1.0j*(sigx*N2+sigxi*N2i+sigeh*N2eh))*tau)
        lf = np.exp(-(tau13)/T1)
        ei =T1*(1-np.exp((tau13-t)/T1))
    tf = np.exp(1.0j*Ot[j]*(tau13-t))
    cp = np.conj(mu[j])*mu[j].dot(u3)*np.conj(mu[i]).dot(np.conj(u1))*mu[i].dot(u2)
    return vecdim(wf*ei*lf*tf,cp)
def Pnw(t,tau):
    sig=0
    for i in range(2):
        sig = sig + PSF*(aterm(t,tau,i,i,0)+bterm(t,tau,i,i,0))
        for j in range(3):
            if(j<2):
                sig = sig+PSF*(aterm(t,tau,i,j,0)+bterm(t,tau,i,j,0))
            sig =sig+EID*A*sigx*(cterm(t,tau,i,j)+dterm(t,tau,i,j))
    return 1*K*-1.0j*A/(c.hbar**3)*np.heaviside(t-tau13,1.0)*sig
def Pimp(t,tau):
    sig = PSF*(2*aterm(t,tau,2,2,1)+2*bterm(t,tau,2,2,1))
    for i in range(3):
        sig=sig + EID*Ai*sigx*(cterm(t,tau,i,2)+dterm(t,tau,i,2))
    return 1*K*-1.0j*Ai/(c.hbar**3)*np.heaviside(t-tau13,1.0)*sig
def P(t,tau):
    return nw*Pnw(t,tau)+imp*Pimp(t,tau)
def FWM(t,tau):
    return np.abs(P(t,tau))[0]
t = np.linspace(7e-12,10e-12,1000)
Pv = np.abs(P(t,0e-12))[0]

fig,axs = plt.subplots(2)
axs[0].set_yscale('log')
axs[1].set_yscale('log')
axs[0].plot(t*1e12-7,Pv)
tau = np.linspace(-1.5e-12,1.5e-12,1000)