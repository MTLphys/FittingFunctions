import scipy.constants as c
import numpy as np 
import matplotlib.pyplot as plt 


u1  = np.asarray([1.0,0.0])
u2  = np.asarray([1.0,0.0])
u3  = np.asarray([1.0,0.0])

N1  = 1.0e15
N2  = 1.0e15
N3  = 1.0e15

N1eh = 1.0e15
N2eh = 1.0e15
N3eh = 1.0e15 

gm21  = 5.0e12
gm61  = 5.0e12
gmi21 = 5.0e12

sigx  = .0002
sigeh = .0000

A = N1 + N2 + N3 
Aeh = N1eh + N2eh + N3eh
K = 1e-117
T1 = 100.0e-12

O21 = c.e*1.5139/c.hbar -1.0j*gm21 
O61 = c.e*(1.5139+.0025)/c.hbar -1.0j*gm61
Oi21 = c.e*(1.5139-0.023)/c.hbar -1.0j*gmi21

Ot21 = O21-1.0j*(sigx*A+sigeh*Aeh)
Ot61 = O61-1.0j*(sigx*A+sigeh*Aeh)
Oti21 = Oi21-1.0j*(sigx*A+sigeh*Aeh)

O = np.asarray([O21,O61,Oi21])
Ot = np.asarray([Ot21,Ot61,Oti21])

tau13 = 7.0e-12 

mu21 = 1.0*np.asarray([1.0,0.0])
mu61 = 0.75*np.asarray([1.0,0.0])
mui21 = 0.4*np.asarray([1.0,0.0])

mu = np.vstack((mu21,mu61,mui21))
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
def aterm(t,tau,i,j):
    wf =np.heaviside(tau,1.0)*np.exp(1.0j*(np.conj(O[i])+1.0j*(sigx*N1+sigeh*N1eh))*tau)
    lf = np.exp(-(tau13-tau)/T1)
    tf = np.exp(1.0j*Ot[j]*(tau13-t))
    cp = np.conj(mu[j])*mu[j].dot(u3)*mu[i].dot(u2)*np.conj(mu[i]).dot(np.conj(u1))
    return vecdim(wf*lf*tf,cp)
def bterm(t,tau,i,j):
    wf =np.heaviside(tau,1.0)*np.exp(1.0j*(O[i]-1.0j*(sigx*N2+sigeh*N2eh))*tau)
    lf = np.exp(-(tau13)/T1)
    tf = np.exp(1.0j*Ot[j]*(tau13-t))
    cp = np.conj(mu[j])*mu[j].dot(u3)*np.conj(mu[i]).dot(np.conj(u1))*mu[i].dot(u2)
    return vecdim(wf*lf*tf,cp)
def cterm(t,tau,i,j):
    wf =np.heaviside(tau,1.0)*np.exp(1.0j*(np.conj(O[i])+1.0j*(sigx*N1+sigeh*N1eh))*tau)
    lf = np.exp(-(tau13-tau)/T1)
    ei =T1*(1-np.exp((tau13-t)/T1))
    tf = np.exp(1.0j*Ot[j]*(tau13-t))
    cp = np.conj(mu[j])*mu[j].dot(u3)*mu[i].dot(u2)*np.conj(mu[i]).dot(np.conj(u1))
    return vecdim(wf*ei*lf*tf,cp)
def dterm(t,tau,i,j):
    wf =np.heaviside(tau,1.0)*np.exp(1.0j*(O[i]-1.0j*(sigx*N2+sigeh*N2eh))*tau)
    lf = np.exp(-(tau13)/T1)
    ei =T1*(1-np.exp((tau13-t)/T1))
    tf = np.exp(1.0j*Ot[j]*(tau13-t))
    cp = np.conj(mu[j])*mu[j].dot(u3)*np.conj(mu[i]).dot(np.conj(u1))*mu[i].dot(u2)
    return vecdim(wf*ei*lf*tf,cp)
def Pnw(t,tau):
    sig=0
    for i in range(2):
        sig = sig + PSF*(aterm(t,tau,i,i)+bterm(t,tau,i,i))
        for j in range(3):
            if(j<2):
                sig = sig+PSF*(aterm(t,tau,i,j)+bterm(t,tau,i,j))
            sig =sig+EID*A*sigx*(cterm(t,tau,i,j)+dterm(t,tau,i,j))
    return 1*K*-1.0j*A/(c.hbar**3)*np.heaviside(t-tau13,1.0)*sig
def Pimp(t,tau):
    sig = PSF*(2*aterm(t,tau,2,2)+2*bterm(t,tau,2,2))
    for i in range(3):
        sig=sig + EID*A*sigx*(cterm(t,tau,i,2)+dterm(t,tau,i,2))
    return 1*K*-1.0j*A/(c.hbar**3)*np.heaviside(t-tau13,1.0)*sig
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