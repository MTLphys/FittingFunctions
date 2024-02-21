import numpy as np 
import matplotlib.pyplot as plt 
import scipy.constants as c

#define some convience functions and variables
e  = c.e 
h = c.h 
hbar = c.hbar
pi = np.pi 
th =np.heaviside

Ehh = 1.5139    #energy of the heavy hole 
Elh = Ehh+.0032 #energy of the light hole
Eib = Ehh-.0125 #energy of the impurity bound exciton 
#assign exciation 
P1 =1 #mW power of first pulse in milliwats 
P2 =1 #mW" " second pulse ""
P3 =1 #mW" " third pulse  ""

u1 = np.asarray([1,0])# polarization vector 1
u2 = np.asarray([1,0])# polarization vector 2 
u3 = np.asarray([1,0])# polarization vector 3 
#conjugate polarization vectors 
u1cj = np.conj(u1)# conj polarization vector 1
u2cj = np.conj(u2)# conj polarization vector 2 
u3cj = np.conj(u3)# conj polarization vector 3

#assign grating parameters  
N1 = P1*2.0e15#1/cm^3 exciton density 1
N2 = P2*2.0e15#1/cm^3 exciton density 2
N3 = P3*2.0e15#1/cm^3 exciton density 3
A = N1+N2+N3 

N1eh = P1*1.0e15#1/cm^3 electron hole pair density 1
N2eh = P2*1.0e15#1/cm^3 electron hole pair density 2
N3eh = P3*1.0e15#1/cm^3 electron hole pair density 3
Aeh= N1eh+N2eh+N3eh

sigx = .0002 #1/cm^3 exciton/exciton scattering parameter for GaAs From literature
sigeh = .000 #1/cm^3 exciton/eh scattering parameter for GaAs From literature

#exciton Dephasing parameters 
g21 =    6.5e12 # 1/s heavy hole dephasing
g61 =    6.5e12 # 1/s light hole dephasing 
gi12=    6.5e12 # 1/s impurity dephasing 

#exciton coupling constants
mu21 =  1.0* np.asarray([1, 0])# Dipole coupling hh (non-directional)
mu61 =  0.0* np.asarray([1, 0])# Dipole coupling lh (non-directional)
mui21=  0.0* np.asarray([1, 0])# impurity coupling (non-directional)

#conjugate coupling constants 
mu21cj = np.conj(mu21)# conjugate hh 
mu61cj = np.conj(mu61)# conjugate lh 
mui21cj  = np.conj(mui21)# conjugate impurity

#angular frequency definitions for available states
Om21 =  Ehh/h*e*2.0*pi - 1.0j*g21
Om61 =  Elh/h*e*2.0*pi - 1.0j*g61
Omi21 = Eib/h*e*2.0*pi - 1.0j*g61

#EID frequencys( compensating for scattering effects)
Om21t = Om21 -1.0j*(sigx*A+sigeh*Aeh)
Om61t = Om61 -1.0j*(sigx*A+sigeh*Aeh)
Omi21t = Omi21 -1.0j*(sigx*A+sigeh*Aeh)

#conjugated angular frequencies
Om21cj = np.conj(Om21)
Om61cj = np.conj(Om61)
Omi21cj = np.conj(Omi21)

#conjugated EID frequencies 
Om21cj = np.conj(Om21t)
Om61cj = np.conj(Om61t)
Omi21cj = np.conj(Omi21t)

PSF = 1  #psf scaling 
EID = 0 #eid scaling  
tau13 =7.0e-12#31 delay 
K = 1e-117 # defining scaling factor
T1 = 100e-12 #1/sradiative Exciton Lifetime 

def Pnw(t,tau):
    #define PSF term
    P1 =PSF*(th(tau,1.0)*np.exp(1.0j*(Om21cj+1j*(sigx*N1+sigeh*N1eh))*tau)*
        np.exp(-(tau13-tau)/T1)*np.exp(1j*Om21t*(tau13-t))*
        2*mu21cj[:,np.newaxis]*mu21.dot(u3)*mu21.dot(u2)*mu21cj.dot(u1cj)+
        
        th(-tau,1.0)*np.exp(1.0j*(Om21-1j*(sigx*N2+sigeh*N2eh))*tau)*
        np.exp(-(tau13)/T1)*np.exp(1j*Om21t*(tau13-t))*
        2*mu21cj[:,np.newaxis]*mu21.dot(u3)*mu21cj.dot(u1cj)*mu21.dot(u2)+
        
        th(tau,1.0)*np.exp(1.0j*(Om61cj+1j*(sigx*N1+sigeh*N1eh))*tau)*
        np.exp(-(tau13-tau)/T1)*np.exp(1j*Om21t*(tau13-t))*
        mu21cj[:,np.newaxis]*mu21.dot(u3)*mu61.dot(u2)*mu61cj.dot(u1cj)+
        
        th(-tau,1.0)*np.exp(1.0j*(Om61-1j*(sigx*N2+sigeh*N2eh))*tau)*
        np.exp(-(tau13)/T1)*np.exp(1j*Om21t*(tau13-t))*
        mu21cj[:,np.newaxis]*mu21.dot(u3)*mu61cj.dot(u1cj)*mu61.dot(u2)+
        
        th(tau,1.0)*np.exp(1.0j*(Om61cj+1j*(sigx*N1+sigeh*N1eh))*tau)*
        np.exp(-(tau13-tau)/T1)*np.exp(1j*Om61t*(tau13-t))*
        2*mu61cj[:,np.newaxis]*mu61.dot(u3)*mu61.dot(u2)*mu61cj.dot(u1cj)+
        
        th(-tau,1.0)*np.exp(1.0j*(Om61-1j*(sigx*N2+sigeh*N2eh))*tau)*
        np.exp(-(tau13)/T1)*np.exp(1j*Om61t*(tau13-t))*
        2*mu61cj[:,np.newaxis]*mu61.dot(u3)*mu61cj.dot(u1cj)*mu61.dot(u2)+
        
        th(tau,1.0)*np.exp(1.0j*(Om21cj+1j*(sigx*N1+sigeh*N1eh))*tau)*
        np.exp(-(tau13-tau)/T1)*np.exp(1j*Om61t*(tau13-t))*
        mu61cj[:,np.newaxis]*mu61.dot(u3)*mu21.dot(u2)*mu21cj.dot(u1cj)+
        
        th(-tau,1.0)*np.exp(1.0j*(Om21-1j*(sigx*N2+sigeh*N2eh))*tau)*
        np.exp(-(tau13)/T1)*np.exp(1j*Om61t*(tau13-t))*
        mu61cj[:,np.newaxis]*mu61.dot(u3)*mu21cj.dot(u1cj)*mu21.dot(u2))
    
    #define EID term
    P2 = A*sigx*EID*(th(tau,1.0)*np.exp(1.0j*(Om21cj+1.0j*(sigx*N1+sigeh*N1eh))*tau)*
        np.exp(-(tau13-tau)/T1)*T1*(1-np.exp(-((t-tau13)/T1)))*np.exp(1.0j*Om21t*(tau13-t))*
        mu21cj[:,np.newaxis]*mu21.dot(u3)*mu21.dot(u2)*mu21cj.dot(u1cj)+

        th(-tau,1.0)*np.exp(1.0j*(Om21-1.0j*(sigx*N2+sigeh*N2eh))*tau)*
        np.exp(-(tau13)/T1)*T1*(1-np.exp(-((t-tau13)/T1)))*np.exp(1.0j*Om21t*(tau13-t))*
        mu21cj[:,np.newaxis]*mu21.dot(u3)*mu21cj.dot(u1)*mu21.dot(u2)+

        th(tau,1.0)*np.exp(1.0j*(Om61cj+1.0j*(sigx*N1+sigeh*N1eh))*tau)*
        np.exp(-(tau13-tau)/T1)*T1*(1-np.exp(-((t-tau13)/T1)))*np.exp(1.0j*Om21t*(tau13-t))*
        mu21cj[:,np.newaxis]*mu21.dot(u3)*mu61.dot(u2)*mu21cj.dot(u1cj)+

        th(-tau,1.0)*np.exp(1.0j*(Om61-1.0j*(sigx*N2+sigeh*N2eh))*tau)*
        np.exp(-(tau13)/T1)*T1*(1-np.exp(-((t-tau13)/T1)))*np.exp(1.0j*Om21t*(tau13-t))*
        mu21cj[:,np.newaxis]*mu21.dot(u3)*mu61cj.dot(u1cj)*mu61.dot(u2)+
        
        th(tau,1.0)*np.exp(1.0j*(Om61cj+1.0j*(sigx*N1+sigeh*N1eh))*tau)*
        np.exp(-(tau13-tau)/T1)*T1*(1-np.exp(-((t-tau13)/T1)))*np.exp(1.0j*Om61t*(tau13-t))*
        mu61cj[:,np.newaxis]*mu61.dot(u3)*mu61.dot(u2)*mu61cj.dot(u1cj)+

        th(-tau,1.0)*np.exp(1.0j*(Om61-1.0j*(sigx*N2+sigeh*N2eh))*tau)*
        np.exp(-(tau13)/T1)*T1*(1-np.exp(-((t-tau13)/T1)))*np.exp(1.0j*Om61t*(tau13-t))*
        mu61cj[:,np.newaxis]*mu61.dot(u3)*mu61cj.dot(u1cj)*mu61.dot(u2)+
        
        th(tau,1.0)*np.exp(1.0j*(Om21cj+1.0j*(sigx*N1+sigeh*N1eh))*tau)*
        np.exp(-(tau13-tau)/T1)*T1*(1-np.exp(-((t-tau13)/T1)))*np.exp(1.0j*Om61t*(tau13-t))*
        mu61cj[:,np.newaxis]*mu61.dot(u3)*mu21.dot(u2)*mu21cj.dot(u1cj)+

        th(-tau,1.0)*np.exp(1.0j*(Om21-1.0j*(sigx*N2+sigeh*N2eh))*tau)*
        np.exp(-(tau13)/T1)*T1*(1-np.exp(-((t-tau13)/T1)))*np.exp(1.0j*Om61t*(tau13-t))*
        mu61cj[:,np.newaxis]*mu61.dot(u3)*mu21cj.dot(u1cj)*mu21.dot(u2)
        )
    return K*-1.0j*(A/(hbar**3))*th(t-tau13,1.0)*(P1+P2)
def Pimp(t,tau):
    P1 = PSF*(th(tau,1.0)*np.exp(1.0j*(Omi21cj+1.0j*(sigx*N1+sigeh*N1eh))*tau)*
        np.exp(-(tau13-tau)/T1)*np.exp(1.0j*Omi21t*(tau13-t))*
        2*mui21cj[:,np.newaxis]*mui21.dot(u3)*mui21.dot(u2)*mui21cj.dot(u1)+
            
        th(-tau,1.0)*np.exp(1.0j*(Omi21-1.0j*(sigx*N2+sigeh*N2eh))*tau)*
        np.exp(-(tau13)/T1)*np.exp(1.0j*Omi21t*(tau13-t))*
        2*mui21cj[:,np.newaxis]*mui21.dot(u3)*mui21cj.dot(u1cj)*mui21.dot(u2)
    )
    P2= EID*A*sigx*(th(tau,1.0)*np.exp(1.0j*(Omi21cj+1.0j*(sigx*N1+sigeh*N1eh))*tau)*
        np.exp(-(tau13-tau)/T1)*T1*(1-np.exp(-(t-tau13)/T1))*np.exp(1.0j*Omi21*(tau13-t))*
        mui21cj[:,np.newaxis]*mui21.dot(u3)*mui21.dot(u2)*mui21cj.dot(u1cj)+
        
        th(-tau,1.0)*np.exp(1.0j*(Om21-1.0j*(sigx*N2+sigeh*N2eh))*tau)*
        np.exp(-(tau13)/T1)*T1*(1-np.exp(-(t-tau13)/T1))*np.exp(1.0j*Omi21*(tau13-t))*
        mui21cj[:,np.newaxis]*mui21.dot(u3)*mui21cj.dot(u1cj)*mui21.dot(u2)+

        th(tau,1.0)*np.exp(1.0j*(Omi21cj+1.0j*(sigx*N1+sigeh*N1eh))*tau)*
        np.exp(-(tau13-tau)/T1)*T1*(1-np.exp(-(t-tau13)/T1))*np.exp(1.0j*Omi21*(tau13-t))*
        mui21cj[:,np.newaxis]*mui21.dot(u3)*mui21.dot(u2)*mui21cj.dot(u1cj)+

        th(-tau,1.0)*np.exp(1.0j*(Omi21-1.0j*(sigx*N2+sigeh*N2eh))*tau)*
        np.exp(-(tau13)/T1)*T1*(1-np.exp(-(t-tau13)/T1))*np.exp(1.0j*Omi21*(tau13-t))*
        mui21cj[:,np.newaxis]*mui21.dot(u3)*mui21cj.dot(u1cj)*mui21.dot(u2)   
    )

    return -K*1.0j*A/hbar**3*th(t-tau13,1.0)*(P1+P2)

def Pfull(t,tau):
    return Pnw(t,tau)+Pimp(t,tau)
fig,ax = plt.subplots(2)
taudef = 0e-12
tdef = 7e-12

tau = np.linspace(-1.5e-12, 1.5e-12, 1000)
t   = np.linspace( 7.0e-12, 8.2e-12, 1000)
print(Pfull(tdef,0))
print(np.abs(Pfull(tdef,tau))[0])

ax[0].plot(tau,np.abs(Pfull(tdef,tau))[0])
ax[1].plot(t, np.abs(Pfull(t,taudef)[0]))
plt.show()