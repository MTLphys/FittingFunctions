import numpy as np 

def epszzOliva(hs,Rc,a_s,a_c):
    Ac =  Rc**2
    As = (hs+Rc)**2-Rc**2
    F = As/(Ac+As)
    f = (a_s-a_c)/a_c
    return F*f 

def epszzPrete(hs,Rc,hcap,a_s,a_c):
    F_top =hs**2/Rc**2 +2*hs/Rc  
    F_bottom = (hs/Rc)**2+ (hcap/Rc)**2+ 2*hcap/Rc*(1+hs/Rc)+2*((hs)/(Rc))+1
    f = (a_s-a_c)/a_c
    #print(f)
    #f=5.15e-4
    F=  F_top/F_bottom
    return F*f 


def DEcbhhS(a,d,nu,epszz):
    H= (1-2*nu)/3 
    return (3*a*H + (np.sqrt(3)/2)*d*(1-H))*epszz
def DEcblhS(a,d,nu,epszz,D0):
    H=(1-2*nu)/3
    return (3*a*H - (np.sqrt(3)/2)*d*(1-H))*epszz-3*(d**2)*(1-H)/(2*D0)*(epszz**2)


def DEcbhhP(a,d,nu,epszz):
    return epszz*(a*(1-2*nu) + (1/np.sqrt(3))*d*(1-nu))
def DEcblhP(a,d,nu,epszz,D0):
    return epszz*(a*(1-2*nu) - (1/np.sqrt(3))*d*(1-nu))-(epszz**2)*(d**2)*(1-nu)/(D0)

def DEcbhhO(epszz,A):
    return epszz*A
def DEcblhO(epszz,A):
    return epszz*A

#strain parameters of GaAs
a = -8.6 # Hydrostatic Deformation potential signorello
d = -5.2  # Shear deformation potentail signorello
D0 = .34 # spin orbit band splitting ioffe 
nu = .16 # poisson Ratio sigornello 
H =  (1 - 2*nu)/3


#dimensions of Nanowires
D=341e-9 #nanowire diameter 
Rc =75e-9 #core radius
hcap = 30e-9#cap Thickness 
hs = (D-2*Rc-2*hcap)/2#shell thickness 
X = .80 #Ga concentration in shell 
print(hs)

#lattice constants of cores and shells
aGaAs   = 5.65325 #A lattice constant in for GaAs Ioffe
aAlAs = 5.66164 #A lattice constant in pure AlAs

a_c = aGaAs 
a_s = X*aGaAs+(1-X)*aAlAs

import matplotlib.pyplot as plt 


#sigornello Calculations
fig, ax = plt.subplots(3)
#set up data for plots vs the oliva paper 
epszz = np.linspace(-.02,.03)
Ehh =  1.482+DEcbhhS(a,d,nu,epszz)
Elh =  1.482+DEcblhS(a,d,nu,epszz,D0)
ax[0].plot(epszz,Ehh)
ax[0].plot(epszz,Elh)

np.savetxt('EnergyForSigornello.csv',np.transpose(np.vstack((epszz,Ehh,Elh))),delimiter=';')

Rc  =  40e-9 #core radius
hcap = 3e-9#cap Thickness 
hs = 30e-9#shell thickness 
X= .7
a_s = X*aGaAs+(1-X)*aAlAs
epszzSig = epszzPrete(hs,Rc,hcap,a_s,a_c)
Ehh = 1.482+ DEcblhP(a,d,nu,epszzSig,D0)


print("estimate values Ehh for Sigornello paper",Ehh)
#Olivia Calculations
rc = 50e-9
rnw= 80e-9
rs = rnw-rc

a_s = .90*aGaAs+(1-.90)*aAlAs
epszz10 = epszzOliva(rs,rc,a_s,a_c)

a_s = .70*aGaAs+(1-.70)*aAlAs
epszz30 = epszzOliva(rs,rc,a_s,a_c)

print("estimate values epszz for Olivia paper(10%,30%)",epszz10*1e-4,epszz30*1e-4)

x = np.linspace(0.0,0.4,50)
a_s = x*aAlAs+(1-x)*aGaAs
epszzx = epszzOliva(rs,rc,a_s,a_c)
EhhO =   -DEcbhhO(epszzx,56)*1e2
EhhS =   DEcbhhS(a,d,nu,epszzx)*1e3
EhhP =   DEcbhhP(a,d,nu,epszzx)*1e3


ax[1].plot(x,EhhO)
ax[1].plot(x,EhhS)
ax[1].plot(x,EhhP)

np.savetxt('EnergyForOliva.csv',np.transpose(np.vstack((x,EhhO,EhhS,EhhP))),delimiter=';')
print("strainrate = ",DEcbhhP(a,d,nu,.1))
print("strainrate off by = ",1+.85/DEcbhhP(a,d,nu,.1))


#Prete Calculations
Eg = 1.510
a = -8.6 # Hydrostatic Deformation potential signorello
d = -5.2  # Shear deformation potentail signorello
D0 = .34 # spin orbit band splitting ioffe 
nu = .16 # poisson Ratio sigornello 
Rc =31.5e-9
hs = np.linspace(0,10*Rc,100)
hcap=Rc*.2
X = .33
a_s= X*aGaAs+(1-X)*aAlAs
a_c=aGaAs
epszzhs = epszzPrete(hs,Rc,hcap,a_s,a_c)
Ehh = Eg + DEcbhhP(a,d,nu,epszzhs)
Elh = Eg + DEcblhP(a,d,nu,epszzhs,D0)
ax[2].plot(hs/Rc,Ehh)
ax[2].plot(hs/Rc,Elh)
np.savetxt('EnergyshiftForPrete.csv',np.transpose(np.vstack((hs/Rc,Ehh,Elh))),delimiter=';')
plt.show()

epszz = epszzPrete()
print('')