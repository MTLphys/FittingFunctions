def epszzPrete(hs,Rc,hcap,a_s,a_c):
    """Z strain coeffecient for GaAs core in AlGaAs nanowires based on publication by Paola Prete  
    
    Args:
        hs (double): Height of the shell
        Rc (double): Radius of the core
        hcap (double): Height of the cap
        a_s (double): shell lattice constant
        a_c (double): core lattice constant

    Returns:
        double: epszz for core
    """
    F_top =hs**2/Rc**2 +2*hs/Rc  
    F_bottom = (hs/Rc)**2+ (hcap/Rc)**2+ 2*hcap/Rc*(1+hs/Rc)+2*((hs)/(Rc))+1
    f = (a_s-a_c)/a_c #lattice parameter
    print(f)
    #f=5.15e-4
    F=  F_top/F_bottom #area ratio
    return F*f #total strain

import numpy as np 
X=.2
T= 69.5 #calculated wire tempurature at 60uW excitation power



#Calculate the lattice coeffiecients for core and shell
CAlGaAs= (5.73-.53*X)*1e-6 #/C Coeffiecient of thermal expansion 
aAlGaAsRT = (5.6533+.0078*X)*1e-10 #m RT lattice constant AlGaAs
aAlGaAs=aAlGaAsRT*(1+CAlGaAs*(T-300))#m at temp lattice constant
CGaAs = (5.73)*1e-6 #/C Coeffiecient of thermal expansion
aGaAsRT = (5.6533)*1e-10 #m RTlattice constant GaAs
aGaAs = aGaAsRT*(1+CGaAs*(T-300))#m at temp lattice constant GaAs

#Core shell cap parameters
a_s = aAlGaAs
a_c = aGaAs

D=341e-9 #nanowire diameter 
Rc =75e-9 #core radius
hcap = 30e-9#cap Thickness 
hs = (D-2*Rc-2*hcap)/2#shell thickness 
print(hs)


#calculate the strain for the nanowire:
epszz = epszzPrete(hs,Rc,hcap,a_s,a_c)
print("percentage strain for Core of nanowires",epszz*1e2)
print("energy shift for the nanowire",epszz*1e2*85.3)
Eg=1.519-5.405e-4*(T**2)/(T+204)# (eV)
print('binding energy',Eg)
print('exciton energy',Eg-epszz*8.53-.0042)