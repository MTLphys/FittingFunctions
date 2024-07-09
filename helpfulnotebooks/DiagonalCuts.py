import numpy as np 
import matplotlib.pyplot as plt 
xmin = -3
xmax = 3 
ymin = -1
ymax = 3
Nx = 121
Ny = 81

import glob as gb 
directory = 'C:/Users/mattl/OneDrive - University of Cincinnati/Desktop/Lab Data/HFWMV2 CS NW InitialTest/830nm/Longboi/'
endmoniker = ''
filetype = '.csv'
outputfile = 'Results/'
integration= .15
files = gb.glob(directory+'*'+endmoniker+filetype)

print(files)
file = files[1]

for file in files : 
    print(file)
    filetag = file[len(directory):-4]
    contents = np.loadtxt(file,delimiter=';')
    t12 = contents[0,1:]
    t3ref = contents[1:,0]

#output = np.vstack((t3ref[1:],t3refcut[1:])).transpose()
#np.savetxt(directory+outputfile+filetag+'T3refcut'+filetype,output,delimiter=';')
#output = np.vstack((t12[1:],t12cut[1:])).transpose()
#np.savetxt(directory+outputfile+filetag+'T12cut'+filetype,output,delimiter=';')



    Nx = len(t12)
    Ny = len(t3ref)
    x = t12
    y = t3ref
    Z = np.zeros((Ny,Nx))
    ZM = np.zeros((Ny,Nx))
    R = np.zeros((Ny,Nx))
    XX,YY= np.meshgrid(x,y)
    Z = contents[1:,1:]
    Ns =5 # number of slices
    dx = x[1]-x[0]
    dy = y[1]-y[0]
    dr = np.sqrt(dx**2+dy**2)
    Rf = np.array([])
    ZMf = np.array([])
    for i in range(Ns):
        pos=dy*(i)-(Ns-Ns%2)*dy/2
        print(pos)
        R[np.abs(XX-YY-pos)<0.01*dx]= ((XX[np.abs(XX-YY-pos)<0.01])+(YY[np.abs(XX-YY-pos)<0.01]))#np.sign(XX[np.abs(XX-YY-pos)<0.01]+YY[np.abs(XX-YY-pos)<0.01]) 
        ZM[np.abs(XX-YY-pos)<0.01*dx]=Z[np.abs(XX-YY-pos)<0.01*dx] 
        Rf=np.append(Rf,R[np.abs(XX-YY-pos)<0.01*dx].flatten())
        ZMf = np.append(ZMf,ZM[np.abs(XX-YY-pos)<0.01*dx].flatten())

    fig,ax = plt.subplots(2,2)
    cset1 = ax[0,0].imshow(np.flipud(Z),extent=(xmin,xmax,ymin,ymax))
    fig.colorbar(cset1)
    cset2 = ax[1,0].imshow(np.flipud(R),extent=(xmin,xmax,ymin,ymax))
    fig.colorbar(cset2)
    cset2 = ax[0,1].imshow(np.flipud(ZM),extent=(xmin,xmax,ymin,ymax))
    Rind = np.argsort(Rf)
    fig.colorbar(cset2)
    r = Rf[Rind]
    z = ZMf[Rind]
    ax[1,1].plot(r,z)
    ax[1,1].set_yscale('log')
    output= np.vstack((r,z)).transpose()
    np.savetxt(directory+outputfile+filetag+'diagonalcut'+filetype,output,delimiter=';')
    

plt.show()
