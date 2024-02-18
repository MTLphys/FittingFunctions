import matplotlib.pyplot as plt 
import numpy as np 
import glob as gb 
directory = 'C:/Users/mattl/OneDrive - University of Cincinnati/Desktop/Lab Data/HFWMV2 CS NW InitialTest/LongCollections/'
endmoniker = '100'
filetype = '.csv'
integration= .15
files = gb.glob(directory+'*'+endmoniker+filetype)
fig,axs = plt.subplots(2)
axs[0].set_yscale('log')
axs[1].set_yscale('log')
print(files)
for file in files : 
    filetag = file[len(directory):-4]
    contents = np.loadtxt(file,delimiter=';')
    indicest12= np.where(np.abs(contents[:,0])<=integration)[0]
    t12cut = np.sum(contents[indicest12[1:],:],axis=0)/(len(indicest12)-1)
    t12 = contents[0,:]
    t12cut = t12cut-np.min(t12cut[1:])+1e-7
    axs[0].plot(t12[1:],t12cut[1:],label = filetag[-6:])
    axs[0].legend()
    indicest3ref= np.where(np.abs(contents[0,:])<=integration)[0]
    t3refcut = np.sum(contents[:,indicest3ref[1:]],axis=1)/(len(indicest3ref)-1)
    t3refcut = t3refcut-np.min(t3refcut[1:])+1e-7
    t3ref = contents[:,0]
    axs[1].plot(t3ref[1:],t3refcut[1:],label = filetag[0:6])
    output = np.vstack((t3ref[1:],t3refcut[1:])).transpose()
    np.savetxt(directory+filetag+'T3refcut'+filetype,output,delimiter=';')
    output = np.vstack((t12[1:],t12cut[1:])).transpose()
    np.savetxt(directory+filetag+'T12cut'+filetype,output,delimiter=';')