
import numpy as np
import matplotlib.pyplot as plt

'''

'''


 ### This is the Gaussian data smoothing function I wrote ###  

def smoothListGaussian(list,degree=5):  

     window=degree*2-1  

     weight=np.array([1.0]*window)  

     weightGauss=[]  

     for i in range(window):  

         i=i-degree+1  

         frac=i/float(window)  

         gauss=1/(np.exp((4*(frac))**2))  

         weightGauss.append(gauss)  

     weight=np.array(weightGauss)*weight  

     smoothed=[0.0]*(len(list)-window)  

     for i in range(len(smoothed)):  

         smoothed[i]=sum(np.array(list[i:i+window])*weight)/sum(weight)  

     return smoothed  

 



#fe_obswave,fe_obsflux   = np.loadtxt('SN2011fe_110907_uvopt_obswave_obsflux.dat',  dtype=float,usecols=(0,1),unpack=True)
fe_obswave,fe_obsflux   = np.loadtxt('SN2011fe_peak072_uvopt_obswave_obsflux.dat',  dtype=float,usecols=(0,1),unpack=True)
by_obswave,by_obsflux   = np.loadtxt('SN2011by_110509_uvopt_obswave_obsflux.dat',  dtype=float,usecols=(0,1),unpack=True)
erp_obswave,erp_obsflux = np.loadtxt('SN2017erp_170702_uvopt_obswave_obsflux.dat', dtype=float,usecols=(0,1),unpack=True)

#plt.plot(erp_obswave,erp_obsflux)
#plt.show(block=True)

fe_smoothflux=smoothListGaussian(fe_obsflux,degree=3)
by_smoothflux=smoothListGaussian(by_obsflux,degree=3)
erp_smoothflux=smoothListGaussian(erp_obsflux,degree=3)

fe_smoothwave=np.array(smoothListGaussian(fe_obswave,degree=3))
by_smoothwave=np.array(smoothListGaussian(by_obswave,degree=3))
erp_smoothwave=np.array(smoothListGaussian(erp_obswave,degree=3))


#print(len(fe_obsflux))


#plt.plot(by_smoothwave,by_smoothflux)
#plt.plot(fe_smoothwave,fe_smoothflux)
#plt.show(block=True)

fe_z=  0.00080400 
by_z=  0.00284300 
erp_z= 0.00617400 


fe_restwave=   fe_smoothwave/(1.0+fe_z)
by_restwave=   by_smoothwave/(1.0+by_z)
erp_restwave= erp_smoothwave/(1.0+erp_z)

#
#print(fe_restwave[0])
#print(by_restwave[0])
#print(erp_restwave[0])

#1778.96970835 # for 110907
# 1598 for 11fe_peak072
#1595.46409558
#1590.182215

#print(fe_restwave[-1])
#print(by_restwave[-1])
#print(erp_restwave[-1])

#24965.3278764
#7977.32047788
#5565.63775252



# sp_ea = np.interp(wavez,filter_lambda,filter_area) ### spectrum effective area         
w1_wave,w1_ea   = np.loadtxt('../filters/UVW1_B11.txt',  dtype=float,usecols=(0,1),unpack=True)


#use w1_wave with wavelength from 1600 to 8000 as the region of interest

wavelength10_array=w1_wave

#erp_normfactor=mean(erp_obsflux(where(erp_restwave gt 5400 and erp_restwave lt 5560))
#erp_normfactor=np.mean(erp_obsflux(where(erp_restwave gt 5400 and erp_restwave lt 5560))

norm_indices=np.logical_and(wavelength10_array < 5560,wavelength10_array > 5400).nonzero()


fe_restflux10 = np.interp(wavelength10_array,fe_restwave, fe_smoothflux)
by_restflux10 = np.interp(wavelength10_array,by_restwave, by_smoothflux) 
erp_restflux10 = np.interp(wavelength10_array,erp_restwave, erp_smoothflux) 


fe_restflux10_norm =  fe_restflux10  / np.mean(fe_restflux10[norm_indices])
by_restflux10_norm =  by_restflux10  / np.mean(by_restflux10[norm_indices])
erp_restflux10_norm = erp_restflux10 / np.mean(erp_restflux10[norm_indices])


feby=fe_restflux10_norm / by_restflux10_norm
feerp=fe_restflux10_norm / erp_restflux10_norm

byfe=by_restflux10_norm / fe_restflux10_norm
erpfe=erp_restflux10_norm / fe_restflux10_norm



optical_indices=(wavelength10_array > 5555).nonzero()
faruv_indices=(wavelength10_array < 1675).nonzero()


feby[optical_indices]=1.0
feerp[optical_indices]=1.0
byfe[optical_indices]=1.0
erpfe[optical_indices]=1.0

feby[faruv_indices]=3.0
feerp[faruv_indices]=10.0
byfe[faruv_indices]=0.3
erpfe[faruv_indices]=0.1



#plt.plot(wavelength10_array,fe_restflux10_norm)
#plt.plot(wavelength10_array,feby)
#plt.plot(wavelength10_array,feby)
#plt.axis([0,2000,0,9])
#plt.show(block=True)

#


plt.plot(wavelength10_array,feby)
plt.plot(wavelength10_array,feerp)
plt.axis([1500,6000,0,9])
plt.show(block=True)


plt.plot(wavelength10_array,byfe)
plt.plot(wavelength10_array,erpfe)
plt.axis([1500,6000,0,1.1])
plt.show(block=True)








#print("done")
