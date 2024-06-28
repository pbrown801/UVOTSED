from pyphot import Filter
from pyphot import unit
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pandas as pd
import extinction
import time
import sys
import os
import seaborn as sns
import matplotlib.colors as colors


#Set up swift filters

Udata = pd.read_csv('filters/U_P08.txt',delim_whitespace=True,comment='#')
Udata.columns = ['Wavelength','Area']

Bdata = pd.read_csv('filters/B_P08.txt',delim_whitespace=True,comment='#')
Bdata.columns = ['Wavelength','Area']

Vdata = pd.read_csv('filters/V_P08.txt',delim_whitespace=True,comment='#')
Vdata.columns = ['Wavelength','Area']


W2data = pd.read_csv('filters/UVW2_B11.txt',delim_whitespace=True,comment='#')
W2data.columns = ['Wavelength','Area']

M2data = pd.read_csv('filters/UVM2_B11.txt',delim_whitespace=True,comment='#')
M2data.columns = ['Wavelength','Area']

W1data = pd.read_csv('filters/UVW1_B11.txt',delim_whitespace=True,comment='#')
W1data.columns = ['Wavelength','Area']


## SWIFT_UVOT_Filter
SWIFT_UVOT_B = Filter(Bdata.Wavelength,Bdata.Area,name = 'SWIFT_UVOT_B',dtype = 'photon',unit = 'Angstrom')
SWIFT_UVOT_U = Filter(Udata.Wavelength,Udata.Area,name = 'SWIFT_UVOT_U',dtype = 'photon',unit = 'Angstrom')
SWIFT_UVOT_V = Filter(Vdata.Wavelength,Vdata.Area,name = 'SWIFT_UVOT_V',dtype = 'photon',unit = 'Angstrom')
SWIFT_UVOT_UVM2 = Filter(M2data.Wavelength,M2data.Area,name = 'SWIFT_UVOT_UVM2',dtype = 'photon',unit = 'Angstrom')
SWIFT_UVOT_UVW1 = Filter(W1data.Wavelength,W1data.Area,name = 'SWIFT_UVOT_UVW1',dtype = 'photon',unit = 'Angstrom')
SWIFT_UVOT_UVW2 = Filter(W2data.Wavelength,W2data.Area,name = 'SWIFT_UVOT_UVW2',dtype = 'photon',unit = 'Angstrom')
FilterVec = [SWIFT_UVOT_B,SWIFT_UVOT_U,SWIFT_UVOT_V,SWIFT_UVOT_UVM2,SWIFT_UVOT_UVW1,SWIFT_UVOT_UVW2]


def kfun(wave,flux):
    f_b = FilterVec[0].get_flux(wave,flux)
    f_b = -2.5*np.log10(f_b) - FilterVec[0].Vega_zero_mag
    
    f_u = FilterVec[1].get_flux(wave,flux)
    f_u = -2.5*np.log10(f_u) - FilterVec[1].Vega_zero_mag
    
    f_v = FilterVec[2].get_flux(wave,flux)
    f_v = -2.5*np.log10(f_v) - FilterVec[2].Vega_zero_mag
    
    f_m2 = FilterVec[3].get_flux(wave,flux)
    f_m2 = -2.5*np.log10(f_m2) - FilterVec[3].Vega_zero_mag
    
    f_w1 = FilterVec[4].get_flux(wave,flux)
    f_w1 = -2.5*np.log10(f_w1) - FilterVec[4].Vega_zero_mag
    
    f_w2 = FilterVec[5].get_flux(wave,flux)
    f_w2 = -2.5*np.log10(f_w2) - FilterVec[5].Vega_zero_mag
    
    return [f_w2,f_m2,f_w1,f_u,f_b] - f_v

# mods = np.asarray(['SN2011by_110509_uvopt_obswave_obsflux.dat','PTF11kly_20110903.obs.dat',
#                      'SN2011fe_110907_uvopt_obswave_obsflux.dat','SN2011iv_111211_full.dat',
#                     'SN2013dy_visit3_hst.flm','SN2015F_20150326_full.dat',
#                    'SN2016ccj_160514_uvopt_obswave_obsflux.dat',
#                    'SN2017erp_170702_uvopt_obswave_obsflux.dat',
#                    'ASASSN14lp_visit3_hst.flm','Hsiao18_.dat',])

mods = np.asarray(['spectra/ASASSN-14LP_peak_11fe_appended.dat','spectra/SN1992A_UV.dat','spectra/SN2011by_peak_11fe_appended.dat',
	'spectra/SN2011fe_peak_11fe_appended.dat','spectra/SN2011iv_peak_11fe_appended.dat','spectra/SN2015F_peak_11fe_appended.dat',
	'spectra/SN2016ccj_peak_11fe_appended.dat','spectra/SN2017erp_peak_11fe_appended.dat','spectra/SN2021fxy_peak_11fe_appended.dat',
	'spectra/SN2022hrs_peak_11fe_appended.dat','spectra/SN2013dy_peak_11fe_appended.dat'])

def modgen(modl):
    
    modname = modl.split('_')[0]
    Tdata = pd.read_csv('kcortemplates/'+modl,delim_whitespace=True,comment='#')
    if len(Tdata.columns) == 3:
        Tdata.columns = ['wavelength','flux','ferr']

    else:
        Tdata.columns = ['wavelength','flux']
        
    SpectrumFlux = np.asarray(Tdata.flux)/max(Tdata.flux)
    SlambdaWave = np.asarray(Tdata.wavelength)
        



    def reddening(AV, RV):
        return extinction.apply(extinction.ccm89(SlambdaWave,AV,RV),SpectrumFlux)
    
    plt.figure(1,figsize=(8,6),dpi=250)
    divnorm = colors.TwoSlopeNorm(vmin=-0.5, vcenter=0.0, vmax=0.5)
    mycmap = plt.get_cmap('coolwarm')(np.linspace(0.0,1.0,11))
    avs = np.round(np.linspace(-0.5,0.5,11),2)
    avs[2] = 0.0
    colst = []
    colnm = ['UVW2-V','UVM2-V','UVW1-V','U-V','B-V']
    plt.subplot(212)
    for i in range(10):
        nf = reddening(avs[i],1.8)
        for k in range(len(nf)):
            if nf[k] < 0.0:
                nf[k] = 0.0
        plt.plot(SlambdaWave, nf,label=str(avs[i]), color=mycmap[i])
        nucol = kfun(SlambdaWave,nf)
        colst.append(nucol)
        redata = {'wavelength':SlambdaWave,'flux':nf}
        nux = pd.DataFrame(data=redata)
        nux.to_csv('outputs/models/'+modname+'_1.8_'+str(avs[i])+'_.csv',index=False)
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Flux (Arbitrary Units)')
    
    #plt.legend(title='AV =',loc='upper left',ncol=2)
    plt.xlim(1500,6000)
    plt.subplot(211)
    for i in range(len(colst)):
        plt.scatter(colnm,colst[i],color=mycmap[i],label=str(avs[i]))
    plt.ylabel('Color')
    plt.legend(title='AV =',loc = 'upper right',ncol=2)
    plt.title(modname+', RV = 1.8')
    plt.tight_layout()
    plt.savefig('outputs/models/'+modname+'_1.8_.png',facecolor='white')
    plt.close()
    
    plt.figure(2,figsize=(8,6),dpi=250)
    divnorm = colors.TwoSlopeNorm(vmin=-0.5, vcenter=0.0, vmax=0.5)
    mycmap = plt.get_cmap('coolwarm')(np.linspace(0.0,1.0,11))
    avs = np.round(np.linspace(-0.5,0.5,11),2)
    avs[2] = 0.0
    colst = []
    colnm = ['UVW2-V','UVM2-V','UVW1-V','U-V','B-V']
    plt.subplot(212)
    for i in range(10):
        nf = reddening(avs[i],2.7)
        for k in range(len(nf)):
            if nf[k] < 0.0:
                nf[k] = 0.0
        plt.plot(SlambdaWave, nf,label=str(avs[i]), color=mycmap[i])
        nucol = kfun(SlambdaWave,nf)
        colst.append(nucol)
        redata = {'wavelength':SlambdaWave,'flux':nf}
        nux = pd.DataFrame(data=redata)
        nux.to_csv('outputs/models/'+modname+'_2.7_'+str(avs[i])+'_.csv',index=False)
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Flux (Arbitrary Units)')

    #plt.legend(title='AV =',loc='upper left',ncol=2)
    plt.xlim(1500,6000)
    plt.subplot(211)
    for i in range(len(colst)):
        plt.scatter(colnm,colst[i],color=mycmap[i],label=str(avs[i]))
    plt.ylabel('Color')
    plt.legend(title='AV =',loc = 'upper right',ncol=2)
    plt.title(modname+', RV = 2.7')
    plt.tight_layout()
    plt.savefig('outputs/models/'+modname+'_2.7_.png',facecolor='white')
    plt.close()
    

    plt.figure(3,figsize=(8,6),dpi=250)
    divnorm = colors.TwoSlopeNorm(vmin=-0.5, vcenter=0.0, vmax=0.5)
    mycmap = plt.get_cmap('coolwarm')(np.linspace(0.0,1.0,11))
    avs = np.round(np.linspace(-0.5,0.5,11),2)
    avs[2] = 0.0
    colst = []
    colnm = ['UVW2-V','UVM2-V','UVW1-V','U-V','B-V']
    plt.subplot(212)
    for i in range(10):
        nf = reddening(avs[i],3.1)
        for k in range(len(nf)):
            if nf[k] < 0.0:
                nf[k] = 0.0
        plt.plot(SlambdaWave, nf,label=str(avs[i]), color=mycmap[i])
        nucol = kfun(SlambdaWave,nf)
        colst.append(nucol)
        redata = {'wavelength':SlambdaWave,'flux':nf}
        nux = pd.DataFrame(data=redata)
        nux.to_csv('outputs/models/'+modname+'_3.1_'+str(avs[i])+'_.csv',index=False)
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Flux (Arbitrary Units)')
    #plt.legend(title='AV =',loc='upper left',ncol=2)
    plt.xlim(1500,6000)
    plt.subplot(211)
    for i in range(len(colst)):
        plt.scatter(colnm,colst[i],color=mycmap[i],label=str(avs[i]))
    plt.ylabel('Color')
    plt.legend(title='AV =',loc = 'upper right',ncol=2)
    plt.title(modname+', RV = 3.1')
    plt.tight_layout()
    plt.savefig('outputs/models/'+modname+'_3.1_.png',facecolor='white')
    plt.close()
    
for t in range(len(mods)):
    modgen(mods[t])
#modgen(mods[5])


#Read in the UV model from Foley et al 2016
SlambdaData = pd.read_csv('spectra/uvmodel.data',delim_whitespace=True,comment='#')
SlambdaData.columns = ['wavelength','f_11','Slambda']
SlambdaWave = np.asarray(SlambdaData.wavelength).astype(float)
Slambda = np.asarray(SlambdaData.Slambda)
f_11 = np.asarray(SlambdaData.f_11)

f11Fun = interp1d(SlambdaWave,f_11,kind = 'cubic')
SlambdaFun = interp1d(SlambdaWave,Slambda ,kind = 'cubic')

#apply a stretch to the flux based on the estimated DB15



dms = np.round(np.linspace(0.8,2.1,27),2)

def fgen(dmb):
    SpectrumFlux = f11Fun(SlambdaWave) + SlambdaFun(SlambdaWave)*(dmb-1.1)
    modname = 'foley'+str(dmb)

    def reddening(AV, RV):
        return extinction.apply(extinction.ccm89(SlambdaWave,AV,RV),SpectrumFlux)
    
    plt.figure(1,figsize=(8,6),dpi=250)
    divnorm = colors.TwoSlopeNorm(vmin=-0.5, vcenter=0.0, vmax=0.5)
    mycmap = plt.get_cmap('coolwarm')(np.linspace(0.0,1.0,11))
    avs = np.round(np.linspace(-0.5,0.5,11),2)
    avs[2] = 0.0
    colst = []
    colnm = ['UVW2-V','UVM2-V','UVW1-V','U-V','B-V']
    plt.subplot(212)
    for i in range(10):
        nf = reddening(avs[i],1.8)
        for k in range(len(nf)):
            if nf[k] < 0.0:
                nf[k] = 0.0
        plt.plot(SlambdaWave, nf,label=str(avs[i]), color=mycmap[i])
        nucol = kfun(SlambdaWave,nf)
        colst.append(nucol)
        redata = {'wavelength':SlambdaWave,'flux':nf}
        nux = pd.DataFrame(data=redata)
        nux.to_csv('outputs/models/'+modname+'_1.8_'+str(avs[i])+'_.csv',index=False)
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Flux (Arbitrary Units)')
    
    #plt.legend(title='AV =',loc='upper left',ncol=2)
    plt.xlim(1500,6000)
    plt.subplot(211)
    for i in range(len(colst)):
        plt.scatter(colnm,colst[i],color=mycmap[i],label=str(avs[i]))
    plt.ylabel('Color')
    plt.legend(title='AV =',loc = 'upper right',ncol=2)
    plt.title(modname+', RV = 1.8')
    plt.tight_layout()
    plt.savefig('outputs/models/'+modname+'_1.8_.png',facecolor='white')
    plt.close()
    
    plt.figure(2,figsize=(8,6),dpi=250)
    divnorm = colors.TwoSlopeNorm(vmin=-0.5, vcenter=0.0, vmax=0.5)
    mycmap = plt.get_cmap('coolwarm')(np.linspace(0.0,1.0,11))
    avs = np.round(np.linspace(-0.5,0.5,11),2)
    avs[2] = 0.0
    colst = []
    colnm = ['UVW2-V','UVM2-V','UVW1-V','U-V','B-V']
    plt.subplot(212)
    for i in range(10):
        nf = reddening(avs[i],2.7)
        for k in range(len(nf)):
            if nf[k] < 0.0:
                nf[k] = 0.0
        plt.plot(SlambdaWave, nf,label=str(avs[i]), color=mycmap[i])
        nucol = kfun(SlambdaWave,nf)
        colst.append(nucol)
        redata = {'wavelength':SlambdaWave,'flux':nf}
        nux = pd.DataFrame(data=redata)
        nux.to_csv('outputs/models/'+modname+'_2.7_'+str(avs[i])+'_.csv',index=False)
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Flux (Arbitrary Units)')

    #plt.legend(title='AV =',loc='upper left',ncol=2)
    plt.xlim(1500,6000)
    plt.subplot(211)
    for i in range(len(colst)):
        plt.scatter(colnm,colst[i],color=mycmap[i],label=str(avs[i]))
    plt.ylabel('Color')
    plt.legend(title='AV =',loc = 'upper right',ncol=2)
    plt.title(modname+', RV = 2.7')
    plt.tight_layout()
    plt.savefig('outputs/models/'+modname+'_2.7_.png',facecolor='white')
    plt.close()
    

    plt.figure(3,figsize=(8,6),dpi=250)
    divnorm = colors.TwoSlopeNorm(vmin=-0.5, vcenter=0.0, vmax=0.5)
    mycmap = plt.get_cmap('coolwarm')(np.linspace(0.0,1.0,11))
    avs = np.round(np.linspace(-0.5,0.5,11),2)
    avs[2] = 0.0
    colst = []
    colnm = ['UVW2-V','UVM2-V','UVW1-V','U-V','B-V']
    plt.subplot(212)
    for i in range(10):
        nf = reddening(avs[i],3.1)
        for k in range(len(nf)):
            if nf[k] < 0.0:
                nf[k] = 0.0
        plt.plot(SlambdaWave, nf,label=str(avs[i]), color=mycmap[i])
        nucol = kfun(SlambdaWave,nf)
        colst.append(nucol)
        redata = {'wavelength':SlambdaWave,'flux':nf}
        nux = pd.DataFrame(data=redata)
        nux.to_csv('outputs/models/'+modname+'_3.1_'+str(avs[i])+'_.csv',index=False)
    plt.xlabel('Wavelength ($\AA$)')
    plt.ylabel('Flux (Arbitrary Units)')
    #plt.legend(title='AV =',loc='upper left',ncol=2)
    plt.xlim(1500,6000)
    plt.subplot(211)
    for i in range(len(colst)):
        plt.scatter(colnm,colst[i],color=mycmap[i],label=str(avs[i]))
    plt.ylabel('Color')
    plt.legend(title='AV =',loc = 'upper right',ncol=2)
    plt.title(modname+', RV = 3.1')
    plt.tight_layout()
    plt.savefig('outputs/models/'+modname+'_3.1_.png',facecolor='white')
    plt.close()
    

for i in range(len(dms)):
    fgen(dms[i])
