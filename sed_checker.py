from pyphot import Filter
from pyphot import unit
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import extinction
import time
import sys
import os
import matplotlib.colors as colors
from snpy.utils.IRSA_dust_getval import get_dust_RADEC, get_dust_sigma_RADEC
import glob
from astropy import units as u
from astropy.coordinates import SkyCoord


 #This is the main program, set up as a function that takes a supernova name as the only input.
def sed_checker(targ_name):

	#Initialize program with the inputs
	start = time.time()
	#data file of the SNe in the sample
	df = pd.read_csv('SNPY_Sample_Decline.csv')
	tp = df.loc[(df.sname == targ_name)]
	tp = tp.reset_index(drop=True)


	#generate folder with for each supernova to hold the outputs
	snm = tp.sname[0]
	#if path exists == true, output path = path. Else, create folder and set output path
	if os.path.exists('outputs/'+snm+'/') == True:
		opath = 'outputs/'+snm+'/'
	else:
		os.mkdir('outputs/'+snm+'/')
		opath = 'outputs/'+snm+'/'

	sys.stdout = open(opath+'output.txt', 'w')

	#Set up Swift UVOT filters in order to generate SEDS

	Udata = pd.read_csv('UVabsmags-master/filters/U_P08.txt',delim_whitespace=True,comment='#')
	Udata.columns = ['Wavelength','Area']

	Bdata = pd.read_csv('UVabsmags-master/filters/B_P08.txt',delim_whitespace=True,comment='#')
	Bdata.columns = ['Wavelength','Area']

	Vdata = pd.read_csv('UVabsmags-master/filters/V_P08.txt',delim_whitespace=True,comment='#')
	Vdata.columns = ['Wavelength','Area']


	W2data = pd.read_csv('UVabsmags-master/filters/UVW2_B11.txt',delim_whitespace=True,comment='#')
	W2data.columns = ['Wavelength','Area']

	M2data = pd.read_csv('UVabsmags-master/filters/UVM2_B11.txt',delim_whitespace=True,comment='#')
	M2data.columns = ['Wavelength','Area']

	W1data = pd.read_csv('UVabsmags-master/filters/UVW1_B11.txt',delim_whitespace=True,comment='#')
	W1data.columns = ['Wavelength','Area']


	SWIFT_UVOT_B = Filter(Bdata.Wavelength,Bdata.Area,name = 'SWIFT_UVOT_B',dtype = 'photon',unit = 'Angstrom')
	SWIFT_UVOT_U = Filter(Udata.Wavelength,Udata.Area,name = 'SWIFT_UVOT_U',dtype = 'photon',unit = 'Angstrom')
	SWIFT_UVOT_V = Filter(Vdata.Wavelength,Vdata.Area,name = 'SWIFT_UVOT_V',dtype = 'photon',unit = 'Angstrom')
	SWIFT_UVOT_UVM2 = Filter(M2data.Wavelength,M2data.Area,name = 'SWIFT_UVOT_UVM2',dtype = 'photon',unit = 'Angstrom')
	SWIFT_UVOT_UVW1 = Filter(W1data.Wavelength,W1data.Area,name = 'SWIFT_UVOT_UVW1',dtype = 'photon',unit = 'Angstrom')
	SWIFT_UVOT_UVW2 = Filter(W2data.Wavelength,W2data.Area,name = 'SWIFT_UVOT_UVW2',dtype = 'photon',unit = 'Angstrom')
	FilterVec = [SWIFT_UVOT_B,SWIFT_UVOT_U,SWIFT_UVOT_V,SWIFT_UVOT_UVM2,SWIFT_UVOT_UVW1,SWIFT_UVOT_UVW2]

	#Setup up the spectral models based on the observed templates
	mods = np.asarray(['ASASSN-14LP_peak_11fe_appended.dat','SN1992A_UV.dat','SN2011by_peak_11fe_appended.dat',
		'SN2011fe_peak_11fe_appended.dat','SN2011iv_peak_11fe_appended.dat','SN2015F_peak_11fe_appended.dat',
		'SN2016ccj_peak_11fe_appended.dat','SN2017erp_peak_11fe_appended.dat','SN2021fxy_peak_11fe_appended.dat',
		'SN2022hrs_peak_11fe_appended.dat','SN2013dy_peak_11fe_appended.dat'])

	mod_dm = [0.8, 1.47, 1.14, 1.05, 1.77, 1.26, 0.67, 1.11, 1.05, 1.41, 0.92]


	#this functinon calculate the filter magnitudes for a given spectrum
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
	    
	    return np.array([f_w2,f_m2,f_w1,f_u,f_b,f_v])

	#calculate the X-V colors from a given spectrum
	def kfun1(wave,flux):
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
	    
	    return [f_w2,f_m2,f_w1,f_u,f_b]-f_v

	#Call in the models, redshift them, calculate color difference, and store the RMSE

	def colorcomp(wl,flux,z,obscol): 
	    
	    SpectrumFlux = np.asarray(flux)#/max(flux)
	    SlambdaWave = np.asarray(wl)*(1+z)
	    modcol = kfun1(SlambdaWave,SpectrumFlux)
	    difcol = obscol-modcol
	    rmse = np.nansum(difcol**2 / 5)
	    return rmse

	#calculate the corrections 

	def mwcor(wl,flux,ebv,obv,dm,z):
		#MW corrections
	    rv = 3.1
	    av = rv*ebv
	    wl = wl*(1+z)
	    mwf = extinction.apply(extinction.ccm89(wl,-1.0*av,rv),flux)
	    col_obs = kfun(wl,flux)
	    col_mw = kfun(wl,mwf)
	    
	    #K corrections
	    nwl = wl/(1+z)
	    kcf = mwf*(1+z)
	    col_z = kfun(nwl,kcf)
	    
	    #Host Corrections
	    #the conversion from intrinsic color from decline from phillips 99
	    BV0 = 0.114*(dm-1.1) - 0.07
	    #obv is the observed b-v, ebv is the MW ebv, BV0 is the intrinsic color from p99, HBV is thus the host excess
	    HBV = obv[4] - ebv - BV0
	    rv = 2.6
	    av = rv*HBV
	    print(BV0,HBV,rv,av)
	    hof = extinction.apply(extinction.ccm89(nwl,-1.0*av,rv),kcf)

	    col_ho = kfun(nwl,hof)
	    
	    cors = np.array([col_obs-col_mw,col_mw-col_z,col_z-col_ho,(col_obs-col_mw)+(col_mw-col_z)+(col_z-col_ho)])
	    return cors


	#Set up the functions that actually run everything all together
	fulcor = pd.DataFrame(columns = ['Template','RMSE','W2_MW','M2_MW','W1_MW','U_MW','B_MW','V_MW','W2_Z','M2_Z','W1_Z','U_Z','B_Z','V_Z',
		'W2_HOST','M2_HOST','W1_HOST','U_HOST','B_HOST','V_HOST','W2_TOT','M2_TOT','W1_TOT','U_TOT','B_TOT','V_TOT'])

	obs_col = [tp.w2mag[0], tp.m2mag[0], tp.w1mag[0], tp.umag[0], tp.bmag[0]] - tp.vmag[0]

	#This function actually calculates the RMSE for each model and saves the outputs
	def runsed1(modat,flflg):
		if flflg == True:
			premod = ['outputs/models/'+modat]
			print(premod)
		else:
			modname = modat.split('_')[0]
			premod = glob.glob('outputs/models/'+modname+'*.csv')
		res_run = []
		tempname = []
		dw1v, dbv = [],[]
		for i in range(len(premod)):
			Tdata = pd.read_csv(premod[i])
			Tdata.columns = ['wavelength','flux']
			print(premod[i])
			SpectrumFlux = np.asarray(Tdata.flux)/max(Tdata.flux)
			SlambdaWave = np.asarray(Tdata.wavelength)
			run_res = colorcomp(SlambdaWave,SpectrumFlux,tp.z[0],obs_col)
			res_run.append(run_res)
			tp1 = premod[i].split('/')[-1]
			tempname.append(tp1)
		
		return res_run, tempname

	#This function finally calculates all of the corrections for a given model
	def runsed2(tempname,res_run):
		radec = pd.read_csv('NewSwiftSNweblist.csv')
		rd = radec.loc[(radec.SNname == targ_name)]
		rd = rd.reset_index(drop=True)
		ratp = rd.SNra[0].split("'")[1]
		dectp = rd.SNdec[0]
		c = SkyCoord(ra=ratp, dec=dectp, frame='icrs',unit=(u.hourangle, u.deg))

		mwreddening,_ = get_dust_RADEC(c.ra.degree, c.dec.degree, calibration='SF11')

		mwreddening = mwreddening[0]

		thiscor= pd.DataFrame(columns = ['Template','RMSE','W2_MW','M2_MW','W1_MW','U_MW','B_MW','V_MW','W2_Z','M2_Z','W1_Z','U_Z','B_Z','V_Z',
		'W2_HOST','M2_HOST','W1_HOST','U_HOST','B_HOST','V_HOST','W2_TOT','M2_TOT','W1_TOT','U_TOT','B_TOT','V_TOT'])

		Tdata = pd.read_csv('outputs/models/'+tempname)
		Tdata.columns = ['wavelength','flux']
		SpectrumFlux = np.asarray(Tdata.flux)/max(Tdata.flux)
		SlambdaWave = np.asarray(Tdata.wavelength)
		corec = mwcor(SlambdaWave,SpectrumFlux,mwreddening,obs_col,tp.dm[0],tp.z[0])
		thiscor = thiscor.append({'Template':tempname, 'RMSE':res_run, 
			'W2_MW':corec[0][0], 'M2_MW':corec[0][1],'W1_MW':corec[0][2],'U_MW':corec[0][3],'B_MW':corec[0][4], 'V_MW':corec[0][5],
			'W2_Z':corec[1][0],'M2_Z':corec[1][1],'W1_Z':corec[1][2],'U_Z':corec[1][3],'B_Z':corec[1][4], 'V_Z':corec[1][5],
			'W2_HOST':corec[2][0],'M2_HOST':corec[2][1],'W1_HOST':corec[2][2],'U_HOST':corec[2][3],'B_HOST':corec[2][4], 'V_HOST':corec[2][5],
			'W2_TOT':corec[3][0],'M2_TOT':corec[3][1],'W1_TOT':corec[3][2],'U_TOT':corec[3][3],'B_TOT':corec[3][4], 'V_TOT':corec[3][5]}, ignore_index=True)


		return thiscor

	#now that all the stuff is set up, this loop pulls any observed template within 0.2 dm15 of the observed SN and calculates RMSE
	miz1, miz2= [],[]
	for r in range(len(mod_dm)):
		if (tp.dm[0]-0.2) <= mod_dm[r] and mod_dm[r] <= (tp.dm[0]+0.2):
			rezzy, reznm = runsed1(mods[r],False)
			miz1 = miz1 + rezzy
			miz2 = miz2 + reznm

	#This does the same as above, but for the foley models
	fols = glob.glob('outputs/models/foley*.csv')
	for r in range(len(fols)):
		fn1 = fols[r].split('/')[-1]
		fn2 = fn1.split('_')[0]
		fn3 = fn2.split('y')[1]
		if (tp.dm[0]-0.2) <= float(fn3) and float(fn3) <= (tp.dm[0]+0.2):
			rezzy, reznm = runsed1(fn1,True)
			miz1 = miz1 + rezzy
			miz2 = miz2 + reznm

	#Once all the models are run, this pulls the top 10 tempaltes and runs the correction functions
	miz3 = sorted(miz1)

	for i in range(len(miz1)):
		if miz1[i] <= miz3[9]:
			driz = runsed2(miz2[i],miz1[i])

	end = time.time()
	print('Time to completion is '+str(end-start))

	#saves the output
	fulcor.to_csv(opath+'output.csv',index=False)
	return


if __name__ == '__main__':
    args = sys.argv
    globals()[args[1]](args[2])
    print(args)