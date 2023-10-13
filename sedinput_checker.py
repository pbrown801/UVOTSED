from pyphot import Filter
from pyphot import unit
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import extinction
import time
import sys
import os
import math
import matplotlib.colors as colors
from snpy.utils.IRSA_dust_getval import get_dust_RADEC, get_dust_sigma_RADEC
import glob
from astropy import units as u
from astropy.coordinates import SkyCoord,Angle
from astroquery.irsa_dust import IrsaDust


###########
#Set up Swift UVOT filters 

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

SWIFT_UVOT_B = Filter(Bdata.Wavelength,Bdata.Area,name = 'SWIFT_UVOT_B',dtype = 'photon',unit = 'Angstrom')
SWIFT_UVOT_U = Filter(Udata.Wavelength,Udata.Area,name = 'SWIFT_UVOT_U',dtype = 'photon',unit = 'Angstrom')
SWIFT_UVOT_V = Filter(Vdata.Wavelength,Vdata.Area,name = 'SWIFT_UVOT_V',dtype = 'photon',unit = 'Angstrom')
SWIFT_UVOT_UVM2 = Filter(M2data.Wavelength,M2data.Area,name = 'SWIFT_UVOT_UVM2',dtype = 'photon',unit = 'Angstrom')
SWIFT_UVOT_UVW1 = Filter(W1data.Wavelength,W1data.Area,name = 'SWIFT_UVOT_UVW1',dtype = 'photon',unit = 'Angstrom')
SWIFT_UVOT_UVW2 = Filter(W2data.Wavelength,W2data.Area,name = 'SWIFT_UVOT_UVW2',dtype = 'photon',unit = 'Angstrom')
FilterVec = [SWIFT_UVOT_B,SWIFT_UVOT_U,SWIFT_UVOT_V,SWIFT_UVOT_UVM2,SWIFT_UVOT_UVW1,SWIFT_UVOT_UVW2]


#this functinon calculate the filter magnitudes for a given spectrum
def uvotmags(wave,flux,FilterVec):
    print('uvotmags')
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
def uvotcolors(wave,flux,FilterVec):
    print('uvotcolors')
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

def colorcomp(wave,flux,z,obscol):
    print('colorcomp')
    SpectrumFlux = np.asarray(flux)#/max(flux)
    SlambdaWave = np.asarray(wave)*(1+z)
    modcol = uvotcolors(SlambdaWave,SpectrumFlux,FilterVec)
    difcol = obscol-modcol
    rmse = np.nansum(difcol**2 / 5)
    return rmse

#calculate the corrections

def mwcor(wave,flux,ebv,obv,dm,z):
    print('mwcor')
    #MW corrections
    rv = 3.1
    av = rv*ebv
    wave = wave*(1+z)
    mwf = extinction.apply(extinction.ccm89(wave,-1.0*av,rv),flux)
    col_obs = uvotmags(wave,flux,FilterVec)
    col_mw = uvotmags(wave,mwf,FilterVec)
    
    #K corrections
    nwave = wave/(1+z)
    kcf = mwf*(1+z)
    col_z = uvotmags(nwave,kcf,FilterVec)
    
    #Host Corrections
    #the conversion from intrinsic color from decline from phillips 99
    BV0 = 0.114*(dm-1.1) - 0.07
    #obv is the observed b-v, ebv is the MW ebv, BV0 is the intrinsic color from p99, HBV is thus the host excess
    HBV = obv[4] - ebv - BV0
    rv = 2.6
    av = rv*HBV
    #print(BV0,HBV,rv,av)
    hof = extinction.apply(extinction.ccm89(nwave,-1.0*av,rv),kcf)
    
    col_ho = uvotmags(nwave,hof,FilterVec)
    
    cors = np.array([col_obs-col_mw,col_mw-col_z,col_z-col_ho,(col_obs-col_mw)+(col_mw-col_z)+(col_z-col_ho)])
    return cors


#This function actually calculates the RMSE for each model and saves the outputs
def runsed1(modat,flflg,z,obs_col):
    print('runsed1')
    if flflg == True:
        premod = ['outputs/models/'+modat]
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
        run_res = colorcomp(SlambdaWave,SpectrumFlux,z,obs_col)
        res_run.append(run_res)
        tp1 = premod[i].split('/')[-1]
        tempname.append(tp1)
    
    return res_run, tempname

#This function finally calculates all of the corrections for a given model
def runsed2(tempname,res_run,z,targ_name,obs_col,dm15,mwebv):
    print('runsed2')

    thiscor= pd.DataFrame(columns = ['Template','RMSE','W2_MW','M2_MW','W1_MW','U_MW','B_MW','V_MW','W2_Z','M2_Z','W1_Z','U_Z','B_Z','V_Z',
    'W2_HOST','M2_HOST','W1_HOST','U_HOST','B_HOST','V_HOST','W2_TOT','M2_TOT','W1_TOT','U_TOT','B_TOT','V_TOT'])

    Tdata = pd.read_csv('outputs/models/'+tempname)
    Tdata.columns = ['wavelength','flux']
    SpectrumFlux = np.asarray(Tdata.flux)/max(Tdata.flux)
    SlambdaWave = np.asarray(Tdata.wavelength)
    corec = mwcor(SlambdaWave,SpectrumFlux,mwebv,obs_col,dm15,z)
    thiscor = thiscor.append({'Template':tempname, 'RMSE':res_run, 
        'W2_MW':corec[0][0], 'M2_MW':corec[0][1],'W1_MW':corec[0][2],'U_MW':corec[0][3],'B_MW':corec[0][4], 'V_MW':corec[0][5],
        'W2_Z':corec[1][0],'M2_Z':corec[1][1],'W1_Z':corec[1][2],'U_Z':corec[1][3],'B_Z':corec[1][4], 'V_Z':corec[1][5],
        'W2_HOST':corec[2][0],'M2_HOST':corec[2][1],'W1_HOST':corec[2][2],'U_HOST':corec[2][3],'B_HOST':corec[2][4], 'V_HOST':corec[2][5],
        'W2_TOT':corec[3][0],'M2_TOT':corec[3][1],'W1_TOT':corec[3][2],'U_TOT':corec[3][3],'B_TOT':corec[3][4], 'V_TOT':corec[3][5]}, ignore_index=True)


    return thiscor

def getAVbest(inputcoordinates):
    print('getAVbest')
    '''
    Coordinates are input as a single string. Output is the recommended Av value for MW reddening, error, and reference
    '''
    
    #inputcoordinates = sys.argv[1]
    testCoords = SkyCoord(inputcoordinates,frame='fk5')
    #print('\n-----\nReading input files...')
    inFile = 'Brown_Walker_table_1.dat'
    inTable = pd.read_csv(inFile,header=None,delimiter=' ')
    #ra = Angle(inTable.iloc[:,1])
    #dec = Angle(inTable.iloc[:,2])
    #sourceCoords = SkyCoord(ra,dec,frame='fk5')
    #ra = Angle(inTable.iloc[:,1])
    #dec = Angle(inTable.iloc[:,2])
    sourceCoords = SkyCoord(inTable.iloc[:,1],inTable.iloc[:,2],frame='fk5')
    
    #print('Calculating separation from table coordinates')
    separations = testCoords.separation(sourceCoords).arcminute
    # compare to the distances in the table
    within = np.less(separations,inTable.iloc[:,3])
    
    # Are any of the input coordinates within the tabulated distance 
    # of the coordinates in the table?
    correctedAV = np.where(within,inTable.iloc[:,4],None) #get calculated value
    fix=any(within)
    #print('fix?',fix)
    
    if fix:
        AV = next((item for item in correctedAV if item is not None),None)
        correctedAVerr = np.where(within,inTable.iloc[:,5],None) #get calculated val
        newAVerr = next((item for item in correctedAVerr if item is not None),None)
        AVerr = math.sqrt((int(float(newAVerr)))**2+(int(AV)*0.1)**2)
        sources=np.where(within,inTable.iloc[:,6],None)
        source = next((item for item in sources if item is not None),None)+",S_F_2011"
        
    if not fix:
        AVtable = IrsaDust.get_extinction_table(testCoords,show_progress = False)
        AV=AVtable['A_SandF'][2]
        AVerr = AV*0.1
        source = 'S_F_2011'

    #print(AV, AVerr, source)
    return (AV, AVerr, source)


 #This is the main program, set up as a function that takes a supernova name and the filter setup
def sedinput_checker(targ_name,FilterVec,obs_col, dm15, z, mwebv):
	print('sed_checker')

	#Initialize program with the inputs
	start = time.time()
	#data file of the SNe in the sample
	#df = pd.read_csv('SNPY_Sample_Decline.csv')
	#tp = df.loc[(df.sname == targ_name)]
	#tp = tp.reset_index(drop=True)
	#z = tp.z[0]
	#dm15 = tp.dm[0]
	#generate folder with for each supernova to hold the outputs
	#snm = tp.sname[0]
	#if path exists == true, output path = path. Else, create folder and set output path
	snm=targ_name
	if os.path.exists('outputs/'+snm+'/') == True:
		opath = 'outputs/'+snm+'/'
	else:
		os.mkdir('outputs/'+snm+'/')
		opath = 'outputs/'+snm+'/'

	# after this, all print statements go to this file
	sys.stdout = open(opath+'output.txt', 'w')
	print("starting at "+str(time.time()))


	#Setup up the spectral models based on the observed templates
	mods = np.asarray(['ASASSN-14LP_peak_11fe_appended.dat','SN1992A_UV.dat','SN2011by_peak_11fe_appended.dat',
		'SN2011fe_peak_11fe_appended.dat','SN2011iv_peak_11fe_appended.dat','SN2015F_peak_11fe_appended.dat',
		'SN2016ccj_peak_11fe_appended.dat','SN2017erp_peak_11fe_appended.dat','SN2021fxy_peak_11fe_appended.dat',
		'SN2022hrs_peak_11fe_appended.dat','SN2013dy_peak_11fe_appended.dat'])

	mod_dm = [0.8, 1.47, 1.14, 1.05, 1.77, 1.26, 0.67, 1.11, 1.05, 1.41, 0.92]

	print("1 "+str(time.time()))


	#Set up the functions that actually run everything all together
	fulcor = pd.DataFrame(columns = ['Template','RMSE','W2_MW','M2_MW','W1_MW','U_MW','B_MW','V_MW','W2_Z','M2_Z','W1_Z','U_Z','B_Z','V_Z',
		'W2_HOST','M2_HOST','W1_HOST','U_HOST','B_HOST','V_HOST','W2_TOT','M2_TOT','W1_TOT','U_TOT','B_TOT','V_TOT'])

	#obs_col = [tp.w2mag[0], tp.m2mag[0], tp.w1mag[0], tp.umag[0], tp.bmag[0]] - tp.vmag[0]


	#now that all the stuff is set up, this loop pulls any observed template within 0.2 dm15 of the observed SN and calculates RMSE
	miz1, miz2= [],[]
	for r in range(len(mod_dm)):
		if (dm15-0.2) <= mod_dm[r] and mod_dm[r] <= (dm15+0.2):
			rezzy, reznm = runsed1(mods[r],False,z,obs_col)
			miz1 = miz1 + rezzy
			miz2 = miz2 + reznm
	print("2 "+str(time.time()))

	#This does the same as above, but for the foley models
	fols = glob.glob('outputs/models/foley*.csv')
	for r in range(len(fols)):
		fn1 = fols[r].split('/')[-1]
		fn2 = fn1.split('_')[0]
		fn3 = fn2.split('y')[1]
		if (dm15-0.2) <= float(fn3) and float(fn3) <= (dm15+0.2):
			rezzy, reznm = runsed1(fn1,True,z,obs_col)
			miz1 = miz1 + rezzy
			miz2 = miz2 + reznm
	print("3 "+str(time.time()))

	#Once all the models are run, this pulls the top 10 tempaltes and runs the correction functions
	miz3 = sorted(miz1)

	for i in range(len(miz1)):
		if miz1[i] <= miz3[9]:
			driz = runsed2(miz2[i],miz1[i],z,targ_name,obs_col,dm15,mwebv)
			fulcor = fulcor.append(driz)

	end = time.time()
	print("Finishing at "+str(time.time()))
	print('Time to completion is '+str(end-start))

	#saves the output
	fulcor.to_csv(opath+'output.csv',index=False)
	sys.stdout.close()
	return


if __name__ == '__main__':
    args = sys.argv
    globals()[args[1]](args[2])(args[3])(args[4])(args[5])(args[6])
    print(args)
