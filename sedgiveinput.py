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
from astropy.coordinates import SkyCoord
import sedinput_checker
from sedinput_checker import getAVbest
import pandas as pd





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


#targ_name='SN2011by'

targetlist  = pd.read_csv('SNPY_Sample_Decline.csv')

for i in range(len(targetlist)):
	targ_name=targetlist.sname[i]


	radec = pd.read_csv('NewSwiftSNweblist.csv')
	rd = radec.loc[(radec.SNname == targ_name)]
	rd = rd.reset_index(drop=True)
	ratp = rd.SNra[0].split("'")[1]
	dectp = rd.SNdec[0]
	c = SkyCoord(ra=ratp, dec=dectp, frame='icrs',unit=(u.hourangle, u.deg))

	mwreddening,_ = get_dust_RADEC(c.ra.degree, c.dec.degree, calibration='SF11')

	mwreddening = mwreddening[0]
	cfk5 = c.transform_to('fk5')
	mwav = getAVbest(cfk5)
	mwebv=mwav[0]/3.1
	#print(mwreddening, mwebv)

	# magnitudes and redshift
	df = pd.read_csv('SNPY_Sample_Decline.csv')
	tp = df.loc[(df.sname == targ_name)]
	tp = tp.reset_index(drop=True)
	z = tp.z[0]
	dm15 = tp.dm[0]
	obs_col = [tp.w2mag[0], tp.m2mag[0], tp.w1mag[0], tp.umag[0], tp.bmag[0]] - tp.vmag[0]


	#run_res = colorcomp(SlambdaWave,SpectrumFlux,tp.z[0],obs_col)

	#sed_checker.sed_checker(targetlist.sname[i])

	print(targ_name)
	sedinput_checker.sedinput_checker(targ_name,FilterVec,obs_col, dm15, z, mwebv)