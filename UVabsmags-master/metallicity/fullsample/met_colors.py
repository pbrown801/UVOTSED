import math
import matplotlib.pyplot as plt
import numpy as np
import linmix
import os
import astropy
from astropy.io import fits
from astropy.table import Table
from scipy import stats
from sklearn.utils import resample
import itertools


def Import_Data(path):
	data=np.genfromtxt(path, delimiter='&', dtype=str, comments='#', missing_values='-99.90', unpack=True)
	return(data)


########################  SN metallicity
# read in a file from Graham et al. 2019
# https://ui.adsabs.harvard.edu/abs/2019arXiv190513197G/abstract

metpath='SNe_Z_table_dumb.txt'

metSNname_array, metredshifts_array, metMB_array, metKK04_array, metT04_array, metD16_array, mettype_array=Import_Data(metpath)

####### remove spaces in name, i.e. SN 2005am --> SN2005am
nospaces=metSNname_array[:]
nospaces = [item.replace(" ", "") for item in metSNname_array]
metSNname_array=nospaces[:]

metKK04_array=[float(i) for i in metKK04_array]
metT04_array =[float(i) for i in metT04_array]
metD16_array =[float(i) for i in metD16_array]



#############  read in the massive csv file of data read in from the spreadsheet 
############# and computed in IDL

hostpath='../../host.sav.190819.csv'

host_SNname_array, host_hostname_array, host_PSNname_array, host_SNname2_array, host_redshift_array, host_redshift_err_array, host_redshifttext_array, host_SNtype_array, host_SNtype2_array, host_Nobs_array, host_TargetIDs_array, host_SNra_array, host_SNdec_array, host_template_status_array, host_max_exptime_array, host_filters_array, host_notes_array, host_sibling_SNe_array, host_RAhost_array, host_DEChost_array, host_gal_lat_host_array, host_gal_long_host_array, host_av_schlafly_array, host_redshift_ref_array, host_velocity_array, host_velocity_err_array,  host_corrvelocity_array,  host_corrvelocity_err_array,  host_dm_cepheid_array, host_dm_sbf_array, host_dm_pnlf_array, host_dm_other_array, host_dm_tf_array, host_dm_cepheid_err_array, host_dm_sbf_err_array, host_dm_pnlf_err_array, host_dm_other_err_array, host_dm_tf_err_array, host_dm_cepheid_ref_array, host_dm_sbf_ref_array, host_dm_pnlf_ref_array, host_dm_other_ref_array, host_dm_tf_ref_array, host_dm_z_array, host_dm_z_err_array, host_dm_z_test_array, host_dm_z_test_err_array, host_dm_z_cor_array, host_dm_z_cor_err_array, host_dm_best_array, host_dm_best_err_array, host_dm_best_ref_array, host_Host_Morphology_array, host_morph_array, host_Hubbletype_array, host_petrosian_array, host_petrosianminor_array, host_C9050_w2_array, host_C9050_m2_array, host_C9050_w1_array, host_C9050_uu_array, host_C9050_bb_array, host_C9050_vv_array, host_galmag_w2_array, host_galmag_m2_array, host_galmag_w1_array, host_galmag_uu_array, host_galmag_bb_array, host_galmag_vv_array, host_galmag_w2_err_array, host_galmag_m2_err_array, host_galmag_w1_err_array, host_galmag_uu_err_array, host_galmag_bb_err_array, host_galmag_vv_err_array, host_appmag_w2_array, host_appmag_m2_array, host_appmag_w1_array, host_appmag_array, host_appmag_bb_array, host_appmag_vv_array, host_appmagerr_w2_array, host_appmagerr_m2_array, host_appmagerr_w1_array, host_appmagerr_uu_array, host_appmagerr_bb_array, host_appmagerr_array, host_bpeakappmag_w2_array, host_bpeakappmag_m2_array, host_bpeakappmag_w1_array, host_bpeakappmag_uu_array, host_bpeakappmag_bb_array, host_bpeakappmag_vv_array, host_bpeakmjd_array, host_bpeakappmagerr_w2_array, host_bpeakappmagerr_m2_array, host_bpeakappmagerr_w1_array, host_bpeakappmagerr_uu_array, host_bpeakappmagerr_bb_array, host_bpeakappmagerr_vv_array, host_peaktime_w2_array, host_peaktime_m2_array, host_peaktime_w1_array, host_peaktime_uu_array, host_peaktime_bb_array, host_peaktime_vv_array, host_firstepochmjd_w2_array, host_firstepochmjd_m2_array, host_firstepochmjd_w1_array, host_firstepochmjd_uu_array, host_firstepochmjd_bb_array, host_firstepochmjd_vv_array, host_firstepochmjd_array, host_absmag_w2_array, host_absmag_m2_array, host_absmag_w1_array, host_absmag_uu_array, host_absmag_bb_array, host_absmag_vv_array, host_absmagerr_w2_array, host_absmagerr_m2_array, host_absmagerr_w1_array, host_absmagerr_uu_array, host_absmagerr_bb_array, host_absmagerr_vv_array, host_bestbabsmag_w2_array, host_bestbabsmag_m2_array, host_bestbabsmag_w1_array, host_bestbabsmag_uu_array, host_bestbabsmag_bb_array, host_bestbabsmag_vv_array, host_bestbabsmagerr_w2_array, host_bestbabsmagerr_m2_array, host_bestbabsmagerr_w1_array, host_bestbabsmagerr_uu_array, host_bestbabsmagerr_bb_array, host_bestbabsmagerr_vv_array, host_bestvabsmag_w2_array, host_bestvabsmag_m2_array, host_bestvabsmag_w1_array, host_bestvabsmag_uu_array, host_bestvabsmag_bb_array, host_bestvabsmag_vv_array, host_bestvabsmagerr_w2_array, host_bestvabsmagerr_m2_array, host_bestvabsmagerr_w1_array, host_bestvabsmagerr_uu_array, host_bestvabsmagerr_bb_array, host_bestvabsmagerr_vv_array, host_dm15_w2_array, host_dm15_m2_array, host_dm15_w1_array, host_dm15_uu_array, host_dm15_bb_array, host_dm15_vv_array, host_dm15err_w2_array, host_dm15err_m2_array, host_dm15err_w1_array, host_dm15err_uu_array, host_dm15err_bb_array, host_dm15err_vv_array, host_v_best_err_array, host_v_best_array, host_b_best_err_array, host_b_best_array, host_dmb15_best_err_array, host_dmb15_best_array, host_velocityclass_array, host_vgclass_array, host_si0velocity_array, host_ebvvelocityclass_array, host_ebvvelocity_array, host_data_array, host_mast_array, host_fit_array, host_galfile_array, host_galphot_array, host_png_array, end_array =np.genfromtxt(hostpath, delimiter=',', dtype=str, comments='#', missing_values=math.nan, unpack=True)


#print(host_bpeakappmag_w2_array)
# host_bpeakappmag_w2_array=[w.replace('nan', ' ') for w in host_bpeakappmag_w2_array]
# print()
#print(host_bpeakappmag_w2_array)


host_bpeakappmag_w2_array=[float(i) for i in host_bpeakappmag_w2_array]


host_bpeakappmag_m2_array=[float(i) for i in host_bpeakappmag_m2_array]
host_bpeakappmag_w1_array=[float(i) for i in host_bpeakappmag_w1_array]
host_bpeakappmag_uu_array=[float(i) for i in host_bpeakappmag_uu_array]
host_bpeakappmag_bb_array=[float(i) for i in host_bpeakappmag_bb_array]
host_bpeakappmag_vv_array=[float(i) for i in host_bpeakappmag_vv_array]




########## List of Swift Supernovae Ia
fullsnlistpath='../../snialist.txt'
fullsnlist, fullsnquality=np.genfromtxt(fullsnlistpath, delimiter=' ', dtype=str, comments='#', missing_values='0.0', unpack=True)

print('full sn list = ', len(fullsnlist))


goodsnlist=fullsnlist.compress((fullsnquality =='good').flat)
print('good sn list = ', len(goodsnlist))


#exit()

#############

#############   Now make arrays the same length as good SN list

#good_m2_bpeakmags_array=np.array([0.0]*len(goodsnlist))
good_w2_bpeakmags_err_array=np.array([math.nan]*len(goodsnlist))
good_m2_bpeakmags_err_array=np.array([math.nan]*len(goodsnlist))
good_w1_bpeakmags_err_array=np.array([math.nan]*len(goodsnlist))
good_uu_bpeakmags_err_array=np.array([math.nan]*len(goodsnlist))
good_bb_bpeakmags_err_array=np.array([math.nan]*len(goodsnlist))
good_vv_bpeakmags_err_array=np.array([math.nan]*len(goodsnlist))
good_w2_bpeakmags_array=np.array([math.nan]*len(goodsnlist))
good_m2_bpeakmags_array=np.array([math.nan]*len(goodsnlist))
good_w1_bpeakmags_array=np.array([math.nan]*len(goodsnlist))
good_uu_bpeakmags_array=np.array([math.nan]*len(goodsnlist))
good_bb_bpeakmags_array=np.array([math.nan]*len(goodsnlist))
good_vv_bpeakmags_array=np.array([math.nan]*len(goodsnlist))
good_KK04_array =np.array([math.nan]*len(goodsnlist))
good_T04_array =np.array([math.nan]*len(goodsnlist))
good_D16_array =np.array([math.nan]*len(goodsnlist))
#good_D16_array =np.array([0.0]*len(goodsnlist))


for j in range(len(goodsnlist)):
	SNname=goodsnlist[j]
	#print(' ')
	#print(SNname)
	metindex=np.where(np.array(metSNname_array) == SNname)[0]
	#print('met index ', metindex)
	if metindex.size > 0 and 0 <= metindex < len(metSNname_array):

		good_T04_array[j]=metT04_array[np.asscalar(metindex)]
		if good_T04_array[j] == 0.0 or good_T04_array[j] == -99.9:
			good_T04_array[j]=math.nan
		
		good_D16_array[j]=metD16_array[np.asscalar(metindex)]
		if good_D16_array[j] == 0.0 or good_D16_array[j] == -99.9:
			good_D16_array[j]=math.nan

		good_KK04_array[j]=metKK04_array[np.asscalar(metindex)]
		if good_KK04_array[j] == 0.0 or good_KK04_array[j] == -99.9:
			good_KK04_array[j]=math.nan

	hostindex=np.where(np.array(host_SNname_array) == SNname)[0]
	#print('host index', hostindex)
	
	if hostindex.size > 0 and 0 <= hostindex < len(host_SNname_array) and float(host_dm15_bb_array[np.asscalar(hostindex)]) < 1.6 and float(host_dm15_bb_array[np.asscalar(hostindex)]) > 1.0 :
				
		good_w2_bpeakmags_array[j]=host_bpeakappmag_w2_array[np.asscalar(hostindex)]
		good_m2_bpeakmags_array[j]=host_bpeakappmag_m2_array[np.asscalar(hostindex)]
		good_w1_bpeakmags_array[j]=host_bpeakappmag_w1_array[np.asscalar(hostindex)]
		good_uu_bpeakmags_array[j]=host_bpeakappmag_uu_array[np.asscalar(hostindex)]
		good_bb_bpeakmags_array[j]=host_bpeakappmag_bb_array[np.asscalar(hostindex)]
		good_vv_bpeakmags_array[j]=host_bpeakappmag_vv_array[np.asscalar(hostindex)]
				
		good_w2_bpeakmags_err_array[j]=host_bpeakappmagerr_w2_array[np.asscalar(hostindex)]
		good_m2_bpeakmags_err_array[j]=host_bpeakappmagerr_m2_array[np.asscalar(hostindex)]
		good_w1_bpeakmags_err_array[j]=host_bpeakappmagerr_w1_array[np.asscalar(hostindex)]
		good_uu_bpeakmags_err_array[j]=host_bpeakappmagerr_uu_array[np.asscalar(hostindex)]
		good_bb_bpeakmags_err_array[j]=host_bpeakappmagerr_bb_array[np.asscalar(hostindex)]
		good_vv_bpeakmags_err_array[j]=host_bpeakappmagerr_vv_array[np.asscalar(hostindex)]

good_EBV_bpeakmags_array=good_bb_bpeakmags_array-good_vv_bpeakmags_array

R_MW=[6.2,8.01,5.43,4.92,4.16,3.16]


good_w2_bpeakmags_excor_array=R_MW[0]*(-1)*good_EBV_bpeakmags_array+good_w2_bpeakmags_array
good_m2_bpeakmags_excor_array=R_MW[1]*(-1)*good_EBV_bpeakmags_array+good_m2_bpeakmags_array
good_w1_bpeakmags_excor_array=R_MW[2]*(-1)*good_EBV_bpeakmags_array+good_w1_bpeakmags_array
good_uu_bpeakmags_excor_array=R_MW[3]*(-1)*good_EBV_bpeakmags_array+good_uu_bpeakmags_array
good_bb_bpeakmags_excor_array=R_MW[4]*(-1)*good_EBV_bpeakmags_array+good_bb_bpeakmags_array
good_vv_bpeakmags_excor_array=R_MW[5]*(-1)*good_EBV_bpeakmags_array+good_vv_bpeakmags_array

#print(good_m2_bpeakmags_excor_array-good_bb_bpeakmags_excor_array)
#print(good_w1_bpeakmags_excor_array-good_bb_bpeakmags_excor_array)
#print(good_uu_bpeakmags_excor_array-good_bb_bpeakmags_excor_array)


#print(good_EBV_bpeakmags_array)
		
#print((good_KK04_array))
#print((good_T04_array))
#print((good_D16_array))
#print((good_w2_bpeakmags_array))

print('T04 ', np.count_nonzero(~np.isnan(good_T04_array)))
print('KK04 ', np.count_nonzero(~np.isnan(good_KK04_array)))
print('D16 ', np.count_nonzero(~np.isnan(good_D16_array)))




colors = (0,0,0)
area = np.pi*3

# Plot

#Set plotting parameters
#https://matplotlib.org/users/customizing.html
plt.rcParams['xtick.labelsize']=8
plt.rcParams['ytick.labelsize']=8
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.rcParams['figure.figsize']=[8.8,10]
plt.rcParams['axes.labelsize']=12
plt.rcParams['axes.titlesize']=14
plt.rcParams['legend.fontsize']=10
plt.rcParams['xtick.major.size']=6
plt.rcParams['xtick.minor.size']=3
plt.rcParams['xtick.minor.visible']=True
plt.rcParams['ytick.major.size']=6
plt.rcParams['ytick.minor.size']=3
plt.rcParams['ytick.minor.visible']=True
plt.rcParams['savefig.format']='png'
plt.rcParams['savefig.dpi']=250
plt.rcParams['figure.dpi']=250
plt.rcParams['lines.markersize']=2.0
plt.rcParams['lines.linewidth']=0.5
plt.rcParams['text.usetex']=False #command to process all strings using Latex
plt.rcParams['figure.subplot.hspace']=0
plt.rcParams['figure.subplot.wspace']=0
plt.rcParams['errorbar.capsize']=0.5



#plt.scatter(np.array(good_D16_array), good_w1_bpeakmags_excor_array-good_bb_bpeakmags_excor_array, s=area, alpha=0.5)
#plt.xlabel('Metallicity D16')
#plt.ylabel('uvw1 - b')
#plt.show()
#plt.savefig('w1-b_D16met.eps', format='eps')



ubrange=[0.2,-0.8]
w1brange=[1.8, 0.6]
m2brange=[5.8,3.0]
T04range=[8.0,9.4]
D16range=[8.4,9.4]

#plot
f, ((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,sharex='col',sharey='row')





ax1.yaxis.set_ticks_position('both')
ax1.xaxis.set_ticks_position('both')
ax1.set_ylim(ubrange)
ax1.set_xlim(D16range)
ax1.errorbar(np.array(good_D16_array), good_uu_bpeakmags_excor_array-good_bb_bpeakmags_excor_array, xerr=0.1, yerr=np.sqrt(np.power(good_w1_bpeakmags_err_array,2.0)+np.power(good_bb_bpeakmags_err_array,2.0)),ls='')
ax1.scatter(np.array(good_D16_array), good_uu_bpeakmags_excor_array-good_bb_bpeakmags_excor_array)
ax1.set_ylabel('u-b Color')
#filter out NAN values
x1=good_D16_array
x1sig=np.array([0.1]*len(good_D16_array))
y1=good_uu_bpeakmags_excor_array-good_bb_bpeakmags_excor_array
y1sig=np.sqrt(np.power(good_w1_bpeakmags_err_array,2.0)+np.power(good_bb_bpeakmags_err_array,2.0))
x=[]
xsig=[]
y=[]
ysig=[]
for i in range(len(x1)):
	if str(x1[i])!='nan':
		if str(y1[i])!='nan':
			if str(x1sig[i])!='nan':
				if str(y1sig[i])!='nan':
					x.append(x1[i])
					y.append(y1[i])
					xsig.append(x1sig[i])
					ysig.append(y1sig[i])
corr_coeff, p_value_pearson=stats.pearsonr(x,y)
tau, p_value_kendall=stats.kendalltau(x,y)
lm= linmix.LinMix(x, y, xsig, ysig, K=2)
lm.run_mcmc(silent=True)
for i in range(0, len(lm.chain), 25):
	xs = np.arange(6,13)
	ys = lm.chain[i]['alpha'] + xs * lm.chain[i]['beta']
	ax1.plot(xs, ys, color='r', alpha=0.02)
intercept_array=[]
slope_array=[]
for i in range(0, len(lm.chain), 25):
	intercept=lm.chain[i]['alpha']
	slope=lm.chain[i]['beta']
	intercept_array.append(intercept)
	slope_array.append(slope)
mean_intercept=np.mean(intercept_array)
std_intercept=np.std(intercept_array)
mean_slope=np.mean(slope_array)
std_slope=np.std(slope_array)
x_mean=np.arange(6,13)
y_mean=mean_intercept+x_mean*mean_slope
ax1.plot(x_mean, y_mean, color='k')
#Record data
filename='MCMC_Regression_For_Color_vs_D16.txt'
text_file=open(filename,'a+')
color='u-b'
slope=str(np.round(mean_slope,3))
slope_err=str(np.round(std_slope,3))
intercept=str(np.round(mean_intercept,3))
intercept_err=str(np.round(std_intercept,3))
K_cor_coe=str(np.round(tau,3))
K_pval=str(np.round(p_value_kendall,3))
if float(K_pval)==0:
	K_pval=str(p_value_kendall)
P_cor_coe=str(np.round(corr_coeff,3))
P_pval=str(np.round(p_value_pearson,3))
if float(P_pval)==0:
	P_pval=str(p_value_pearson)
text_file.write('\n'+color+','+slope+'$\pm$'+slope_err+','+intercept+'$\pm$'+intercept_err+','+K_cor_coe+','+K_pval+','+P_cor_coe+','+P_pval)
    
ax3.yaxis.set_ticks_position('both')
ax3.xaxis.set_ticks_position('both')
ax3.set_ylim(w1brange)
ax3.set_xlim(D16range)
ax3.errorbar(np.array(good_D16_array), good_w1_bpeakmags_excor_array-good_bb_bpeakmags_excor_array, xerr=0.1, yerr=np.sqrt(np.power(good_w1_bpeakmags_err_array,2.0)+np.power(good_bb_bpeakmags_err_array,2.0)),ls='')
ax3.scatter(np.array(good_D16_array), good_w1_bpeakmags_excor_array-good_bb_bpeakmags_excor_array)
ax3.set_ylabel('w1-b Color')
#filter out NAN values
x1=good_D16_array
x1sig=np.array([0.1]*len(good_D16_array))
y1=good_w1_bpeakmags_excor_array-good_bb_bpeakmags_excor_array
y1sig=np.sqrt(np.power(good_w1_bpeakmags_err_array,2.0)+np.power(good_bb_bpeakmags_err_array,2.0))


x=[]
xsig=[]
y=[]
ysig=[]
for i in range(len(x1)):
	if str(x1[i])!='nan':
		if str(y1[i])!='nan':
			if str(x1sig[i])!='nan':
				if str(y1sig[i])!='nan':
					x.append(x1[i])
					y.append(y1[i])
					xsig.append(x1sig[i])
					ysig.append(y1sig[i])
corr_coeff, p_value_pearson=stats.pearsonr(x,y)
tau, p_value_kendall=stats.kendalltau(x,y)
lm= linmix.LinMix(x, y, xsig, ysig, K=2)
lm.run_mcmc(silent=True)
for i in range(0, len(lm.chain), 25):
	xs = np.arange(6,13)
	ys = lm.chain[i]['alpha'] + xs * lm.chain[i]['beta']
	ax3.plot(xs, ys, color='r', alpha=0.02)
intercept_array=[]
slope_array=[]
for i in range(0, len(lm.chain), 25):
	intercept=lm.chain[i]['alpha']
	slope=lm.chain[i]['beta']
	intercept_array.append(intercept)
	slope_array.append(slope)
mean_intercept=np.mean(intercept_array)
std_intercept=np.std(intercept_array)
mean_slope=np.mean(slope_array)
std_slope=np.std(slope_array)
x_mean=np.arange(6,13)
y_mean=mean_intercept+x_mean*mean_slope
ax3.plot(x_mean, y_mean, color='k')
#Record data
color='w1-b'
slope=str(np.round(mean_slope,3))
slope_err=str(np.round(std_slope,3))
intercept=str(np.round(mean_intercept,3))
intercept_err=str(np.round(std_intercept,3))
K_cor_coe=str(np.round(tau,3))
K_pval=str(np.round(p_value_kendall,3))
if float(K_pval)==0:
	K_pval=str(p_value_kendall)
P_cor_coe=str(np.round(corr_coeff,3))
P_pval=str(np.round(p_value_pearson,3))
if float(P_pval)==0:
	P_pval=str(p_value_pearson)
text_file.write('\n'+color+','+slope+'$\pm$'+slope_err+','+intercept+'$\pm$'+intercept_err+','+K_cor_coe+','+K_pval+','+P_cor_coe+','+P_pval)
    
ax5.set_xlabel('log(OH) Metallicity D16')
ax5.yaxis.set_ticks_position('both')
ax5.xaxis.set_ticks_position('both')
ax5.set_ylim(m2brange)
ax5.set_xlim(D16range)
ax5.errorbar(np.array(good_D16_array), good_m2_bpeakmags_excor_array-good_bb_bpeakmags_excor_array, xerr=0.1, yerr=np.sqrt(np.power(good_m2_bpeakmags_err_array,2.0)+np.power(good_bb_bpeakmags_err_array,2.0)),ls='')
ax5.scatter(np.array(good_D16_array), good_m2_bpeakmags_excor_array-good_bb_bpeakmags_excor_array)
ax5.set_ylabel('m2-b Color')
ax5.set_xlabel('Host Metallicity D16')
#filter out NAN values
x1=good_D16_array
x1sig=np.array([0.1]*len(good_D16_array))
y1=good_m2_bpeakmags_excor_array-good_bb_bpeakmags_excor_array
y1sig=np.sqrt(np.power(good_m2_bpeakmags_err_array,2.0)+np.power(good_bb_bpeakmags_err_array,2.0))


x=[]
xsig=[]
y=[]
ysig=[]
for i in range(len(x1)):
	if str(x1[i])!='nan':
		if str(y1[i])!='nan':
			if str(x1sig[i])!='nan':
				if str(y1sig[i])!='nan':
					x.append(x1[i])
					y.append(y1[i])
					xsig.append(x1sig[i])
					ysig.append(y1sig[i])
corr_coeff, p_value_pearson=stats.pearsonr(x,y)
tau, p_value_kendall=stats.kendalltau(x,y)
lm= linmix.LinMix(x, y, xsig, ysig, K=2)
lm.run_mcmc(silent=True)
for i in range(0, len(lm.chain), 25):
	xs = np.arange(6,13)
	ys = lm.chain[i]['alpha'] + xs * lm.chain[i]['beta']
	ax5.plot(xs, ys, color='r', alpha=0.02)
intercept_array=[]
slope_array=[]
for i in range(0, len(lm.chain), 25):
	intercept=lm.chain[i]['alpha']
	slope=lm.chain[i]['beta']
	intercept_array.append(intercept)
	slope_array.append(slope)
mean_intercept=np.mean(intercept_array)
std_intercept=np.std(intercept_array)
mean_slope=np.mean(slope_array)
std_slope=np.std(slope_array)
x_mean=np.arange(6,13)
y_mean=mean_intercept+x_mean*mean_slope
ax5.plot(x_mean, y_mean, color='k')
#Record data
color='m2-b'
slope=str(np.round(mean_slope,3))
slope_err=str(np.round(std_slope,3))
intercept=str(np.round(mean_intercept,3))
intercept_err=str(np.round(std_intercept,3))
K_cor_coe=str(np.round(tau,3))
K_pval=str(np.round(p_value_kendall,3))
if float(K_pval)==0:
	K_pval=str(p_value_kendall)
P_cor_coe=str(np.round(corr_coeff,3))
P_pval=str(np.round(p_value_pearson,3))
if float(P_pval)==0:
	P_pval=str(p_value_pearson)
text_file.write('\n'+color+','+slope+'$\pm$'+slope_err+','+intercept+'$\pm$'+intercept_err+','+K_cor_coe+','+K_pval+','+P_cor_coe+','+P_pval)
    


ax2.yaxis.set_ticks_position('both')
ax2.xaxis.set_ticks_position('both')
ax2.set_ylim(ubrange)
ax2.set_xlim(T04range)

for q in range(len(good_T04_array)):
	if str(good_T04_array[q])!='nan': print(goodsnlist[q], good_T04_array[q],good_uu_bpeakmags_excor_array[q]-good_bb_bpeakmags_excor_array[q], good_T04_array[q],good_w1_bpeakmags_excor_array[q]-good_bb_bpeakmags_excor_array[q]) 
ax2.errorbar(np.array(good_T04_array), good_uu_bpeakmags_excor_array-good_bb_bpeakmags_excor_array, xerr=0.1, yerr=np.sqrt(np.power(good_uu_bpeakmags_err_array,2.0)+np.power(good_bb_bpeakmags_err_array,2.0)),ls='')
ax2.scatter(np.array(good_T04_array), good_uu_bpeakmags_excor_array-good_bb_bpeakmags_excor_array)
#ax2.set_ylabel('u-b Color')
#filter out NAN values
x1=good_T04_array
x1sig=np.array([0.1]*len(good_T04_array))
y1=good_uu_bpeakmags_excor_array-good_bb_bpeakmags_excor_array
y1sig=np.sqrt(np.power(good_uu_bpeakmags_err_array,2.0)+np.power(good_bb_bpeakmags_err_array,2.0))
x=[]
xsig=[]
y=[]
ysig=[]
for i in range(len(x1)):
	if str(x1[i])!='nan':
		if str(y1[i])!='nan':
			if str(x1sig[i])!='nan':
				if str(y1sig[i])!='nan':
					x.append(x1[i])
					y.append(y1[i])
					xsig.append(x1sig[i])
					ysig.append(y1sig[i])
corr_coeff, p_value_pearson=stats.pearsonr(x,y)
tau, p_value_kendall=stats.kendalltau(x,y)
lm= linmix.LinMix(x, y, xsig, ysig, K=2)
lm.run_mcmc(silent=True)
for i in range(0, len(lm.chain), 25):
	xs = np.arange(8,10)
	ys = lm.chain[i]['alpha'] + xs * lm.chain[i]['beta']
	if i==0:
		ax2.plot(xs, ys, color='r', alpha=0.3, label='MCMC Fit Lines')
	else:
		ax2.plot(xs, ys, color='r', alpha=0.02)
intercept_array=[]
slope_array=[]
for i in range(0, len(lm.chain), 25):
	intercept=lm.chain[i]['alpha']
	slope=lm.chain[i]['beta']
	intercept_array.append(intercept)
	slope_array.append(slope)
mean_intercept=np.mean(intercept_array)
std_intercept=np.std(intercept_array)
mean_slope=np.mean(slope_array)
std_slope=np.std(slope_array)
x_mean=np.arange(8,10)
y_mean=mean_intercept+x_mean*mean_slope
ax2.plot(x_mean, y_mean, color='k', label='average fit')
ax2.legend()
#Record data
filename='MCMC_Regression_For_Color_vs_T04.txt'
text_file=open(filename,'a+')
color='u-b'
slope=str(np.round(mean_slope,3))
slope_err=str(np.round(std_slope,3))
intercept=str(np.round(mean_intercept,3))
intercept_err=str(np.round(std_intercept,3))
K_cor_coe=str(np.round(tau,3))
K_pval=str(np.round(p_value_kendall,3))
if float(K_pval)==0:
	K_pval=str(p_value_kendall)
P_cor_coe=str(np.round(corr_coeff,3))
P_pval=str(np.round(p_value_pearson,3))
if float(P_pval)==0:
	P_pval=str(p_value_pearson)
text_file.write('\n'+color+','+slope+'$\pm$'+slope_err+','+intercept+'$\pm$'+intercept_err+','+K_cor_coe+','+K_pval+','+P_cor_coe+','+P_pval)
    
ax4.yaxis.set_ticks_position('both')
ax4.xaxis.set_ticks_position('both')
ax4.set_xlim(T04range)
ax4.set_ylim(w1brange)
ax4.errorbar(np.array(good_T04_array), good_w1_bpeakmags_excor_array-good_bb_bpeakmags_excor_array, xerr=0.1, yerr=np.sqrt(np.power(good_w1_bpeakmags_err_array,2.0)+np.power(good_bb_bpeakmags_err_array,2.0)),ls='')
ax4.scatter(np.array(good_T04_array), good_w1_bpeakmags_excor_array-good_bb_bpeakmags_excor_array)
#ax4.set_ylabel('w1-b Color')
#filter out NAN values
x1=good_T04_array
x1sig=np.array([0.1]*len(good_T04_array))
y1=good_w1_bpeakmags_excor_array-good_bb_bpeakmags_excor_array
y1sig=np.sqrt(np.power(good_w1_bpeakmags_err_array,2.0)+np.power(good_bb_bpeakmags_err_array,2.0))
x=[]
xsig=[]
y=[]
ysig=[]
for i in range(len(x1)):
	if str(x1[i])!='nan':
		if str(y1[i])!='nan':
			if str(x1sig[i])!='nan':
				if str(y1sig[i])!='nan':
					x.append(x1[i])
					y.append(y1[i])
					xsig.append(x1sig[i])
					ysig.append(y1sig[i])
corr_coeff, p_value_pearson=stats.pearsonr(x,y)
tau, p_value_kendall=stats.kendalltau(x,y)
lm= linmix.LinMix(x, y, xsig, ysig, K=2)
lm.run_mcmc(silent=True)
for i in range(0, len(lm.chain), 25):
	xs = np.arange(8,10)
	ys = lm.chain[i]['alpha'] + xs * lm.chain[i]['beta']
	ax4.plot(xs, ys, color='r', alpha=0.02)
intercept_array=[]
slope_array=[]
for i in range(0, len(lm.chain), 25):
	intercept=lm.chain[i]['alpha']
	slope=lm.chain[i]['beta']
	intercept_array.append(intercept)
	slope_array.append(slope)
mean_intercept=np.mean(intercept_array)
std_intercept=np.std(intercept_array)
mean_slope=np.mean(slope_array)
std_slope=np.std(slope_array)
x_mean=np.arange(8,10)
y_mean=mean_intercept+x_mean*mean_slope
ax4.plot(x_mean, y_mean, color='k')
#Record data
color='w1-b'
slope=str(np.round(mean_slope,3))
slope_err=str(np.round(std_slope,3))
intercept=str(np.round(mean_intercept,3))
intercept_err=str(np.round(std_intercept,3))
K_cor_coe=str(np.round(tau,3))
K_pval=str(np.round(p_value_kendall,3))
if float(K_pval)==0:
        K_pval=str(p_value_kendall)
P_cor_coe=str(np.round(corr_coeff,3))
P_pval=str(np.round(p_value_pearson,3))
if float(P_pval)==0:
	P_pval=str(p_value_pearson)
text_file.write('\n'+color+','+slope+'$\pm$'+slope_err+','+intercept+'$\pm$'+intercept_err+','+K_cor_coe+','+K_pval+','+P_cor_coe+','+P_pval)
    
ax6.yaxis.set_ticks_position('both')
ax6.xaxis.set_ticks_position('both')
ax6.set_ylim(m2brange)
ax6.set_xlim(T04range)
ax6.errorbar(np.array(good_T04_array), good_m2_bpeakmags_excor_array-good_bb_bpeakmags_excor_array, xerr=0.1, yerr=np.sqrt(np.power(good_m2_bpeakmags_err_array,2.0)+np.power(good_bb_bpeakmags_err_array,2.0)),ls='')
ax6.scatter(np.array(good_T04_array), good_m2_bpeakmags_excor_array-good_bb_bpeakmags_excor_array)
#ax6.set_ylabel('m2-b Color')
#filter out NAN values
x1=good_T04_array
x1sig=np.array([0.1]*len(good_T04_array))
y1=good_m2_bpeakmags_excor_array-good_bb_bpeakmags_excor_array
y1sig=np.sqrt(np.power(good_m2_bpeakmags_err_array,2.0)+np.power(good_bb_bpeakmags_err_array,2.0))
ax6.set_xlabel('log(OH) Metallicity T04')
x=[]
xsig=[]
y=[]
ysig=[]
for i in range(len(x1)):
	if str(x1[i])!='nan':
		if str(y1[i])!='nan':
			if str(x1sig[i])!='nan':
				if str(y1sig[i])!='nan':
					x.append(x1[i])
					y.append(y1[i])
					xsig.append(x1sig[i])
					ysig.append(y1sig[i])
corr_coeff, p_value_pearson=stats.pearsonr(x,y)
tau, p_value_kendall=stats.kendalltau(x,y)
lm= linmix.LinMix(x, y, xsig, ysig, K=2)
lm.run_mcmc(silent=True)
for i in range(0, len(lm.chain), 25):
	xs = np.arange(8,10)
	ys = lm.chain[i]['alpha'] + xs * lm.chain[i]['beta']
	ax6.plot(xs, ys, color='r', alpha=0.02)
intercept_array=[]
slope_array=[]
for i in range(0, len(lm.chain), 25):
	intercept=lm.chain[i]['alpha']
	slope=lm.chain[i]['beta']
	intercept_array.append(intercept)
	slope_array.append(slope)
mean_intercept=np.mean(intercept_array)
std_intercept=np.std(intercept_array)
mean_slope=np.mean(slope_array)
std_slope=np.std(slope_array)
x_mean=np.arange(8,10)
y_mean=mean_intercept+x_mean*mean_slope
ax6.plot(x_mean, y_mean, color='k')


#Record data
color='m2-b'
slope=str(np.round(mean_slope,3))
slope_err=str(np.round(std_slope,3))
intercept=str(np.round(mean_intercept,3))
intercept_err=str(np.round(std_intercept,3))
K_cor_coe=str(np.round(tau,3))
K_pval=str(np.round(p_value_kendall,3))
if float(K_pval)==0:
	K_pval=str(p_value_kendall)
P_cor_coe=str(np.round(corr_coeff,3))
P_pval=str(np.round(p_value_pearson,3))
if float(P_pval)==0:
	P_pval=str(p_value_pearson)
text_file.write('\n'+color+','+slope+'$\pm$'+slope_err+','+intercept+'$\pm$'+intercept_err+','+K_cor_coe+','+K_pval+','+P_cor_coe+','+P_pval)
    

#####################################

print('T04')
for q in range(len(good_T04_array)):
	if str(good_T04_array[q])!='nan': print(goodsnlist[q], good_D16_array[q], good_T04_array[q],good_uu_bpeakmags_excor_array[q]-good_bb_bpeakmags_excor_array[q], good_w1_bpeakmags_excor_array[q]-good_bb_bpeakmags_excor_array[q]) 



print('D16')
for q in range(len(good_D16_array)):
	if str(good_D16_array[q])!='nan': print(goodsnlist[q], good_D16_array[q], good_T04_array[q],good_uu_bpeakmags_excor_array[q]-good_bb_bpeakmags_excor_array[q], good_w1_bpeakmags_excor_array[q]-good_bb_bpeakmags_excor_array[q]) 



# f.suptitle('Color vs Mass and Metallicity T04',y=0.91)
plt.savefig('Color_vs_Mass_and_Metallicity_T04')
    




#plt.scatter(np.array(good_D16_array), good_w1_bpeakmags_excor_array-good_bb_bpeakmags_excor_array, s=area, alpha=0.5)
#plt.xlabel('Metallicity D16')
#plt.ylabel('uvw1 - b')
plt.show()
plt.savefig('w1-b_D16met.eps', format='eps')

