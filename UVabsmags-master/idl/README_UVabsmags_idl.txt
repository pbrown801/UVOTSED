README_UVabsmags_idl.txt -- this file -- a list of all the files in this directory 
Last updated Sep 20, 2018

SN2011fe_redbolmags161.sav  -- an array of UVOT magnitudes from the SN2011fe peak spectrum with different reddening laws and amounts of reddening applied

SN_abs_all6_forcebv.eps -- output of plots_snabs.pro showing the absolute magnitudes if you fix the distance and E(B-V) color excess by requiring a constant mag and color

bpeak_m2w1vw1v_spectra.eps
bpeak_w1vbv_spectra.eps

host.sav -- IDL save file containing all of the info from my spreadsheet as well as reading in the results of light curve fitting for the type Ia supernovae and any galaxy magnitudes.  This should be versioned when making a final analysis.

hstspectralist.txt -- a list of supernovae Ia with HST UV spectroscopy

m2w1_w1v_bpeak_nuvbr.eps

pjb_phot_array.pro -- reads in the SOUSA data files and lines up the magnitudes across filters within a given time bin for making colors
pjb_phot_array_B141.pro  -- the version of the above which I have been using for a while.  Given here for historical reproducability.

plots_snabs.pro --  reads in magnitudes, distances, and light curve shape parameters from host.sav in order to make some Phillips-style absolute magnitude versus deltaM15(B) plots

snabs6_loop.eps -- absolute magnitude versus deltaM15(B) plot with all 6 UVOT filters

snialist.txt -- list of all type Ia supernovae observed by Swift up to the start of the Swift-ASASSN volume-limited project

swiftspectralist.txt -- list of supernovae Ia with Swift UV grism spectroscopy (no cut on quality)

w1bv.pro
w1v_bluebpeak.eps
w1v_bpeak_colors_abs_color.eps
w1v_bpeak_colors_abs_dmb15.eps
w1v_bpeak_colors_abs_dmb15_nocolorbar.eps
