
import numpy as np
import matplotlib.pyplot as plt
from speccounts import *
from plotsizing import setFigSizeInches
#from photoncountperwv import *
#from photonstack import spectophotstackplots
from photoncount import photoncount
#from a list of 3 spectra 2022hrs, 2016ccj, ptf11kly
# for each filter in 
filterfiles = ['../filters/UVW2_2010.txt','../filters/UVM2_2010.txt','../filters/UVW1_2010.txt','../filters/U_UVOT.txt','../filters/B_UVOT.txt', '../filters/V_UVOT.txt']
filternames = ['UVW2_2010','UVM2_2010','UVW1_2010','U_UVOT','B_UVOT','V_UVOT']
spectralfiles = ['../spectra/SN2022hrs_muv_20220426.3_10.dat','../spectra/ptf11kly_20110907.obs.dat','../spectra/SN2016ccj_hst_20160514.dat']
labellist = ['SN2016ccj_hst_20160514','ptf11kly_20110907','SN2022hrs_muv_20220426.3_10']
# make a stacked plot 3 plots high
# top plot: spectra and flux
# middle plot: effective area of the filter
# bottom plot: photon count of all 3 files through that filter

# how to execute:
# 1. for loop (for each file):
#  2. for loop(for each filter)
#   3. plot spectra, put in top area of figure
#   4. plot effective area curve of filter put in middle area of figure
#   5. run spectophotstackplots on the filelist, and plot all the photon counts on bottom area of figure with different colors
#   6. profit


def rangebin(rangemin,rangemax,wavelist, fluxlist):
  #print("Creating bins from ", rangemin, " to ", rangemax)
  
  amendedwave = []
  amendedflux = []
  binwave = []
  binflux = []
  sum = [0,0]
  numtimes = 0
  #range selection
  
  for x in range(len(wavelist)):
    if((wavelist[x] >= rangemin) and wavelist[x] <= rangemax  ): 
      amendedwave.append(wavelist[x])
      amendedflux.append(fluxlist[x])
  #print(len(amendedflux))
  #binning
  binsize = (amendedwave[-1]-amendedwave[0])/641
  binmin = amendedwave[0]
  binmax = amendedwave[0] + binsize
  #need to have 641 wavelength indicies
  #special case if already binned and greater binsteps 
  #if((amendedwave[1] - amendedwave[0]) > binsize):
  if(len(amendedwave)<641):
    binwave = amendedwave
    binflux = amendedflux
    abmin = min(amendedwave)
    for x in range(641-len(binwave)):
      binwave.append(abmin-1)
      binflux.append(0)
  else:    
   for x in range(len(amendedwave)):
    if(binmin<amendedwave[x]<binmax):
      sum[0]+=amendedflux[x]
      sum[1]+=1 
    elif(amendedwave[x]>binmax):
      binwave.append(binmin)
      alf = 0
      if(sum[1] < 1):
        sum[1] = 1
      if(sum[0] == 0):
        fzero = x
        try:
           sum[0] = prevsum
        except:
          pass   
      numtimes +=1  
      binflux.append((sum[0]/sum[1]))
      prevsum = sum[0]/sum[1]
      sum[0] = 0
      sum[1] = 0
      binmin = binmax
      binmax = binmin + binsize
  if(len(binwave)<641):
    abmin = min(binwave)
    for x in range(641-len(binwave)):
      binwave.insert(0,abmin-1)
      binflux.insert(0,0)  
  if(len(binwave)>641):
    while(len(binwave)>641):
        del(binwave[640])
        del(binflux[640])

  return binwave, binflux    

def loadtext(file,wavelist=[],fluxlist=[],range=5950):
  #x = 0
  #for file in filelist:
  f = open('' +file,'r') #open file
  for line in f:
        line = line.rstrip()
        column = line.split() # split line into a list 
        wavelength = float(column[0]) #putting data into named variable for readability
        flux = float(column[1])
        #print('wavelist', x)
        if(float(wavelength) > range): #if the wavelength is out of range, it will ignore it, this is so the photoncount function can work
          continue
        wavelist.append(float(wavelength))
        fluxlist.append(float(flux))
  f.close()
   #x+=1     
  return wavelist, fluxlist

x = 0
for sfile in spectralfiles:
 y =0
 for filter in filterfiles:
    # Creating empty lists to contain  data
    print("Inner Loop Begin")
    wavelength_spectra, wavelength_filters, flux_spectra, area_filters = [],[],[],[]
    wavelength_spectra, flux_spectra = loadtext(sfile,wavelist=wavelength_spectra,fluxlist=flux_spectra)
    try:
     print(wavelength_spectra==wavelength_filters)
    except:
        pass 
    
     

    #print("Plot 1")
    fig, axes = plt.subplots(3,1,figsize=setFigSizeInches('journal', subplots=(3,1)))                                                
    axes[0].plot(wavelength_spectra, flux_spectra)                    # Limit is 6000     
    #print(wavelength_spectra==wavelength_filters)                               
    axes[0].set_xlabel('wavelength spectra')
    axes[0].set_ylabel('flux spectra (ergs/s/cm^2)')
    axes[0].set_xlim(1600,6000)
    axes[0].tick_params(direction="in")
   
    wavelength_filters, area_filters = loadtext(filter,wavelist=wavelength_filters,fluxlist=area_filters)
    try:
     print(wavelength_spectra==wavelength_filters)
    except:
        pass 
    #print("Plot 2")
    #fig, axes = plt.subplots()                                                 
    axes[1].plot(wavelength_filters, area_filters)
    axes[1].set_xlabel('wavelength_filters')
    #axes[1].set_xlabel('wavelength_filters')
    axes[1].set_ylabel('area (cm)')
    axes[1].set_xlim(1600,6000)
    axes[1].tick_params(direction="in")
    try:
     print(wavelength_spectra==wavelength_filters)
    except:
        pass 
   

  # Third plot, plots wavelength_spectra against photon count
    specifiedfilter = filternames[y]
    #print(specifiedfilter)
    wavelength_spectra, flux_spectra = rangebin(1600,5650,wavelength_spectra, flux_spectra)
    counts_array, sumcounts, amend, countsindex = photoncount(wavelength_spectra,flux_spectra, specificfilter = specifiedfilter)
    axes[2].plot(wavelength_spectra, counts_array[countsindex])
    axes[2].set_xlabel('wavelength')
    axes[2].set_ylabel('photon count (photons)')
    axes[2].set_xlim(1600,6000)
    axes[2].tick_params(direction="in")
    try:
     print(wavelength_spectra==wavelength_filters)
    except:
        pass 

    for i in range(2):
      #print("Third For")
      if(i==x):
        continue
      else:
        wavelength_spectra, flux_spectra = loadtext(spectralfiles[x])
        wavelength_spectra, flux_spectra = rangebin(1600,5650,wavelength_spectra, flux_spectra)
        counts_array, sumcounts, amend, countsindex = photoncount(wavelength_spectra,flux_spectra, specificfilter = specifiedfilter)
        axes[2].plot(wavelength_spectra, counts_array[countsindex])
    plt.tight_layout()
    plt.show()
    plt.close('all')
    try:
     print(wavelength_spectra==wavelength_filters)
    except:
        pass 
    wavelength_spectra, wavelength_filters, flux_spectra, area_filters = [],[],[],[]
    y+=1
x+=1




