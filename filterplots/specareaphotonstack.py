
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

spectralfiles = ['../spectra/SN2016ccj_peak_11fe_appended.dat', '../spectra/SN2011fe_peak_11fe_appended.dat', '../spectra/SN2022hrs_peak_11fe_appended.dat']
labellist = ['SN2016ccj','SN2011fe','SN2022hrs']

spectralfiles = ['../spectra/SN2011fe_peak_11fe_appended.dat', '../spectra/SN2016ccj_peak_11fe_appended.dat', '../spectra/SN2022hrs_peak_11fe_appended.dat']
labellist = ['SN2011fe','SN2016ccj','SN2022hrs']

spectralfiles = ['../spectra/SN2011fe_peak_11fe_appended.dat',  '../spectra/SN2022hrs_peak_11fe_appended.dat']
labellist = ['SN2011fe','SN2022hrs']


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
  binsize = (amendedwave[-1]-amendedwave[0])/641   # %%%%%%%%% generalize this for a different range  641 from 1600 to 8000 in bins of 10
  binmin = amendedwave[0]
  binmax = amendedwave[0] + binsize
  #need to have 641 wavelength indicies
  #special case if already binned and greater binsteps 
  #if((amendedwave[1] - amendedwave[0]) > binsize):
  if(len(amendedwave)<641):
    binwave = amendedwave
    binflux = amendedflux
    abmin = min(amendedwave)
    abmax = max(amendedwave)

    for x in range(641-len(binwave)):
      binwave.append(abmax+10)
      binflux.append(binflux[-1])
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
#  print(wavelist)
#  print(binwave)
#  wait = input("Press Enter to continue.")
  return binwave, binflux    

def loadtext(file,wavelist=[],fluxlist=[],range=8000):
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



############################

x = 0
#for sfile in spectralfiles:
# print(labellist[x])
y =0
for filter in filterfiles:
    print(filter)
    maxcounts1 = []
    maxcounts2 = []
    maxflux1 = []
    maxflux2 =[]
    # Creating empty lists to contain  data
    print("Inner Loop Begin")
    wavelength_spectra, wavelength_filters, flux_spectra, area_filters = [],[],[],[]
    wavelength_spectra, flux_spectra = loadtext(spectralfiles[0],wavelist=wavelength_spectra,fluxlist=flux_spectra)
 
    #print("Plot 1")
    fig, axes = plt.subplots(3,1,figsize=setFigSizeInches('journal', subplots=(3,1)))                                                
    axes[0].plot(wavelength_spectra, flux_spectra)                    # Limit is 6000     
    #print(wavelength_spectra==wavelength_filters)                               
    axes[0].set_xlabel('Wavelength [$\AA$ngstroms]')
    axes[0].set_ylabel('flux spectra [ergs/s/cm^2/$\AA$]')
    axes[0].set_xlim(1600,6000)
    axes[0].tick_params(direction="in")
    maxflux1 = np.nanmax(flux_spectra)
    plotlabelname=labellist[0]+'_'+filternames[y]+'.png'


    if len(spectralfiles) > 1:
       wavelength_spectra, wavelength_filters, flux_spectra, area_filters = [],[],[],[]
       wavelength_spectra, flux_spectra = loadtext(spectralfiles[1],wavelist=wavelength_spectra,fluxlist=flux_spectra)
       maxflux2 = np.nanmax(flux_spectra)
       scale = (maxflux1/maxflux2)
       flux_spectra=np.array(flux_spectra)
       scaled_flux = (flux_spectra*scale)
       axes[0].plot(wavelength_spectra, scaled_flux,linestyle='dashed')                   
       plotlabelname=labellist[0]+'_'+labellist[1]+'_'+filternames[y]+'.png'
     

       if len(spectralfiles) > 2:
          wavelength_spectra, wavelength_filters, flux_spectra, area_filters = [],[],[],[]
          wavelength_spectra, flux_spectra = loadtext(spectralfiles[2],wavelist=wavelength_spectra,fluxlist=flux_spectra)
          maxflux3 = np.nanmax(flux_spectra)
          scale = (maxflux1/maxflux3)
          scaled_flux = (np.array(flux_spectra)*scale)
          axes[0].plot(wavelength_spectra, scaled_flux,linestyle='dotted')                    # Limit is 6000  
          plotlabelname=labellist[0]+'_'+labellist[1]+'_'+labellist[2]+'_'+filternames[y]+'.png'
   


    wavelength_filters, area_filters = loadtext(filter,wavelist=wavelength_filters,fluxlist=area_filters)

    #print("Plot 2")
    #fig, axes = plt.subplots()                                                 
    axes[1].plot(wavelength_filters, area_filters, label=filternames[y])
    axes[1].set_xlabel('Wavelength [$\AA$ngstroms]')
    axes[1].set_ylabel('Effective Area [cm]')
    axes[1].set_xlim(1600,6000)
    axes[1].tick_params(direction="in")
    axes[1].legend()
    lowerbound = 1600
    upperbound = 8000 #5650


  # Third plot, plots wavelength_spectra against photon count
    specifiedfilter = filternames[y]
   
    print(specifiedfilter)
    wavelength_spectra, wavelength_filters, flux_spectra, area_filters = [],[],[],[]
    wavelength_spectra, flux_spectra = loadtext(spectralfiles[0],wavelist=wavelength_spectra,fluxlist=flux_spectra)
    wavelength_spectra, flux_spectra = rangebin(lowerbound,upperbound,wavelength_spectra, flux_spectra)
    counts_array, sumcounts, amend, countsindex = photoncount(wavelength_spectra,flux_spectra, specificfilter = specifiedfilter)
    maxcounts1 = np.nanmax(counts_array[countsindex])
    
    axes[2].plot(wavelength_spectra, counts_array[countsindex]/maxcounts1,label=labellist[0])
    axes[2].set_xlabel('Wavelength [$\AA$ngstroms]')
    axes[2].set_ylabel('Count Rate [photons/sec/$\AA$]')
    axes[2].set_xlim(1600,6000)
    axes[2].tick_params(direction="in")

#    for i in range(2):
#      print(i, x)
#      if(i==x):
#        continue
#      else:

    if len(spectralfiles) > 1:

       wavelength_spectra, wavelength_filters, flux_spectra, area_filters = [],[],[],[]
       wavelength_spectra, flux_spectra = loadtext(spectralfiles[1],wavelist=wavelength_spectra,fluxlist=flux_spectra)
       wavelength_spectra, flux_spectra = rangebin(lowerbound,upperbound,wavelength_spectra, flux_spectra)
       counts_array, sumcounts, amend, countsindex = photoncount(wavelength_spectra,flux_spectra, specificfilter = specifiedfilter)
       maxcounts2 = np.nanmax(counts_array[countsindex])
       axes[2].plot(wavelength_spectra, counts_array[countsindex]/maxcounts2,label=labellist[1], linestyle='dashed')

       if len(spectralfiles) > 2:
 
          wavelength_spectra, wavelength_filters, flux_spectra, area_filters = [],[],[],[]
          wavelength_spectra, flux_spectra = loadtext(spectralfiles[2],wavelist=wavelength_spectra,fluxlist=flux_spectra)
          wavelength_spectra, flux_spectra = rangebin(lowerbound,upperbound,wavelength_spectra, flux_spectra)
          counts_array, sumcounts, amend, countsindex = photoncount(wavelength_spectra,flux_spectra, specificfilter = specifiedfilter)
          maxcounts3 = np.nanmax(counts_array[countsindex])
          axes[2].plot(wavelength_spectra, counts_array[countsindex]/maxcounts3,label=labellist[2],linestyle='dotted')


    plt.legend()
    plt.savefig(plotlabelname, bbox_inches='tight', dpi=300)
    plt.show()
    y+=1
    x+=1
plt.close('fig')

   




