#make a separate program which takes as input a list of spectra and outputs a skinny, stacked plot with the photon count plot of each spectrum as a separate row
import numpy as np
import matplotlib.pyplot as plt
from plotsizing import *
from photoncount import photoncount


def rangebin(rangemin,rangemax,binsize, wavelist, fluxlist):
  global t
  print("Creating bins from ", rangemin, " to ", rangemax)
  
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
      #if(t==1):
        #print('in between')
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
  t = 0
  return binwave, binflux    

#function for file reading because np.loadtxt is being fussy
def ltext(filelist,wavelist,fluxlist,range=5950):
  x = 0
  for file in filelist:
   f = open('' +file,'r') #open file
   for line in f:
        line = line.rstrip()
        column = line.split() # split line into a list 
        wavelength = float(column[0]) #putting data into named variable for readability
        flux = float(column[1])
        #print('wavelist', x)
        if(float(wavelength) > range): #if the wavelength is out of range, it will ignore it, this is so the photoncount function can work
          continue
        wavelist[x].append(float(wavelength))
        fluxlist[x].append(float(flux))
   f.close()
   x+=1     
  return wavelist, fluxlist
#making function, takes a list of spectra files and optionally a list of observing times
def spectophotstackplots(spectralist,label,obstimelist=[]):
  print("stacking begin")
  #getting the number of plots to stack
  plotstackheight = len(spectralist)
  #loading in the wavelength and flux of each spectra
  print(plotstackheight, 'height')
  wavelist,fluxlist = [[] for x in spectralist],[[]for x in spectralist]
  print(wavelist, 'start')
  wavelist,fluxlist= ltext(spectralist,wavelist,fluxlist)
  #print(wavelist, '!!! wavelist')
  
  #creating the figure to plot on, using setFigSizeInches to make the sizing right for journals
  fig, axes = plt.subplots(plotstackheight,1,figsize=setFigSizeInches('journal', subplots=(plotstackheight,1)))  
  #doing the plotting
  for x in range(plotstackheight):
    #getting photon count using photoncount function
    print('x',x)
    #if(x == 0):
    
    #  print(fluxlist[x])
    wavelist[x], fluxlist[x] = rangebin(1600,5650, 10, wavelist[x], fluxlist[x])
    

    photons, sumphotons,amendwavelist,photonfilt  = photoncount(wavelist[x], fluxlist[x])
    
    #plotting the subplot
    if(len(amendwavelist)>2):
      #print(len(amendwavelist))
      print("amended?>>")
      axes[x].plot(amendwavelist,sumphotons)
    elif(len(sumphotons)>len(wavelist[x])):
      print("2nd")
      print('len: ', len(sumphotons), len(wavelist[x]) )
      print(wavelist[x])
      fixedphotons = [x in sumphotons for x in range(len(wavelist[x]))]
      axes[x].plot(wavelist[x],fixedphotons)
    else:
     print("3rd")
     '''
     print(wavelist[x])
     print(len((wavelist[x])), 'wave', len(sumphotons))
     lenfr = len((wavelist[x]))
     print(lenfr, 'fr')
     #try:
     f = open("photonizeddata.dat", "x")
      
     for i in range(0, lenfr):
       print(wavelist[x][i], sumphotons[i]) 
       f.write(str(wavelist[x][i]) + " " + str(float(sumphotons[i]))+ " ") #+ str(float(weightedunc[x])) + " " + str(int(weighteddq[x])))
       f.write("\n")  
     f.close()
     '''
     #except:
      #  pass  
     axes[x].plot(wavelist[x],sumphotons)
    axes[x].set_ylabel('photon count of ' +str(label[x][:9]))
    axes[x].set_yscale('log')
    axes[x].tick_params(direction="in")
  #making sure the xlabel is at the bottom
  axes[-1].set_xlabel('wavelength(Angstroms)')
  #making layout nice
  plt.tight_layout()
  #show
  plt.show()











t = 1
#then, to get it, just call function with list of files
#filelist = ['C:/Users/Evan/Desktop/Portfolio/Code/DrBProj/aggienova/spectra/SN2022hrs_muv_20220426.3_10.dat','C:/Users/Evan/Desktop/Portfolio/Code/DrBProj/aggienova/spectra/ptf11kly_20110907.obs.dat','C:/Users/Evan/Desktop/Portfolio/Code/DrBProj/aggienova/spectra/SN2016ccj_hst_20160514.dat']
#labellist= ['SN2016ccj_hst_20160514','ptf11kly_20110907','SN2022hrs_muv_20220426.3_10']
#spectophotstackplots(filelist, labellist)

file = ['C:\\Users\\Evan\\Desktop\\Portfolio\\Code\\DrBProj\\SN2011by_0509_full.dat']
label = ['attempt1']
spectophotstackplots(file,label)



#'C:/Users/Evan/Desktop/Portfolio/Code/DrBProj/aggienova/spectra/ptf11kly_20110907.obs.dat','C:/Users/Evan/Desktop/Portfolio/Code/DrBProj/aggienova/spectra/SN2016ccj_hst_20160514.dat','C:/Users/Evan/Desktop/Portfolio/Code/DrBProj/aggienova/spectra/SN2022hrs_muv_20220426.3_10.dat'
#Current problem, file is out of range of the filters, thus the # of wv datapoints in file far excceds the photon counts the photoncount funciton can create
#If WV greater than .... ignore or zero it