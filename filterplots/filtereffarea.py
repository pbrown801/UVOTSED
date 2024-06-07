import numpy as np
import matplotlib.pyplot as plt
from plotsizing import *
#program to plot normalized effective filter area
#filterfilelist = []
def ltext(filelist,wavelist,arealist):
  x = 0
  for file in filelist:
   f = open('' +file,'r')
   for line in f:
        line = line.rstrip()
        column = line.split()
        wavelength = column[0]
        flux = column[1]
        #print('wavelist', x)
        wavelist[x].append(float(wavelength))
        arealist[x].append(float(flux))
   f.close()
   x+=1     
  return wavelist, arealist

def feffarea(filterfilelist):
 
 wavelist, arealist = [[] for x in filterfilelist], [[] for x in filterfilelist]
 wavelist, arealist = ltext(filterfilelist, wavelist,arealist)
 for x in range(len(arealist)):
  for y in range(len(arealist[x])):
   arealist[x][y] = arealist[x][y]/max(arealist[x])
   
 fig, axes = plt.subplots(len(filterfilelist),1,figsize=setFigSizeInches('journal', subplots=(len(filterfilelist),1)))  
 for x in range(len(wavelist)):
    axes[x].plot(wavelist[x],arealist[x])
    #axes[x].set_ylabel('Scaled Effective Area(%)')
    axes[x].set_xlabel(flist[x][-13:-4])
    axes[x].tick_params(direction="in")
    
 axes[0].set_title('Scaled Effective Area(%)')
 plt.tight_layout()
 plt.show()   

flist = ['C:/Users/Evan/Desktop/Portfolio/Code/DrBProj/aggienova/filters/UVW1_2010.txt','C:/Users/Evan/Desktop/Portfolio/Code/DrBProj/aggienova/filters/UVM2_2010.txt']
feffarea(flist)