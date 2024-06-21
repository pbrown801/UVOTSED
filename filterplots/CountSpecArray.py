
import numpy as np
import matplotlib.pyplot as plt
from speccounts import *
from plotsizing import setFigSizeInches
from photoncountperwv import *

 #option for inputing the file you want through input, rathering than hardcoding
#spectra_file_input = input("What spectra file do you want to use?:  (use full path) ")       # User inputs spectra file name including extension
#spectra_file_name = ('%s' % (spectra_file_input))                # Path is relative, might need to be changed for different computers 

spectra_file_name = 'C:/Users/Evan/Desktop/Portfolio/Code/DrBProj/aggienova/spectra/SN2011fe_uv.dat'
spectra_file_load = np.loadtxt(spectra_file_name)

 #option for inputing the file you want through input, rathering than hardcoding
#filter_file_input = input("What filter file do you want to use?: (use name in ../filters/) ")         # User inputs filter file name including extension
#filter_file_name = ('../filters/%s' % (filter_file_input))                  # Path is relative again, change as needed 

filter_file_name = ('C:/Users/Evan/Desktop/Portfolio/Code/DrBProj/aggienova/filters/UVW1_2010.txt')                  # Path is relative again, change as needed 
filter_file_load = np.loadtxt(filter_file_name)

 # Loading text in from the file, first 2 columns

wave,flux = np.loadtxt(spectra_file_name,dtype=float,usecols=(0,1),unpack=True)

 #function from speccounts to get photon counts
  #specin_countsout(wave,flux)

 # Print debugging
#print('mag array')
#print(mag_array)

 # Creating empty lists to contain  data
 
wavelength_spectra = []
flux_spectra = []
wavelength_filters = []
area_filters = []    

 # Reads in the data from the files, appends to the empty lists above

for i in range(0, len(spectra_file_load)):                                  
    wavelength_spectra.append(spectra_file_load[i][0])
for i in range(0,len(spectra_file_load)):
    flux_spectra.append(spectra_file_load[i][1])
for i in range(0, len(filter_file_load)):
    wavelength_filters.append(filter_file_load[i][0])
for i in range(0, len(filter_file_load)):
    area_filters.append(filter_file_load[i][1])

  # First plot, plots wavelength of the spectra against the flux of the spectra

fig, axes = plt.subplots(3,1,figsize=setFigSizeInches('journal', subplots=(3,1)))                                                
#plt.figure(figsize=setFigSizeInches('journal')) #use a function to set the size of the figure
axes[0].plot(wavelength_spectra, flux_spectra)
plt.xlim(1600,6000)                    # Limit is 6000                                   
plt.xlabel('wavelength spectra')
plt.ylabel('flux spectra')
#axes[0].set_xlabel('wavelength spectra')
axes[0].set_ylabel('flux spectra (ergs/s/cm^2)')
axes[0].set_xlim(1600,6000)
axes[0].tick_params(direction="in")
plt.tick_params(direction="in")
plt.tight_layout()


 # Second plot, plots wavelength filters against the area filters 

#fig, axes = plt.subplots()                                                 
plt.figure(figsize=setFigSizeInches('journal'))
axes[1].plot(wavelength_filters, area_filters)
plt.xlim(1600,6000)
plt.xlabel('wavelength_filters')
plt.ylabel('area_filters')
#axes[1].set_xlabel('wavelength_filters')
axes[1].set_ylabel('area (cm)')
axes[1].set_xlim(1600,6000)
axes[1].tick_params(direction="in")
plt.tick_params(direction="in")
plt.tight_layout()

 # Third plot, plots wavelength_spectra against photon count
counts_array, sumcounts, amend, countsindex = photoncount(wavelength_spectra,flux_spectra)
plt.figure(figsize=setFigSizeInches('journal'))
axes[2].plot(wavelength_filters, counts_array[countsindex])
plt.xlim(1600,6000)
axes[2].set_xlabel('wavelength')
axes[2].set_ylabel('photon count (photons)')
axes[2].set_xlim(1600,6000)
axes[2].tick_params(direction="in")

plt.tight_layout()



 # Interpolates wavelength and flux using the wavelngth filters, then calculates ergs using that data * the area filters

interpolated_data = np.interp(wavelength_filters, wavelength_spectra, flux_spectra)
ergs = np.array(interpolated_data) * np.array(area_filters)

 # Integrates ergs(wavelength_filters) using trapezodial rule

integral = np.trapz(ergs, wavelength_filters)

 # Sets values for percents of the area, e.g ten_percent is 10 percent of the total area, etc.

ten_percent = (integral) * 0.1
fifty_percent = (integral) * 0.5
ninety_percent = (integral) * 0.9

 #The following is a series of for loops in which values (of wvfilters and ergs) for 10%, 50%, and 90% of the area are put into lists so they can later be compared
ten_percent_y_vals = []
ten_percent_x_vals = []
integral_test = 0

 # Calculates 10% of the total area 

for i in range(0, len(wavelength_spectra)):                                
    if integral_test < ten_percent:
        ten_percent_y_vals.append(ergs[i])
        ten_percent_x_vals.append(wavelength_filters[i])
        integral_test = np.trapz(ten_percent_y_vals, ten_percent_x_vals)
    else:
        break

fifty_percent_y_vals = []
fifty_percent_x_vals = []
integral_test = 0

  # Calculates 50% of the total area 

for i in range(0, len(wavelength_spectra)):                               
    if integral_test < fifty_percent:
        fifty_percent_y_vals.append(ergs[i])
        fifty_percent_x_vals.append(wavelength_filters[i])
        integral_test = np.trapz(fifty_percent_y_vals, fifty_percent_x_vals)
    else:
        break
        
ninety_percent_y_vals = []
ninety_percent_x_vals = []
integral_test = 0

   # Calculates 90% of the total area 

for i in range(0, len(wavelength_spectra)):                              
    if integral_test < ninety_percent:
        ninety_percent_y_vals.append(ergs[i])
        ninety_percent_x_vals.append(wavelength_filters[i])
        integral_test = np.trapz(ninety_percent_y_vals, ninety_percent_x_vals)
    else:
        break
'''

  # Plots wavelength filters vs the ergs

#fig, axes = plt.subplots()                                                
plt.figure(figsize=setFigSizeInches('journal'))
axes[2].plot(wavelength_filters, ergs)
plt.xlabel('wavelength_filters')
plt.ylabel('ergs/s/angstrom')
plt.tick_params(direction="in")
plt.tight_layout()
'''
# Start pulling data from a second spectrum, not required for the file to work 
 #Again, an option to input a file rather than hardcode
#spectra_file_input_2 = input("What is the second spectra file that you want to use? ")
spectra_file_input_2 = '../spectra/vega.dat'
spectra_file_name_2 = ('../spectra/%s' % (spectra_file_input_2))
spectra_file_load_2 = np.loadtxt(spectra_file_name_2)

#Repeat everything previously done with the new data
 #empty lists to store the data in
wavelength_spectra_2 = []
flux_spectra_2 = []

for i in range(0, len(spectra_file_load_2)):
    wavelength_spectra_2.append(spectra_file_load_2[i][0])
for i in range(0,len(spectra_file_load_2)):
    flux_spectra_2.append(spectra_file_load_2[i][1])

 #interpolate the new wv and flux with the filters 
interpolated_data_2 = np.interp(wavelength_filters, wavelength_spectra_2, flux_spectra_2)
ergs_2 = np.array(interpolated_data_2) * np.array(area_filters)

fig2 = plt.figure(figsize=setFigSizeInches('journal'))
#Stacking the first and second ergs on top of each other
ax=fig2.add_subplot(111, label=1)
ax2=fig2.add_subplot(111, label=2, frame_on=False)

#1st line on plot
ax.plot(wavelength_filters, ergs, color = 'C0', label=spectra_file_name)
#setting labels for 1st line on plot
ax.set_xlabel('wavelength_filter', color = 'C0')
ax.set_ylabel('ergs_1', color = 'C0')
#setting tick paramaters for 1st line on plot
ax.tick_params(axis='x', colors="C0", direction="in")
ax.tick_params(axis='y', colors="C0", direction="in")

#2nd line on plot
ax2.plot(wavelength_filters, ergs_2, color = 'C3', label=spectra_file_name_2)
#setting labels for 2nd line on plot
ax2.set_xlabel('wavelength_filter', color = 'C3')
ax2.set_ylabel('ergs_2', color = 'C3')
ax2.xaxis.set_label_position('top') #change these label positions or sides, and also size of tick text, they are squished
ax2.yaxis.set_label_position('right')
#setting tick params for 2nd line on plot
ax2.xaxis.tick_top()
ax2.yaxis.tick_right()
ax2.tick_params(axis='x', colors="C3", direction="in", pad = -2)
ax2.tick_params(axis='y', colors="C3", direction="in", pad= -10)

# mpatches not working  plot1_legend = mpatches.Patch(color='C0', label=spectra_file_input)
#plot2_legend = mpatches.Patch(color='C3', label=spectra_file_input_2)
# plt.legend(handles=[plot1_legend,plot2_legend])

# Showing the plot
plt.tight_layout()
plt.show()