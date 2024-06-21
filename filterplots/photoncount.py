
#modifided this code to instead of finding the photon count per filter, find the photon per wv by multiplying by each element area of the filter per wavelength
def photoncount(wavez, fluxz,obstime=1,specificfilter='U_UVOT'):
    ''' 
     takes input of a wavelength and flux file, as well as observation time (default to 1 so its actaully photons/s) and a specific filter if wanted (default of U_UVOT)
    '''

    import numpy as np
    import matplotlib.pyplot as plt

    #plank constant and speed of light
    h = 6.6260755e-27
    c = 2.99792458e18
    hc = h*c


    
    #filter files
    files = ['../filters/UVW2_2010.txt','../filters/UVM2_2010.txt','../filters/UVW1_2010.txt','../filters/U_UVOT.txt','../filters/B_UVOT.txt', '../filters/V_UVOT.txt']
    
    x = 0
    filter_lambda = []
    filter_area = []
    #getting data from files
    for item in files:
        #Necessary to have "../" when running in /python/ directory
        f = open("../" + item,'r')

#	print(item)
        #seperating each filter's values so i can later apply them seperately
        filter_lambda.append([])
        filter_area.append([])
        for line in f:
         line = line.rstrip()
         column = line.split()
#		print(column)
         wavelen = column[0]
         area = column[1]
         filter_lambda[x].append(float(wavelen))
         filter_area[x].append(float(area))
        x+=1

        f.close()



    ##########################################


    filtercurves = ['UVW2_2010','UVM2_2010','UVW1_2010','U_UVOT','B_UVOT','V_UVOT'] ### STRING LIST
    filternums   = [0,1,2,3,4,5]

    filterdict = {'UVW2_2010':0,'UVM2_2010':1,'UVW1_2010':2, 'U_UVOT':3,'B_UVOT':4,'V_UVOT':5 } #dictory so you know which index (of counts_array) to access for which filter you want to see photon count through
    #leftovers
    zeropoints = [17.38, 16.85, 17.44, 18.34, 19.11, 17.89] ### PHOTOMETRIC ZEROPOINTS BASED ON VEGA

    #leftovers
    filtereffwavelength=[2030,2231,2634,3501,4329,5402] ### EFFECTIVE VEGA WAVELENGTH FOR EACH FILTER (IN SAME ORDER)
    
    counts_array = []
    #for loop in which first, a new list is created inside the counts array list, and the photon counts for each wavelength value of one filter is added, then a new filter with a new list does the same thing, for all filters
    #after for loop terminantes, add all of the individual lists together to get a total photon count over all filters
    
    #instead of for x in range of filter lambda, can do x in range len(wavez), that ought to solve the overproduction problem
    for x in range(len(filter_lambda)):
     counts_array.append([])
     counts_array[x] = np.zeros(len(filter_lambda[x]))
     for y in range(len(filter_lambda[x])):
      #print("Before", type(fluxz[y]))    
      #print(y)
      if(y>(len(wavez)-1)):    
        break
      else:          
       counts_array[x][y] = (filter_area[x][y] * fluxz[y] * obstime) / (((h*c)/wavez[y]*(10**-10))*10**7) #np.trapz( , wavez) #area of the filter at this wv, times flux at this wv, times wv at this wv / plank * c
    counts_array = np.array(counts_array)
    #sumcounts = np.array(sumcounts)
    sumcounts = np.zeros(len(counts_array[0]))
    for x in range(len(counts_array)):
       sumcounts = sumcounts + counts_array[x]
    admendedwavez = []
    if(len(sumcounts)<len(wavez)):
      for i in range(len(sumcounts)):
        admendedwavez.append(wavez[i])   
    return counts_array, sumcounts,admendedwavez, filterdict[specificfilter]   #amendedwavez shouldnt need to be used, filterdict takes input of filter name to return an index for counts_array
