These codes are from Evan Batteas Feb 2023

used for making the plots -- specareaphotonstack contains a function for binning, a function for loading text from files, and a nested for loop which will plot the spectra of a file, the effective area of a filter, and the photon count through that filter, for each filter.
If you want to compare multiple spectra, just use multiple spectra in the spectrafiles and labellist.



photoncount contains the photoncount function which intakes a wavelength list, a flux list, and has options to intake a length of observation time, and a specific filter you want. It will output counts_array which contains 6 arrays which each contains the photon count (photons/s if no observation time is input *because the inputs are in ergs/s*), sumcounts which is the sum of all those arrays, amendedwavez which I mostly used for debugging purposes I think, and an index which can be used in counts_array to access the array of the filter you specified in the input.

photonstack contains 3 functions, 1 for binning data within a certain range, 1 for loading text from files, and 1 called spectophotstackplots which takes in a list of spectra files, a list of labels for them, and an optional list of observation times for those files. It should output a stacked plot containing however many plots as files you put in, and uses the figure sizing function I made to make it journal sized.

plotsizing contains the function setFigSizeInches which has parameters width (in # pts desired), fraction, which is optional and affects the width, and subplots which is also optional. It outputs the figure size needed for that # of pts in inches, which is what matPlotLib uses. You can also input 'journal' to width for it to automatically go to 242.26653 pts.


CountSpecArray does that but for 1 file, 1 filter. 

filtereffarea contains the function ltext which is used for loading text from files, and feffarea which creates a plot of the effective area of a filter given its wavelength list and area at each wavelength. 

normalizedata contains a ltext function and then normalizes data so the peak of the y axis is 1.
