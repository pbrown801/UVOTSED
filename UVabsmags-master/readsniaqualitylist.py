import numpy as np
import os

f = open('snnamesquality.txt', 'r')

header1 = f.readline()
header2 = f.readline()


SNnames = []
coverages = []
comments = []
fittingstatus=[]
datafiles=[]

for line in f:
    line=line.strip()
    columns=line.split()
    number = len(columns)
    name = columns[0]
    coverage = columns[1]
    if number > 2:
        comment = columns[2]
        comments.append(comment)
    SNnames.append(name)
    coverages.append(coverage)

    if os.path.exists('/Users/pbrown/SN/github/SOUSA/data/'+name+'_uvotB15.1.dat'):
        datafiles.append('done')
    if not os.path.exists('/Users/pbrown/Desktop/SN/github/SOUSA/data/'+name+'_uvotB15.1.dat'):
        datafiles.append('___')
#        print(name+' needs SOUSA photometry')
#    if coverage == 'check':
#        os.system('open $SOUSA/lcplots/'+name+'_lightcurve.jpg')
    if os.path.exists('/Users/pbrown/Desktop/Dropbox/SN/SOUSA/fitting/'+name+'_6fits16.sav'):
        fittingstatus.append('done')
    else:
        fittingstatus.append(' ')
#
#    if not os.path.exists('/Users/pbrown/Desktop/Dropbox/SN/SOUSA/fitting/'+name+'_6fits16.sav') and os.path.exists('/Users/pbrown/Desktop/SN/github/SOUSA/data/'+name+'_uvotB15.1.dat':
#        print(name+' with coverage '+coverage+' needs to be fit')
    print(name, ' ', coverage, fittingstatus)

f.close()

#print(SNnames)
ngood=coverages.count('good')
print('total ',len(SNnames))
print('number of good ', coverages.count('good'))
print('number of peak ', coverages.count('peak'))
print('number of peakdecline ', coverages.count('peakdecline'))
print('number of decline ', coverages.count('decline'))
print('number to check ', coverages.count('check'))
print(fittingstatus)
