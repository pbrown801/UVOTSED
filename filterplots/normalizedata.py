def ltext(file,wavelist,arealist):
  f = open(file,'r')
  for line in f:
        line = line.rstrip()
        column = line.split()
        wavelength = column[0]
        flux = column[1]
        #print('wavelist', x)
        wavelist.append(float(wavelength))
        arealist.append(float(flux))
  f.close()   
  return wavelist, arealist

combfiles = ['C:\\Users\\Evan\\Desktop\\Portfolio\\Code\\DrBProj\\SN_2011by_2011-05-10_05-48-28_Lick-3m_KAST_UCB-SNDB_0.flm','C:\\Users\\Evan\\Desktop\\Portfolio\\Code\\DrBProj\\aggienova\\pythoncode\\photonizeddata.dat']    #use g430 and g230, just a little overlap

file1 = combfiles[0]
print(file1, '!')
print('????????')
file2 = combfiles[1]
wavelist, arealist = [],[]
wavelist, arealist = ltext(file2, wavelist,arealist)
supermax = max(arealist)
for x in range(len(arealist)):
   arealist[x] = arealist[x]/supermax

try:
     f = open("photonnormaldata2.dat", "x")
     for i in range(0, len(wavelist)):
       #print(wavelist[x][i], sumphotons[i]) 
       f.write(str(wavelist[i]) + " " + str(float(arealist[i]))+ " ") #+ str(float(weightedunc[x])) + " " + str(int(weighteddq[x])))
       f.write("\n")  
     f.close()
    
except:
      pass 