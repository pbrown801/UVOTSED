import sed_checker
import pandas as pd




targetlist  = pd.read_csv('SNPY_Sample_Decline_color2.csv')

for i in range(len(targetlist)):
	sed_checker.sed_checker(targetlist.sname[i])

# sed_checker.sed_checker('SN2011by')