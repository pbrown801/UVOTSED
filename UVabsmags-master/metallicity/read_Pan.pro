pro read_Pan


readcol, 'Pan_2019_tab1.tex', SNNames, z, dm15, dm15err, Phase, Hostnames,    Ebv_host, Morph, logM, logMperr, logMmerr,   logOH, logOHperr, logOHmerr,   AGNflag, $
format='A, F, F, F, F, A, F, A, F, F, F, F, F, F, A'

order=sort(logOH)


for i=0, n_elements(logM)-1 do print, SNNames[order[i]], z[order[i]], dm15[order[i]], dm15err[order[i]], Phase[order[i]], ' ', Hostnames[order[i]],    Ebv_host[order[i]], ' ', Morph[order[i]], logM[order[i]], logMperr[order[i]], logMmerr[order[i]],   logOH[order[i]], logOHperr[order[i]], logOHmerr[order[i]],   ' ', AGNflag[order[i]]


readcol, 'Pan_2019_mass_f2535.dat', mass, f2535
readcol, 'Pan_2019_mass_f2535_err.dat', mass, f2535_err

for i=0, n_elements(mass)-1 do print, mass[i], f2535[i], f2535_err[i]-f2535[i]

readcol, 'Pan_2019_mass_f3025.dat', mass, f3025
readcol, 'Pan_2019_mass_f3025_err.dat', mass, f3025_err

for i=0, n_elements(mass)-1 do print, mass[i], f3025[i], f3025_err[i]-f3025[i]

readcol, 'Pan_2019_tab1_fadded.dat', SNNames, z, dm15, dm15err, Phase, Hostnames,    Ebv_host, Morph, logM, logMperr, logMmerr,   logOH, logOHperr, logOHmerr,   AGNflag, f2535, f2535_err, f3025, f3025_err, $
format='A, F, F, F, F, A, F, A, F, F, F, F, F, F, A, F, F, F, F'



print, 'final stop'
stop
end