pro peaklist
restore, 'host.sav'  
readcol, '../snnamesquality.txt', SNnamelist, quality, format='A,A'

peaksneia=SNnamelist[where(quality eq 'good' or quality eq 'peakdecline' or quality eq 'risepeak')]

cmagicsneia=SNnamelist[where(quality eq 'good' or quality eq 'peakdecline')]
cmagicsneiadmb15s=fltarr(n_elements(cmagicsneia))
cmagicsneiapeaktimes=fltarr(n_elements(peaksneia))


;peaksneia=host.SNname_array

peaktimes=fltarr(n_elements(peaksneia))
dmb15s=fltarr(n_elements(peaksneia))
m2_w1s=fltarr(n_elements(peaksneia))
w1_bs=fltarr(n_elements(peaksneia))
z=fltarr(n_elements(peaksneia))


for i=0,n_elements(peaksneia)-1 do begin
	SNindex=where(host.SNname_array eq peaksneia[i])
	peaktimes[i]=host.PEAKTIME_ARRAY[SNindex,4]
	dmb15s[i]=host.dm15_ARRAY[SNindex,4]
	m2_w1s[i]=host.BPEAKAPPMAG_ARRAY[SNindex,1]-host.BPEAKAPPMAG_ARRAY[SNindex,2]
	w1_bs[i] =host.BPEAKAPPMAG_ARRAY[SNindex,2]-host.BPEAKAPPMAG_ARRAY[SNindex,4]
;	print, peaksneia[i], peaktimes[i]

endfor

	slowest=sort(dmb15s)
	for i=0,n_elements(peaksneia)-1 do print, peaksneia[slowest[i]], dmb15s[slowest[i]],w1_bs[slowest[i]],m2_w1s[slowest[i]]
	print, "that was the slowest"

	print, " "

for i=0,n_elements(peaksneia)-1 do begin
	SNindex=where(host.SNname_array eq peaksneia[i])
	peaktimes[i]=host.PEAKTIME_ARRAY[SNindex,4]
	dmb15s[i]=host.dm15_ARRAY[SNindex,4]
	m2_w1s[i]=host.BPEAKAPPMAG_ARRAY[SNindex,1]-host.BPEAKAPPMAG_ARRAY[SNindex,2]
	w1_bs[i] =host.BPEAKAPPMAG_ARRAY[SNindex,2]-host.BPEAKAPPMAG_ARRAY[SNindex,4]
;	print, peaksneia[i], peaktimes[i]


	print, peaksneia[i], ' ',host.redshift_array[SNindex], ' ',host.PEAKTIME_ARRAY[SNindex,0], ' ',host.PEAKTIME_ARRAY[SNindex,1], ' ',host.PEAKTIME_ARRAY[SNindex,2], ' ',host.PEAKTIME_ARRAY[SNindex,3], ' ',host.PEAKTIME_ARRAY[SNindex,4], ' ',host.PEAKTIME_ARRAY[SNindex,5], ' ',host.APPMAG_ARRAY[SNindex,0], ' ',host.APPMAG_ARRAY[SNindex,1], ' ', host.APPMAG_ARRAY[SNindex,2], ' ',host.APPMAG_ARRAY[SNindex,3], ' ',host.APPMAG_ARRAY[SNindex,4], ' ', host.APPMAG_ARRAY[SNindex,5], format='(A,A,F5.3,A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,F8.2,A,F5.2,A,F5.2,A,F5.2,A,F5.2,A,F5.2,A,F5.2)'

endfor





for i=0,n_elements(cmagicsneia)-1 do begin
	SNindex=where(host.SNname_array eq cmagicsneia[i])
	cmagicsneiapeaktimes[i]=host.PEAKTIME_ARRAY[SNindex,4]
	cmagicsneiadmb15s[i]=host.dm15_ARRAY[SNindex,4]
endfor


for i=0,n_elements(peaksneia)-1 do begin
	SNindex=where(host.SNname_array eq peaksneia[i])
	peaktimes[i]=host.PEAKTIME_ARRAY[SNindex,4]
	dmb15s[i]=host.dm15_ARRAY[SNindex,4]
	m2_w1s[i]=host.BPEAKAPPMAG_ARRAY[SNindex,1]-host.BPEAKAPPMAG_ARRAY[SNindex,2]
	w1_bs[i] =host.BPEAKAPPMAG_ARRAY[SNindex,2]-host.BPEAKAPPMAG_ARRAY[SNindex,4]
;	print, peaksneia[i], peaktimes[i]


	print, peaksneia[i], ' ',host.APPMAG_ARRAY[SNindex,0], ' ',host.APPMAG_ARRAY[SNindex,1], ' ', host.APPMAG_ARRAY[SNindex,2], ' ',host.APPMAG_ARRAY[SNindex,3], ' ',host.APPMAG_ARRAY[SNindex,4], ' ', host.APPMAG_ARRAY[SNindex,5], format='(A,A,F5.2,A,F5.2,A,F5.2,A,F5.2,A,F5.2,A,F5.2)'

endfor



stop
end


