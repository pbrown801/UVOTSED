;;;; array part needs other columns added
;;; THIS IS THE HISTORICAL ONE AND SHOULD NOT BE MODIFIED!!

PRO pjb_phot_read_B141, infile, x

;;;; this function reads in the data from my SN_uvotB14.1.dat files
;;;; with data separated by filter and listed with columns 
;;;; # Filter MJD[days] Mag MagErr 3SigMagLim 0.98SatLim[mag] Rate[c/s] RateErr[c/s] Ap[arcsec] Frametime[s] Exp[s] Telapse[s]
;;;;                                     
;;;;
;;;; The output structure dt has four different formats for the output
;;;; dt.mag_array has the magnitudes from all 6 filters grouped by mjd epoch dt.time_array
;;;; dt.magerr_array, dt.counts_array, dt.countserr_array, dt.maglimits_array are similar
;;;; each filter has an array dt.w2 with mjd, mag, magerr, maglimits, counts, countserr
;;;; as well as dt.w2time, dt.w2mags, dt.w2mags_err, dt.w2maglimits,  dt.w2counts, dt.w2counts_err
;;;; the epoch array can also be output as a fits file and maybe eventually as a tex table

GET_LUN, u
OPENR, u, infile

printfirst=0

b=''
i=0
x = DBLARR(10000,12)
x = make_array(10000,12,/DOUBLE, value=!values.f_nan)
WHILE NOT EOF(u) DO BEGIN
READF, u, b
a1 = ''
a2 = ''
a3 = ''
a4 = ''
a5 = ''
a6 = ''
a7 = ''
a8 = ''
a9 = ''
a10 = ''
a11 = ''
a12 = ''
a13 = ''

; # Filter MJD[days] Mag MagErr 3SigMagLim 0.98SatLim[mag] Rate[c/s] RateErr[c/s] Ap[arcsec] Frametime[s] Exp[s] Telapse[s]
; #                                                                                           
; UVW2     53525.0428  17.790   0.114  20.313  11.085   0.685   0.072  3.0  0.0110     283.23     287.76


;;; skip lines beginning with #
IF STRPOS(b, '#') lt 0 THEN BEGIN
READS, b, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, format='(a9,a10,a8, a8, a8, a8, a8, a8, a5, a11, a11, a11)'
if (printfirst eq 0) then print, b
printfirst=1

;;; replace NULLs with NaNs

   if strpos(a3, 'NULL') ge 0 then a3=!VALUES.F_NAN
   if strpos(a4, 'NULL') ge 0 then a4=!VALUES.F_NAN
   if strpos(a7, 'NULL') ge 0 then a7=!VALUES.F_NAN
   if strpos(a8, 'NULL') ge 0 then a8=!VALUES.F_NAN


if a1 eq 'UVW2     ' then    x[i,0] = 0
if a1 eq 'UVM2     ' then    x[i,0] = 1
if a1 eq 'UVW1     ' then    x[i,0] = 2
if a1 eq 'U        ' then    x[i,0] = 3
if a1 eq 'B        ' then    x[i,0] = 4
if a1 eq 'V        ' then    x[i,0] = 5
;print, a1, x[i,0]
;; if you want to convert to Julian Date from MJD
;   x[i,1] = a2 + 2400000.5
   x[i,1] = a2 
   x[i,2] = a3
   x[i,3] = a4
   x[i,4] = a5
   x[i,5] = a6
   x[i,6] = a7
   x[i,7] = a8
   x[i,8] = a9
   x[i,9] = a10
   x[i,10] = a11
   x[i,11] = a12
;   x[i,12] = a13
;   x[i,13] = a7
;   x[i,7] = f

if x[i,2] eq 0.0 then    x[i,2] = !VALUES.F_NAN
if x[i,3] eq 0.0 then    x[i,3] = !VALUES.F_NAN
if x[i,6] eq  0.0 then    x[i,6] = !VALUES.F_NAN
if x[i,7] eq  0.0 then    x[i,7] = !VALUES.F_NAN


endif
   i = i+1
ENDWHILE

CLOSE, u
FREE_LUN, u
w = WHERE(x[*,1] GT 0.0)
x = x[w,*]

;RETURN, x



END

;;;;;;;;;;;;;;;;;;

PRO pjb_phot_array_B141, infile, dt=dt, fitsfile=fitsfile, delta_t=delta_t, tex6file=tex6file

;This program turns a list of photometry times, magnitudes and errors
;into a binary fits table where photometric points from multiple
;filters are grouped together by time.  The size of these time bins
;can be adjusted with DELTA_T.  This program calls SNPHOT_READ

if keyword_set(delta_t) ne 1 then delta_t = 0.15                ;Half width of time bin in days

pjb_phot_read_B141, infile,x

n = N_ELEMENTS(x[*,0])

;;; sort on mjd and then select which rows correspond to which filter
s = SORT(x[*,1])
x = x[s,*]

w2spots=where(x[*,0] eq 0.0)
m2spots=where(x[*,0] eq 1.0)
w1spots=where(x[*,0] eq 2.0)
uuspots=where(x[*,0] eq 3.0)
bbspots=where(x[*,0] eq 4.0)
vvspots=where(x[*,0] eq 5.0)
whspots=where(x[*,0] eq 6.0)


; # Filter MJD[days] Mag MagErr 3SigMagLim 0.98SatLim[mag] Rate[c/s] RateErr[c/s] Ap[arcsec] Frametime[s] Exp[s] Telapse[s]
; #                                                                                           
; UVW2     53525.0428  17.790   0.114  20.313  11.085   0.685   0.072  3.0  0.0110     283.23     287.76

ncolumns=11 ; not counting filter

dtarray = CREATE_STRUCT('TIME_ARRAY', make_array(n,/DOUBLE, value=!values.f_nan),  $
'MAG_ARRAY',make_array(6,n,/FLOAT, value=!values.f_nan), 'MAGERR_ARRAY', make_array(6,n,/FLOAT, value=!values.f_nan), $
 'COUNTS_ARRAY', make_array(6,n,/FLOAT, value=!values.f_nan),  'COUNTSERR_ARRAY', make_array(6,n,/FLOAT, value=!values.f_nan), $
 'MAGUPLIMIT_ARRAY',  make_array(6,n,/FLOAT, value=!values.f_nan), $
 'MAGSATLIMIT_ARRAY',  make_array(6,n,/FLOAT, value=!values.f_nan) , $
 'AP_ARRAY',  make_array(6,n,/FLOAT, value=!values.f_nan) , $
 'Frametime_ARRAY',  make_array(6,n,/FLOAT, value=!values.f_nan) , $
 'exptime_ARRAY',  make_array(6,n,/FLOAT, value=!values.f_nan) , $
 'telapse_ARRAY',  make_array(6,n,/FLOAT, value=!values.f_nan)  )

i = 0
j = 0

WHILE i LT n DO BEGIN


   t = x[i,1]
   w = WHERE(ABS(x[*,1] - t) LT delta_t, nw)

   dtarray.time_array[j] = MEAN(x[w,1])

   FOR k=0, nw-1 DO BEGIN

	for f=0,5 do begin

      IF x[w[k],0] EQ f THEN BEGIN
	 dtarray.mag_array[f,j] = float(x[w[k],2])
         dtarray.magerr_array[f,j] = x[w[k],3]
         dtarray.maguplimit_array[f,j] = x[w[k],4]
         dtarray.magsatlimit_array[f,j] = x[w[k],5]
         dtarray.counts_array[f,j] = x[w[k],6]
         dtarray.countserr_array[f,j] = x[w[k],7]
          dtarray.ap_array[f,j] = x[w[k],8]
         dtarray.frametime_array[f,j] = x[w[k],9]
         dtarray.exptime_array[f,j] = x[w[k],10]
         dtarray.telapse_array[f,j] = x[w[k],11]
      ENDIF

	ENDFOR

   ENDFOR

i=w[nw-1]+1
   j = j+1
ENDWHILE


dt = CREATE_STRUCT('TIME', make_array(n,/DOUBLE, value=!values.f_nan), $
'W2', make_array(ncolumns,n_elements(w2spots),/DOUBLE, value=!values.f_nan), $
'M2', make_array(ncolumns,n_elements(m2spots),/DOUBLE, value=!values.f_nan),  $
'W1', make_array(ncolumns,n_elements(w1spots),/DOUBLE, value=!values.f_nan),  $
'UU', make_array(ncolumns,n_elements(uuspots),/DOUBLE, value=!values.f_nan), $
'BB', make_array(ncolumns,n_elements(bbspots),/DOUBLE, value=!values.f_nan),  $
'VV', make_array(ncolumns,n_elements(vvspots),/DOUBLE, value=!values.f_nan),  $
'WH', make_array(ncolumns,n_elements(whspots),/DOUBLE, value=!values.f_nan),  $ 
'W2TIME', make_array(n_elements(w2spots),/DOUBLE, value=!values.f_nan), $
'M2TIME', make_array(n_elements(m2spots),/DOUBLE, value=!values.f_nan),  $
'W1TIME', make_array(n_elements(w1spots),/DOUBLE, value=!values.f_nan),  $
'UUTIME', make_array(n_elements(uuspots),/DOUBLE, value=!values.f_nan), $
'BBTIME', make_array(n_elements(bbspots),/DOUBLE, value=!values.f_nan),  $
'VVTIME', make_array(n_elements(vvspots),/DOUBLE, value=!values.f_nan),  $
'WHTIME', make_array(n_elements(whspots),/DOUBLE, value=!values.f_nan),  $ 
'W2COUNTS', make_array(n_elements(w2spots),/FLOAT, value=!values.f_nan), 'W2COUNTS_ERR', make_array(n_elements(w2spots),/FLOAT, value=!values.f_nan), $
'M2COUNTS', make_array(n_elements(m2spots),/FLOAT, value=!values.f_nan), 'M2COUNTS_ERR', make_array(n_elements(m2spots),/FLOAT, value=!values.f_nan),  $
'W1COUNTS', make_array(n_elements(w1spots),/FLOAT, value=!values.f_nan),  'W1COUNTS_ERR', make_array(n_elements(w1spots),/FLOAT, value=!values.f_nan), $
'UUCOUNTS', make_array(n_elements(uuspots),/FLOAT, value=!values.f_nan), 'UUCOUNTS_ERR', make_array(n_elements(uuspots),/FLOAT, value=!values.f_nan), $
'BBCOUNTS', make_array(n_elements(bbspots),/FLOAT, value=!values.f_nan), 'BBCOUNTS_ERR', make_array(n_elements(bbspots),/FLOAT, value=!values.f_nan),  $
'VVCOUNTS', make_array(n_elements(vvspots),/FLOAT, value=!values.f_nan), 'VVCOUNTS_ERR', make_array(n_elements(vvspots),/FLOAT, value=!values.f_nan),  $
'WHCOUNTS', make_array(n_elements(whspots),/FLOAT, value=!values.f_nan), 'WHCOUNTS_ERR', make_array(n_elements(whspots),/FLOAT, value=!values.f_nan),  $
'W2MAGS', make_array(n_elements(w2spots),/FLOAT, value=!values.f_nan), 'W2MAGS_ERR', make_array(n_elements(w2spots),/FLOAT, value=!values.f_nan), $
'M2MAGS', make_array(n_elements(m2spots),/FLOAT, value=!values.f_nan), 'M2MAGS_ERR', make_array(n_elements(m2spots),/FLOAT, value=!values.f_nan),  $
'W1MAGS', make_array(n_elements(w1spots),/FLOAT, value=!values.f_nan),  'W1MAGS_ERR', make_array(n_elements(w1spots),/FLOAT, value=!values.f_nan), $
'UUMAGS', make_array(n_elements(uuspots),/FLOAT, value=!values.f_nan), 'UUMAGS_ERR', make_array(n_elements(uuspots),/FLOAT, value=!values.f_nan), $
'BBMAGS', make_array(n_elements(bbspots),/FLOAT, value=!values.f_nan), 'BBMAGS_ERR', make_array(n_elements(bbspots),/FLOAT, value=!values.f_nan),  $
'VVMAGS', make_array(n_elements(vvspots),/FLOAT, value=!values.f_nan), 'VVMAGS_ERR', make_array(n_elements(vvspots),/FLOAT, value=!values.f_nan),  $
'WHMAGS', make_array(n_elements(whspots),/FLOAT, value=!values.f_nan), 'WHMAGS_ERR', make_array(n_elements(whspots),/FLOAT, value=!values.f_nan),  $
'W2UPLIMIT', make_array(n_elements(w2spots),/FLOAT, value=!values.f_nan),  'M2UPLIMIT', make_array(n_elements(m2spots),/FLOAT, value=!values.f_nan),  $
'W1UPLIMIT', make_array(n_elements(w1spots),/FLOAT, value=!values.f_nan), 'UUUPLIMIT', make_array(n_elements(uuspots),/FLOAT, value=!values.f_nan),  $
'BBUPLIMIT', make_array(n_elements(bbspots),/FLOAT, value=!values.f_nan), 'VVUPLIMIT',  make_array(n_elements(vvspots),/FLOAT, value=!values.f_nan),  $
'WHUPLIMIT', make_array(n_elements(whspots),/FLOAT, value=!values.f_nan), $
 'W2SATLIMIT', make_array(n_elements(w2spots),/FLOAT, value=!values.f_nan),  'M2SATLIMIT', make_array(n_elements(m2spots),/FLOAT, value=!values.f_nan),  $
'W1SATLIMIT', make_array(n_elements(w1spots),/FLOAT, value=!values.f_nan), 'UUSATLIMIT', make_array(n_elements(uuspots),/FLOAT, value=!values.f_nan),  $
'BBSATLIMIT', make_array(n_elements(bbspots),/FLOAT, value=!values.f_nan), 'VVSATLIMIT',  make_array(n_elements(vvspots),/FLOAT, value=!values.f_nan),  $
'WHSATLIMIT', make_array(n_elements(whspots),/FLOAT, value=!values.f_nan), $ 
'TIME_ARRAY', make_array(j,/DOUBLE, value=!values.f_nan),  $
'MAG_ARRAY',make_array(6,j,/FLOAT, value=!values.f_nan), 'MAGERR_ARRAY', make_array(6,j,/FLOAT, value=!values.f_nan), $
 'COUNTS_ARRAY', make_array(6,j,/FLOAT, value=!values.f_nan),  'COUNTSERR_ARRAY', make_array(6,j,/FLOAT, value=!values.f_nan), $
 'MAGUPLIMIT_ARRAY',  make_array(6,j,/FLOAT, value=!values.f_nan), 'MAGSATLIMIT_ARRAY',  make_array(6,j,/FLOAT, value=!values.f_nan)  )



if w2spots[0] ge 0 then for q=0,n_elements(w2spots)-1 do dt.W2[*,q]=x[w2spots[q],1:11]
if m2spots[0] ge 0 then for q=0,n_elements(m2spots)-1 do dt.M2[*,q]=x[m2spots[q],1:11]
if w1spots[0] ge 0 then for q=0,n_elements(w1spots)-1 do dt.W1[*,q]=x[w1spots[q],1:11]
if uuspots[0] ge 0 then for q=0,n_elements(uuspots)-1 do dt.UU[*,q]=x[uuspots[q],1:11]
if bbspots[0] ge 0 then for q=0,n_elements(bbspots)-1 do dt.BB[*,q]=x[bbspots[q],1:11]
if vvspots[0] ge 0 then for q=0,n_elements(vvspots)-1 do dt.VV[*,q]=x[vvspots[q],1:11]
if whspots[0] ne -1 then for q=0,n_elements(whspots)-1 do dt.WH[*,q]=x[whspots[q],1:11]

dt.w2time=dt.w2[0,*]
dt.w2mags=dt.w2[1,*]
dt.w2mags_err=dt.w2[2,*]
dt.w2uplimit=dt.w2[3,*]
dt.w2satlimit=dt.w2[4,*]
dt.w2counts=dt.w2[5,*]
dt.w2counts_err=dt.w2[6,*]

dt.m2time=dt.m2[0,*]
dt.m2mags=dt.m2[1,*]
dt.m2mags_err=dt.m2[2,*]
dt.m2uplimit=dt.m2[3,*]
dt.m2satlimit=dt.m2[4,*]
dt.m2counts=dt.m2[5,*]
dt.m2counts_err=dt.m2[6,*]

dt.w1time=dt.w1[0,*]
dt.w1mags=dt.w1[1,*]
dt.w1mags_err=dt.w1[2,*]
dt.w1uplimit=dt.w1[3,*]
dt.w1satlimit=dt.w1[4,*]
dt.w1counts=dt.w1[5,*]
dt.w1counts_err=dt.w1[6,*]

dt.uutime=dt.uu[0,*]
dt.uumags=dt.uu[1,*]
dt.uumags_err=dt.uu[2,*]
dt.uuuplimit=dt.uu[3,*]
dt.uusatlimit=dt.uu[4,*]
dt.uucounts=dt.uu[5,*]
dt.uucounts_err=dt.uu[6,*]

dt.bbtime=dt.bb[0,*]
dt.bbmags=dt.bb[1,*]
dt.bbmags_err=dt.bb[2,*]
dt.bbuplimit=dt.bb[3,*]
dt.bbsatlimit=dt.bb[4,*]
dt.bbcounts=dt.bb[5,*]
dt.bbcounts_err=dt.bb[6,*]

dt.vvtime=dt.vv[0,*]
dt.vvmags=dt.vv[1,*]
dt.vvmags_err=dt.vv[2,*]
dt.vvuplimit=dt.vv[3,*]
dt.vvsatlimit=dt.vv[4,*]
dt.vvcounts=dt.vv[5,*]
dt.vvcounts_err=dt.vv[6,*]

dt.whtime=dt.wh[0,*]
dt.whmags=dt.wh[1,*]
dt.whmags_err=dt.wh[2,*]
dt.whuplimit=dt.wh[3,*]
dt.whsatlimit=dt.wh[4,*]
dt.whcounts=dt.wh[5,*]
dt.whcounts_err=dt.wH[6,*]

dt.time_array=dtarray.time_array[0:j-1]
for f=0,5 do dt.mag_array[f,*]=dtarray.mag_array[f,0:j-1]
for f=0,5 do dt.magerr_array[f,*]=dtarray.magerr_array[f,0:j-1]
for f=0,5 do dt.counts_array[f,*]=dtarray.counts_array[f,0:j-1]
for f=0,5 do dt.countserr_array[f,*]=dtarray.countserr_array[f,0:j-1]
for f=0,5 do dt.maguplimit_array[f,*]=dtarray.maguplimit_array[f,0:j-1]
for f=0,5 do dt.magsatlimit_array[f,*]=dtarray.magsatlimit_array[f,0:j-1]

;;; create fits table of epoch matched photometry
if n_elements(fitsfile) ne 0 then begin
	cmd = 'rm ' + fitsfile
	if (file_test(fitsfile) eq 1) then SPAWN, cmd
	
	MWRFITS, dt, fitsfile
endif


if keyword_set(tex6file) then begin

	openw, lun, tex6file, /get_lun, width=600

	print, 'printing ', tex6file
	;for k=0,j-1 do print, FORMAT = '((F8.2), 12( " & ", F5.2, :), (A) )', dt.time_array[k], dt.mag_array[0,k], dt.magerr_array[0,k], dt.mag_array[1,k],  dt.magerr_array[1,k], dt.mag_array[2,k],  dt.magerr_array[2,k],  dt.mag_array[3,k],  dt.magerr_array[3,k],  dt.mag_array[4,k], dt.magerr_array[4,k],  dt.mag_array[5,k],  dt.magerr_array[5,k], ' \\'

	for k=0,j-1 do printf, lun, FORMAT = '((F8.2), 12( " & ", F5.2, :), (A) )', dt.time_array[k], dt.mag_array[0,k], dt.magerr_array[0,k], dt.mag_array[1,k],  dt.magerr_array[1,k], dt.mag_array[2,k],  dt.magerr_array[2,k],  dt.mag_array[3,k],  dt.magerr_array[3,k],  dt.mag_array[4,k], dt.magerr_array[4,k],  dt.mag_array[5,k],  dt.magerr_array[5,k], ' \\'

	close, lun
	free_lun, lun

endif


END
