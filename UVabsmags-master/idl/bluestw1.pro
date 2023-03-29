pro bluestw1

fontsize=16

restore, '~/Desktop/Dropbox/SN/www/SwiftSN/host.sav'  

index11fe=where(host.snname_array eq 'SN2011fe')
print, host.dm15_array[index11fe,4] 

normal=where(host.sntype_array eq 'Ia' and host.dm15_array[*,4] gt 1.0 and host.dm15_array[*,4] lt 1.4 and (host.bpeakappmag_array[*,4]-host.bpeakappmag_array[*,5]) lt 0.3)

ia=where(host.sntype_array eq 'Ia')

w1vbluest_ia_array = make_array(6,n_elements(ia),value=!Values.F_NAN)
w1vbluest_ia_err_array = make_array(6,n_elements(ia),value=!Values.F_NAN)


milnenuvb=['SN2011ia','SN2011fe','SN2011by','SN2008hv','SN2008Q']
milnenuvr=['SN2015F','SN2013gy','SN2013gs','SN2013ex','SN2013cs','SN2012hr', 'SN2011im','SN2008ec','SN2007co','SN2007af','SN2005df','SN2005cf']
;;; adding in those too red
milnenuvr=['SN2015F','SN2013gy','SN2013gs','SN2013ex','SN2013cs','SN2012hr', 'SN2011im','SN2008ec','SN2007co','SN2007af','SN2005df','SN2005cf', 'SN2013eu','SN2010kg','SN2010gp','SN2010ev']

milnemuvb=['SN2012ht', 'SN2006ej','SN2006dm','SN2010gn']


;;;;;;;;;;;;;;;;;;;;;   w1-v


for n=0, n_elements(host.snname_array[ia]) -1 do begin

	SNname= host.snname_array[ia[n]]
	print, SNname
	if file_test('$SOUSA/data/'+SNname+'_uvotB15.1.dat') eq 1 then  begin

	pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat',   dt=dt

	fitsfilename='~/Desktop/Dropbox/SN/SOUSA/fitting/'+SNname+'_6fitsB15.sav'

	if file_test(fitsfilename) eq 1 then  begin
		restore, fitsfilename, verbose=0
		Bpeaktime=combbpeakmjd
		if SNname eq 'SN2005cf' or SNname eq 'SN2005df' or  SNname eq 'SN2011fe' then Bpeaktime=groundbpeakmjd
	endif else begin
		Bpeaktime=dt.bb[0,where(dt.bb[1,*] eq min(dt.bb[1,*],/nan))]
	endelse

	if SNname eq 'SN2005cf' or SNname eq 'SN2005df' or  SNname eq 'SN2011fe' then begin

		groundbbfile='$SOUSA/grounddata/'+SNname+'_bb_ground.dat'

		readcol,groundbbfile,bbmjd, bbmag, bbmagerr,/silent

		dt.mag_array[4,*]=interpol(bbmag,bbmjd,dt.time_array)
		dt.magerr_array[4,*]=interpol(bbmagerr,bbmjd,dt.time_array)

		groundvvfile='$SOUSA/grounddata/'+SNname+'_vv_ground.dat'

		readcol,groundvvfile,vvmjd, vvmag, vvmagerr,/silent

		dt.mag_array[5,*]=interpol(vvmag,vvmjd,dt.time_array)
		dt.magerr_array[5,*]=interpol(vvmagerr,vvmjd,dt.time_array)

	endif



	;if SNname eq 'SN2011B' then stop

	w1vbv=where( finite(dt.mag_array[2,*]) eq 1 and finite(dt.mag_array[4,*]) eq 1 and finite(dt.mag_array[5,*]) eq 1 and dt.time_array-Bpeaktime[0] lt 50.0,w1vbvcount)

	firstepoch = dt.time_array[0]
	
	if w1vbvcount gt 0  then w1vbluest=where( (dt.mag_array[2,w1vbv]-dt.mag_array[5,w1vbv]) eq min(dt.mag_array[2,w1vbv]-dt.mag_array[5,w1vbv]) ,w1vbluestcount)

	if w1vbvcount gt 0  then w1vbluest_ia_array[*,n]=dt.mag_array[*,w1vbv[w1vbluest]]
	if w1vbvcount gt 0  then w1vbluest_ia_err_array[*,n]=dt.magerr_array[*,w1vbv[w1vbluest]]

	endif

endfor

; from http://www.iluvatar.org/~dwijn/idlfigures
!p.font = 1
!p.thick = 2
!x.thick = 2
!y.thick = 2
!z.thick = 2
xsize = 8.8
wall = 0.04
margin=0.16
a = xsize/8.8 - (margin + wall)
b = a * 2d / (1 + sqrt(5))

ysize = (margin + b + wall )*8.8
ticklen = 0.01
xticklen = ticklen/b
yticklen = ticklen/a

nxticks=10

x1 = margin*8.8/xsize
x2 = x1 + a*8.8/xsize
xc = x2 + wall*8.8/xsize
y1 = margin*8.8/ysize
y2 = y1 + b*8.8/ysize
y3 = y2 + wall*8.8/ysize
y4 = y3 + b*8.8/ysize
yc = y4 + wall*8.8/ysize

xdata=[0,1,2,3]
ydata=[2,3,4,5]


plotfilename = 'bluest_ia_w1vbv.eps'

SET_PLOT, 'PS'

device, filename= plotfilename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

restore, 'SN2011fe_redbolmags161.sav'
; filters, ebv, epochs, laws
;  extinction laws model=0 rv=3.1, 1 1.7  , model 2 smc, model 3 CSLMC
; plot, epoch, feredmags[1,0,*,0]

febluest=where(min(feredmags[2,0,*,0]-feredmags[5,0,*,0]) eq feredmags[2,0,*,0]-feredmags[5,0,*,0],count)
febpeak=where(min(feredmags[4,0,*,0]) eq feredmags[4,0,*,0],count)

cgplot, charsize=1, feredmags[4,*,febluest,3]-feredmags[5,*,febluest,3], feredmags[2,*,febluest,3]-feredmags[5,*,febluest,3], xrange=[-0.2,0.2], yrange=[-0.2,2.5], xtitle='!A(b-v)!Nbluest       ', ytitle='(w1-v)!Bbluest', position=[x1,y1,x2,y2], linestyle=0, color=black
oplot, feredmags[4,*,febluest,2]-feredmags[5,*,febluest,2], feredmags[2,*,febluest,2]-feredmags[5,*,febluest,2], linestyle=1
oplot, feredmags[4,*,febluest,1]-feredmags[5,*,febluest,1], feredmags[2,*,febluest,1]-feredmags[5,*,febluest,1], linestyle=2
oplot, feredmags[4,*,febluest,0]-feredmags[5,*,febluest,0], feredmags[2,*,febluest,0]-feredmags[5,*,febluest,0], linestyle=3

;;;;;;;;;;;;;;;;;
oploterror, w1vbluest_ia_array[4,*]-w1vbluest_ia_array[5,*], w1vbluest_ia_array[2,*]-w1vbluest_ia_array[5,*], sqrt(w1vbluest_ia_err_array[4,*]^2.0+w1vbluest_ia_err_array[5,*]^2.0), sqrt(w1vbluest_ia_err_array[2,*]^2.0+w1vbluest_ia_err_array[5,*]^2.0), symsize=0.3, psym=15
 

device, /close
SET_PLOT, 'X'
$open bluest_ia_w1vbv.eps


plotfilename = 'bluest_ia_w1vbv_names.eps'

SET_PLOT, 'PS'

device, filename= plotfilename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

restore, 'SN2011fe_redbolmags161.sav'
; filters, ebv, epochs, laws
;  extinction laws model=0 rv=3.1, 1 1.7  , model 2 smc, model 3 CSLMC
; plot, epoch, feredmags[1,0,*,0]

febluest=where(min(feredmags[2,0,*,0]-feredmags[5,0,*,0]) eq feredmags[2,0,*,0]-feredmags[5,0,*,0],count)
febpeak=where(min(feredmags[4,0,*,0]) eq feredmags[4,0,*,0],count)

cgplot, charsize=1, feredmags[4,*,febluest,3]-feredmags[5,*,febluest,3], feredmags[2,*,febluest,3]-feredmags[5,*,febluest,3], xrange=[-0.3,0.1], yrange=[0.5,0.8], xtitle='!A(b-v)!Nbluest       ', ytitle='(w1-v)!Bbluest', position=[x1,y1,x2,y2], linestyle=0, color=black
oplot, feredmags[4,*,febluest,2]-feredmags[5,*,febluest,2], feredmags[2,*,febluest,2]-feredmags[5,*,febluest,2], linestyle=1
oplot, feredmags[4,*,febluest,1]-feredmags[5,*,febluest,1], feredmags[2,*,febluest,1]-feredmags[5,*,febluest,1], linestyle=2
oplot, feredmags[4,*,febluest,0]-feredmags[5,*,febluest,0], feredmags[2,*,febluest,0]-feredmags[5,*,febluest,0], linestyle=3

;;;;;;;;;;;;;;;;;
oploterror, w1vbluest_ia_array[4,*]-w1vbluest_ia_array[5,*], w1vbluest_ia_array[2,*]-w1vbluest_ia_array[5,*], sqrt(w1vbluest_ia_err_array[4,*]^2.0+w1vbluest_ia_err_array[5,*]^2.0), sqrt(w1vbluest_ia_err_array[2,*]^2.0+w1vbluest_ia_err_array[5,*]^2.0), symsize=0.3, psym=15

xyouts, w1vbluest_ia_array[4,*]-w1vbluest_ia_array[5,*] - 0.04, w1vbluest_ia_array[2,*]-w1vbluest_ia_array[5,*]+0.02, host.snname_array[ia], charsize=0.5

device, /close
SET_PLOT, 'X'
$open bluest_ia_w1vbv_names.eps



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




print, 'final stop'
stop

end
