pro stritblue




;restore, '~/Desktop/Dropbox/SN/www/SwiftSN/host.sav'  

;  this is the version from Brown et al. 2017 paper with Nancy Landez
restore, 'nancyhost.sav'  


index11fe=where(host.snname_array eq 'SN2011fe')
print, host.dm15_array[index11fe,4] 

normal=where(host.sntype_array eq 'Ia' and host.dm15_array[*,4] gt 1.0 and host.dm15_array[*,4] lt 1.4 and (host.bpeakappmag_array[*,4]-host.bpeakappmag_array[*,5]) lt 0.3)
 
bluest_array = make_array(6,n_elements(normal),value=!Values.F_NAN)
bluest_err_array = make_array(6,n_elements(normal),value=!Values.F_NAN)

w1vbluest_array = make_array(6,n_elements(normal),value=!Values.F_NAN)
w1vbluest_err_array = make_array(6,n_elements(normal),value=!Values.F_NAN)

bpeak_array=make_array(6,n_elements(normal),value=!Values.F_NAN)
bpeak_err_array=make_array(6,n_elements(normal),value=!Values.F_NAN)
firstepoch_array=make_array(n_elements(normal),value=!Values.F_NAN)
hstspectra_array=make_array(n_elements(normal),value=' ')
swiftspectra_array=make_array(n_elements(normal),value=' ')
nuvb_array=make_array(n_elements(normal),value=' ')
nuvr_array=make_array(n_elements(normal),value=' ')
milne_array=make_array(n_elements(normal),value=' ')

readcol,"swiftspectralist.txt",swiftspectra,format='A',/silent
readcol,"hstspectralist.txt",hstspectra,format='A',/silent


milnenuvb=['SN2011ia','SN2011fe','SN2011by','SN2008hv','SN2008Q']
milnenuvr=['SN2015F','SN2013gy','SN2013gs','SN2013ex','SN2013cs','SN2012hr', 'SN2011im','SN2008ec','SN2007co','SN2007af','SN2005df','SN2005cf']
;;; adding in those too red
milnenuvr=['SN2015F','SN2013gy','SN2013gs','SN2013ex','SN2013cs','SN2012hr', 'SN2011im','SN2008ec','SN2007co','SN2007af','SN2005df','SN2005cf', 'SN2013eu','SN2010kg','SN2010gp','SN2010ev']

milnemuvb=['SN2012ht', 'SN2006ej','SN2006dm','SN2010gn']


;;;;;;;;;;;;
symbol=[15,17,14,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6]



;;;;;;;;;;;;



;;;;;;;;;;;;;;;;;;;;;   w1-v


for n=0, n_elements(host.snname_array[normal]) -1 do begin

	SNname= host.snname_array[normal[n]]
	print, SNname

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

	if firstepoch - Bpeaktime lt -5.0 and w1vbv[0] ne -1 and w1vbluestcount gt 0 then w1vbluest_array[*,n]=dt.mag_array[*,w1vbv[w1vbluest]]
	if firstepoch - Bpeaktime lt -5.0 and w1vbv[0] ne -1 and w1vbluestcount gt 0 then w1vbluest_err_array[*,n]=dt.magerr_array[*,w1vbv[w1vbluest]]
	if firstepoch - Bpeaktime lt -5.0 and w1vbv[0] ne -1 and w1vbluestcount gt 0 then w1vbluest_err_array[*,n]=dt.magerr_array[*,w1vbv[w1vbluest]]

	if firstepoch - Bpeaktime gt -5.0 and w1vbvcount gt 0 then print, 'not early enough data for ', SNname

	firstepoch_array[n]=firstepoch-Bpeaktime[0]

	bpeak_array[*,n]=host.bpeakappmag_array[normal[n],*]
	bpeak_err_array[*,n]=host.bpeakappmagerr_array[normal[n],*]
	;for n=0,n_elements(normal)-1 do print, host.snname_array[normal[n]], w1vbluest_array[*,n], host.dm15_array[normal[n],4], firstepoch_array[n]


	swiftindex=where(swiftspectra eq SNname)
	hstindex  =where(hstspectra eq SNname)
	if swiftindex[0] ne -1 then swiftspectra_array[n]='swift'
	if hstindex[0] ne -1 then hstspectra_array[n]='hst'




	muvbindex  =where(milnemuvb eq SNname)
	nuvbindex  =where(milnenuvb eq SNname)
	nuvrindex  =where(milnenuvr eq SNname)
	if nuvbindex[0] ne -1 then milne_array[n]='nuvb'
	if muvbindex[0] ne -1 then milne_array[n]='muvb'
	if nuvrindex[0] ne -1 then milne_array[n]='nuvr'


endfor




restore, 'host.sav'  


readcol, 'Stritzinger_2018_redblue.txt', strit_SNname, strit_Host, strit_Redshift, strit_EBV_MW, strit_EBV_host, pm3, strit_ebvhosterr, strit_t_first, strit_t_rise, strit_DeltamB15, pm, dmerr, strit_M_B, pm2, strit_mberr, strit_Spectralcode, strit_SpectralType, strit_Color, References, format='(A, A, A, A, A,A,A,A,A,A, A, A, A, A, A, A, A, A, A)', comment='#'

stritbpeak_array=make_array(6,n_elements(strit_SNname),value=!Values.F_NAN)
stritbpeak_err_array=make_array(6,n_elements(strit_SNname),value=!Values.F_NAN)
stritcolor_array=make_array(n_elements(strit_SNname),value=' ')


for n=0, n_elements(strit_SNname) -1 do begin

	SNname= strit_SNname[n]
	print, SNname
	stritindex=where(host.snname_array eq SNname)	

	Bpeaktime=host.bpeakmjd_array[stritindex]

	stritbpeak_array[*,n]=host.bpeakappmag_array[stritindex,*]
	stritbpeak_err_array[*,n]=host.bpeakappmagerr_array[stritindex,*]	
	if stritindex[0] ne -1 then stritcolor_array[n]=strit_color[n]
endfor


;;;;


muvb=where(milne_array eq 'muvb')
nuvb=where(milne_array eq 'nuvb')
nuvr=where(milne_array eq 'nuvr')

stritred=where(stritcolor_array eq 'red')
stritblue=where(stritcolor_array eq 'blue')


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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


plotfilename = 'bpeak_w1vbv_strit.eps'

SET_PLOT, 'PS'

device, filename= plotfilename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

restore, 'SN2011fe_redbolmags161.sav'
; filters, ebv, epochs, laws
;  extinction laws model=0 rv=3.1, 1 1.7  , model 2 smc, model 3 CSLMC
; plot, epoch, feredmags[1,0,*,0]
febpeak=where(min(feredmags[4,0,*,0]) eq feredmags[4,0,*,0],count)


cgplot, charsize=1, feredmags[4,*,febpeak,3]-feredmags[5,*,febpeak,3], feredmags[2,*,febpeak,3]-feredmags[5,*,febpeak,3], xrange=[-0.15,0.4], yrange=[0.8,2.2], ystyle=1, xstyle=1, ytitle='(w1-v)!BB!Lpeak', $ 
xtitle=' !A (b-v)!NBpeak   ', $
; double subscripts falling off page
; xtitle='!S!U (b-v) !N B !R!I peak', $
position=[x1,y1,x2,y2], linestyle=0, color=black
;  not getting this to work
;xyouts, 0.02, 0.5, '(b-v) !R!I B peak'

oplot, feredmags[4,*,febpeak,2]-feredmags[5,*,febpeak,2], feredmags[2,*,febpeak,2]-feredmags[5,*,febpeak,2], linestyle=1
oplot, feredmags[4,*,febpeak,1]-feredmags[5,*,febpeak,1], feredmags[2,*,febpeak,1]-feredmags[5,*,febpeak,1], linestyle=2
oplot, feredmags[4,*,febpeak,0]-feredmags[5,*,febpeak,0], feredmags[2,*,febpeak,0]-feredmags[5,*,febpeak,0], linestyle=3

;;;;;;;;;;;;;;;;;
oploterror, bpeak_array[4,*]-bpeak_array[5,*], bpeak_array[2,*]-bpeak_array[5,*], sqrt(bpeak_err_array[4,*]^2.0+bpeak_err_array[5,*]^2.0), sqrt(bpeak_err_array[2,*]^2.0+bpeak_err_array[5,*]^2.0), symsize=0.3, psym=15



cgoplot, stritbpeak_array[4,stritred]-stritbpeak_array[5,stritred], stritbpeak_array[2,stritred]-stritbpeak_array[5,stritred], psym=16, symsize=1, color='red'

cgoplot, stritbpeak_array[4,stritblue]-stritbpeak_array[5,stritblue], stritbpeak_array[2,stritblue]-stritbpeak_array[5,stritblue], psym=46, symsize=1.2, color='blue'

;xyouts, stritbpeak_array[4,*]-bpeak_array[5,*] - 0.06, bpeak_array[2,*]-bpeak_array[5,*]+0.02, strit_SNname, charsize=0.5

al_legend, ['Stritzinger+2018 red','Stritzinger+2018 blue'], psym=[16,46], color=['red','blue'], symsize=[1,1.2], $
pos=[0.15,1.2], charsize=0.8, box=1

device, /close
SET_PLOT, 'X'
spawn, 'open bpeak_w1vbv_strit.eps'

;;;;;;;;;;;;;;;;;;;;;;;;
xrange=[-5,45]


nxticks=10


xdata=[2,4,6]
ydata=[6,8,9]


fontsize=12

SET_PLOT, 'PS'

device, filename='earlysne_strit.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize, bits_per_pixel=8, /color

x=0
plot, xdata,ydata, /nodata, /noerase, position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle=' Days from Maximum Light',   ytitle='uvw1 mag - uvw1 max', charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xrange=[-15,15],yrange=[6,-1], ystyle=1, xstyle=1

color=strit_Color

print, n_elements(n_elements(strit_SNname))

for n=0, n_elements(strit_SNname) -1 do begin

	SNname= strit_SNname[n]
	print, SNname
	stritindex=where(host.snname_array eq SNname)	

	if stritindex gt -1 then begin

	Bpeaktime=host.bpeakmjd_array[stritindex]

pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat', dt=dt

good=where(finite(dt.w1[1,*]) eq 1)
maxindex=where(dt.w1[1,good] eq min(dt.w1[1,good]))
w1maxtime=dt.w1[0,good[maxindex]]
w1maxmag=dt.w1[1,good[maxindex]]


oploterror, dt.w1[0,where(finite(dt.w1[1,*]) eq 1)]-w1maxtime[0],dt.w1[1,where(finite(dt.w1[1,*]) eq 1)]-w1maxmag[0],dt.w1[2,where(finite(dt.w1[1,*]) eq 1)], psym=-symbol[n] , color=color[n], errcolor=color[n], linestyle=0, symsize=0.5

	endif
endfor

stop

al_legend, strit_SNname, color=color[0:n_elements(strit_SNname)-1], psym=symbol[0:n_elements(strit_SNname)-1], position=[10,-1], box=1, charsize=0.6, linestyle=1, background='white', symsize=0.5



device, /close
SET_PLOT, 'X'
$open earlysne_strit.eps





print, 'final stop'
stop

end
