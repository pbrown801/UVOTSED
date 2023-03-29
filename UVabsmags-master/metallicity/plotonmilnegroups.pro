pro plotonmilnegroups, SNnames, plotfilename=plotfilename, absplotfilename=absplotfilename, peakdates=peakdates

;plotonmilnegroups, ['SN2005cf','SN2005df', 'SN2008Q','SN2009Y', 'SN2009an', 'SN2009ig', 'SN2011ao','SN2011by', 'SN2011fe',  'SN2012cg', 'SN2012fr', 'SN2012ht','SN2013aa', 'SN2013dy', 'iPTF14bdn', 'ASASSN-14lp','SN2015F', 'SN2016ccz','SN2016ekg','SN2017cbv'], plotfilename='PanSample_colors.eps', absplotfilename='PanSample_abs.eps' 

;plotonmilnegroups, ['SN2013aa', 'SN2017cbv', 'SN2011by', 'SN2011fe'], plotfilename='IaUV_twins_colors.eps', absplotfilename='IaUV_twins_abs.eps' 
;$open IaUV_twins_colors.eps

if keyword_set(plotfilename) eq 0 then filename=SNnames[0]+'_milne_colors.eps'
if keyword_set(absplotfilename) eq 0 then filename=SNnames[0]+'_milne_abs.eps'

restore,  filename='milnecolors.sav'
;restore, '~/Desktop/Dropbox/SN/www/SwiftSN/host.sav'  
restore, 'host.sav'  

redSNlist=['SN2005cf', 'SN2007af', 'SN2005df', 'SN2007co', 'SN2007cv', 'SN2007sr', 'SN2008ec', 'SN2009an']
blueSNlist=['SN2011fe','SN2011by', 'SN2008Q', 'SN2008hv', 'SNF20080514002','SN2006ej']

shortnormallist=['SN2005cf','SN2007af', 'SN2011by','SN2011fe']
shortcolors=['red','pink','powder blue','royal blue']

xrange=[-20.0,30.0]
yrangew1v=[0.5,3.5]
yrangem2w1=[0,4.5]
yrangem2vv=[2.5,6.5]

w1absrange=[-12.0,-21.0]
m2absrange=[-11.0,-20.0]
vvabsrange=[-14.0,-23.0]




colors=['blue','red','dark green','blue','red','dark green','blue','red','dark green', $
'blue','red','dark green','blue','red','dark green','blue','red','dark green','blue','red','dark green','blue','red','dark green','blue','red','dark green','blue','red','dark green','blue','red','dark green' ] 
colors=['black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black','black' ] 
symnum=[14,15,16,4,5,6,9,7,11,23,25,27,29,35,37,41,14,15,16,4,5,6,9,7,11,23,25,27,29,35,37,41,14,15,16,4,5,6,9,7,11,23,25,27,29,35,37,41]
symnum=[14,15,16,14,15,16,14,15,16,14,15,16,14,15,16,14,15,16,14,15,16,14,15,16,14,15,16,14,15,16,14,15,16,14,15,16,14,15,16,14,15,16,14,15,16]
sizes=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
sizes=[0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5]

lines=[0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3]

nSNe=n_elements(SNnames)
symbols=symnum[0:nSNe-1]
sizes=sizes[0:nSNe-1]
colors=colors[0:nSNe-1]
lines=lines[0:nSNe-1]

loadct, 11

colortable=intarr(nSNe)


; from http://www.iluvatar.org/~dwijn/idlfigures
!p.font = 1
!p.thick = 2
!x.thick = 2
!y.thick = 2
!z.thick = 2
xsize = 8.8
xsize = 16
wall = 0.03
margin=0.16
a = xsize/8.8 - (margin + wall)
b = a * 2d / (1 + sqrt(5))

ysize = (margin + b + wall + b + wall)*8.8
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

;  this was to make room for legend  x2=0.8
xdata=[0,1,2,3]
ydata=[2,3,4,5]

fontsize=16
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;   NORMAL COLOR   ;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

SET_PLOT, 'PS'

device, filename=plotfilename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

plot, xdata, ydata, /nodata, /noerase, position=[x1,y2,x2,y2+y2-y1], $
xtitle=' ',   ytitle='uvw1 - v', charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrangew1v, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1)


cgcolorfill, redw1vxpoly, redw1vypoly, color='pink'

cgcolorfill, redwhitew1vxpoly, redwhitew1vypoly, color=cgcolor('white')

cgcolorfill, bluew1vxpoly, bluew1vypoly, color='powder blue'

cgcolorfill, bluewhitew1vxpoly, bluewhitew1vypoly, color=cgcolor('white')



for n=0,n_elements(SNnames)-1 do begin
	SNname=SNnames[n]
	print, SNname

	pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat',   dt=dt


	fitsfilename='~/Desktop/Dropbox/SN/SOUSA/fitting/'+SNname+'_6fits16.sav'

	if file_test(fitsfilename) eq 1 then  begin
		restore, fitsfilename, verbose=0
		Bpeaktime=combbpeakmjd
	endif else begin
		Bpeaktime=dt.bb[0,where(dt.bb[1,*] eq min(dt.bb[1,*],/nan))]
	endelse

;;;;;;


		if SNname eq 'SN2005cf' or SNname eq 'SN2005df' or  SNname eq 'SN2011fe'  or  SNname eq 'SN2013aa'   or  SNname eq 'SN2011iv' then begin

			groundbbfile='$SOUSA/grounddata/'+SNname+'_bb_ground.dat'

			readcol,groundbbfile,bbmjd, bbmag, bbmagerr,/silent

			dt.mag_array[4,*]=interpol(bbmag,bbmjd,dt.time_array)
			dt.magerr_array[4,*]=interpol(bbmagerr,bbmjd,dt.time_array)

			groundvvfile='$SOUSA/grounddata/'+SNname+'_vv_ground.dat'

			readcol,groundvvfile,vvmjd, vvmag, vvmagerr,/silent

			dt.mag_array[5,*]=interpol(vvmag,vvmjd,dt.time_array)
			dt.magerr_array[5,*]=interpol(vvmagerr,vvmjd,dt.time_array)

		endif

;;;;;;;





	if finite(Bpeaktime) eq 0 then Bpeaktime=dt.vv[0,where(dt.vv[1,*] eq min(dt.vv[1,*],/nan))]


	if keyword_set(peakdates) then Bpeaktime=peakdates[n]
;  not sure if this works
	w1v=where( finite(dt.mag_array[2,*]) eq 1 and finite(dt.mag_array[5,*]) eq 1 )

	colortable[n]=128+floor( (min(dt.mag_array[2,w1v]-dt.mag_array[5,w1v])-0.5)/3.0*256)


	oploterror, dt.time_array[w1v]-Bpeaktime[0], dt.mag_array[2,w1v]-dt.mag_array[5,w1v], sqrt(dt.magerr_array[2,w1v]^2.0+dt.magerr_array[5,w1v]^2.0), symsize=0.1, linestyle=1, hatlength=0
	cgoplot, dt.time_array[w1v]-Bpeaktime[0], dt.mag_array[2,w1v]-dt.mag_array[5,w1v], psym=symbols[n], color=colortable[n], symsize=sizes[n], linestyle=1


;if SNname eq 'SN2012dn' then stop

endfor


;;;;;;;;;;;;;;;;;

plot, xdata, ydata, /nodata, /noerase, position=[x1,y2,x2,y2+y2-y1], $
xtitle=' ',   ytitle='uvw1 - v', charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrangew1v, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1)



;axis, xaxis=1, xminor=1, xticklen=xticklen,  xrange=xrange, xstyle=1, xticks=nxticks, xtickv=xtickvalues ;, xtickname=replicate(' ',nxticks+1)


; al_legend, SNnames, psym=symbols, color=colors, position=[0.8,y2+y2-y1], /norm, background='white', charsize=0.8

al_legend, ['NUV-red','NUV-blue'], linestyle=[0,0], thick=20, color=['pink','powder blue'], position=[5,1.2], box=0, charsize=1.0


;al_legend, SNnames,   psym=symbols, symsize=0.5, color=colortable, $
;pos=[30,3.5], charsize=0.75, box=1, background='white'

;;;;;;;         m2 - vv

print, 'moving on to m2-vv'

plot, xdata, ydata, /nodata, /noerase, position=[x1,y1,x2,y2], $
xtitle='Days from Bpeak',   ytitle='uvm2 - v', charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrangem2vv, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues



for n=0,n_elements(SNnames)-1 do begin
	SNname=SNnames[n]
	print, SNname

	pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat',   dt=dt


	fitsfilename='~/Desktop/Dropbox/SN/SOUSA/fitting/'+SNname+'_6fits16.sav'

	if file_test(fitsfilename) eq 1 then  begin
		restore, fitsfilename, verbose=0
		Bpeaktime=combbpeakmjd
	endif else begin
		Bpeaktime=dt.bb[0,where(dt.bb[1,*] eq min(dt.bb[1,*],/nan))]
	endelse

;;;;;;


		if SNname eq 'SN2005cf' or SNname eq 'SN2005df' or  SNname eq 'SN2011fe'  or  SNname eq 'SN2013aa'   or  SNname eq 'SN2011iv' then begin

			groundbbfile='$SOUSA/grounddata/'+SNname+'_bb_ground.dat'

			readcol,groundbbfile,bbmjd, bbmag, bbmagerr,/silent

			dt.mag_array[4,*]=interpol(bbmag,bbmjd,dt.time_array)
			dt.magerr_array[4,*]=interpol(bbmagerr,bbmjd,dt.time_array)

			groundvvfile='$SOUSA/grounddata/'+SNname+'_vv_ground.dat'

			readcol,groundvvfile,vvmjd, vvmag, vvmagerr,/silent

			dt.mag_array[5,*]=interpol(vvmag,vvmjd,dt.time_array)
			dt.magerr_array[5,*]=interpol(vvmagerr,vvmjd,dt.time_array)

		endif

;;;;;;;





	if finite(Bpeaktime) eq 0 then Bpeaktime=dt.vv[0,where(dt.vv[1,*] eq min(dt.vv[1,*],/nan))]


	if keyword_set(peakdates) then Bpeaktime=peakdates[n]

	m2w1=where( finite(dt.mag_array[2,*]) eq 1 and finite(dt.mag_array[1,*]) eq 1 )
	m2vv=where( finite(dt.mag_array[5,*]) eq 1 and finite(dt.mag_array[1,*]) eq 1 )
	oploterror, dt.time_array[m2vv]-Bpeaktime[0], dt.mag_array[1,m2vv]-dt.mag_array[5,m2vv], sqrt(dt.magerr_array[5,m2vv]^2.0+dt.magerr_array[1,m2vv]^2.0), symsize=0.1, linestyle=1, hatlength=0

	cgoplot, dt.time_array[m2vv]-Bpeaktime[0], dt.mag_array[1,m2vv]-dt.mag_array[5,m2vv], psym=symbols[n], color=colortable[n], symsize=sizes[n], linestyle=1


;if SNname eq 'SN2012dn' then print, 'SN2012dn'
;if SNname eq 'SN2012dn' then stop

endfor



;;;;;;;;;;;;;;;;;


device, /close
SET_PLOT, 'X'


spawn, 'open '+plotfilename




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;; set values of constants
H_0=72.0     ;
H_0_err=5.0;
vel_thermal=200.0     ; uncertainty in velocity of galaxies due to thermal/gravitational 
                      ; motion after large scale structure has been removed

SET_PLOT, 'PS'

device, filename=absplotfilename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

cgplot, xdata, ydata, /nodata, /noerase, position=[x1,y1,x2,1.0*(y2+y2-y1-y1)/3.0+y1], $
xtitle='Days from Bpeak',   ytitle='v - $\mu$', charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=vvabsrange, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues




for n=0,n_elements(shortnormallist)-1 do begin
	SNname=shortnormallist[n]
	print, SNname

	pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat',   dt=dt


	fitsfilename='~/Desktop/Dropbox/SN/SOUSA/fitting/'+SNname+'_6fits16.sav'

	if file_test(fitsfilename) eq 1 then  begin
		restore, fitsfilename, verbose=0
		Bpeaktime=combbpeakmjd
	endif else begin
		Bpeaktime=dt.bb[0,where(dt.bb[1,*] eq min(dt.bb[1,*],/nan))]
	endelse

;;;;;;


		if SNname eq 'SN2005cf' or SNname eq 'SN2005df' or  SNname eq 'SN2011fe'  or  SNname eq 'SN2013aa'   or  SNname eq 'SN2011iv' then begin

			groundbbfile='$SOUSA/grounddata/'+SNname+'_bb_ground.dat'

			readcol,groundbbfile,bbmjd, bbmag, bbmagerr,/silent

			dt.mag_array[4,*]=interpol(bbmag,bbmjd,dt.time_array)
			dt.magerr_array[4,*]=interpol(bbmagerr,bbmjd,dt.time_array)

			groundvvfile='$SOUSA/grounddata/'+SNname+'_vv_ground.dat'

			readcol,groundvvfile,vvmjd, vvmag, vvmagerr,/silent

			dt.mag_array[5,*]=interpol(vvmag,vvmjd,dt.time_array)
			dt.magerr_array[5,*]=interpol(vvmagerr,vvmjd,dt.time_array)

		endif

;;;;;;;





	if finite(Bpeaktime) eq 0 then Bpeaktime=dt.vv[0,where(dt.vv[1,*] eq min(dt.vv[1,*],/nan))]


	if keyword_set(peakdates) then Bpeaktime=peakdates[n]

index=where(host.snname_array eq SNname)


	mu_best=host.dm_best_array[index]
	mu_best_err=host.dm_best_err_array[index]
	mu_best_ref=host.dm_best_ref_array[index]


	print, SNname, ' mu ', mu_best

vvgood=where( finite(dt.vv[1,*]) eq 1 )

	cgoplot, dt.vv[0,vvgood]-Bpeaktime[0], dt.vv[1,vvgood]-mu_best[0], psym=0, color=shortcolors[n], symsize=0, thick=10

endfor


for n=0,n_elements(shortnormallist)-1 do begin
	SNname=shortnormallist[n]
	print, SNname
index=where(host.snname_array eq SNname)

	pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat',   dt=dt

	fitsfilename='~/Desktop/Dropbox/SN/SOUSA/fitting/'+SNname+'_6fits16.sav'

	if file_test(fitsfilename) eq 1 then  begin
		restore, fitsfilename, verbose=0
		Bpeaktime=combbpeakmjd
	endif else begin
		Bpeaktime=dt.bb[0,where(dt.bb[1,*] eq min(dt.bb[1,*],/nan))]
	endelse
	if finite(Bpeaktime) eq 0 then Bpeaktime=dt.vv[0,where(dt.vv[1,*] eq min(dt.vv[1,*],/nan))]

	if keyword_set(peakdates) then Bpeaktime=peakdates[n]

	index=where(host.snname_array eq SNname)

	mu_best=host.dm_best_array[index]
	mu_best_err=host.dm_best_err_array[index]
	mu_best_ref=host.dm_best_ref_array[index]



vvgood=where( finite(dt.vv[1,*]) eq 1 )

oploterror, dt.vv[0,vvgood]-Bpeaktime[0], dt.vv[1,vvgood]-mu_best[0], dt.vv[2,vvgood], linestyle=2, psym=-symcat(symbols[n]), symsize=sizes[0], color=colors[n], errcolor=cgcolor(colors[n])

	cgoplot, dt.vv[0,vvgood]-Bpeaktime[0], dt.vv[1,vvgood]-mu_best[0], psym=symbols[n], color=colors[n], symsize=1

endfor


al_legend, shortnormallist[0:1],   psym=0, symsize=1.0, color=shortcolors[0:1], thick=10, $
pos=[0.25,0.5], /norm, charsize=0.8, box=0
al_legend, shortnormallist[2:3],   psym=0, symsize=1.0, color=shortcolors[2:3], thick=10, $
pos=[0.25,0.24], /norm, charsize=0.8, box=0



cgplot, xdata, ydata,  /nodata, /noerase, position=[x1,1.0*(y2+y2-y1-y1)/3.0+y1,x2,2.0*(y2+y2-y1-y1)/3.0+y1], $
 ytitle='uvw1 - $\mu$', charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=w1absrange, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1)




for n=0,n_elements(shortnormallist)-1 do begin
	SNname=shortnormallist[n]
	print, SNname

	pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat',   dt=dt


	fitsfilename='~/Desktop/Dropbox/SN/SOUSA/fitting/'+SNname+'_6fits16.sav'

	if file_test(fitsfilename) eq 1 then  begin
		restore, fitsfilename, verbose=0
		Bpeaktime=combbpeakmjd
	endif else begin
		Bpeaktime=dt.bb[0,where(dt.bb[1,*] eq min(dt.bb[1,*],/nan))]
	endelse

;;;;;;


		if SNname eq 'SN2005cf' or SNname eq 'SN2005df' or  SNname eq 'SN2011fe'  or  SNname eq 'SN2013aa'   or  SNname eq 'SN2011iv' then begin

			groundbbfile='$SOUSA/grounddata/'+SNname+'_bb_ground.dat'

			readcol,groundbbfile,bbmjd, bbmag, bbmagerr,/silent

			dt.mag_array[4,*]=interpol(bbmag,bbmjd,dt.time_array)
			dt.magerr_array[4,*]=interpol(bbmagerr,bbmjd,dt.time_array)

			groundvvfile='$SOUSA/grounddata/'+SNname+'_vv_ground.dat'

			readcol,groundvvfile,vvmjd, vvmag, vvmagerr,/silent

			dt.mag_array[5,*]=interpol(vvmag,vvmjd,dt.time_array)
			dt.magerr_array[5,*]=interpol(vvmagerr,vvmjd,dt.time_array)

		endif

;;;;;;;



	if finite(Bpeaktime) eq 0 then Bpeaktime=dt.vv[0,where(dt.vv[1,*] eq min(dt.vv[1,*],/nan))]


	if keyword_set(peakdates) then Bpeaktime=peakdates[n]

index=where(host.snname_array eq SNname)

	mu_best=host.dm_best_array[index]
	mu_best_err=host.dm_best_err_array[index]
	mu_best_ref=host.dm_best_ref_array[index]

w1good=where( finite(dt.w1[1,*]) eq 1 )

	cgoplot, dt.w1[0,w1good]-Bpeaktime[0], dt.w1[1,w1good]-mu_best[0], psym=0, color=shortcolors[n], symsize=0, thick=10

endfor


for n=0,n_elements(SNnames)-1 do begin
	SNname=SNnames[n]
	print, SNname
index=where(host.snname_array eq SNname)

	pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat',   dt=dt

	fitsfilename='~/Desktop/Dropbox/SN/SOUSA/fitting/'+SNname+'_6fits16.sav'

	if file_test(fitsfilename) eq 1 then  begin
		restore, fitsfilename, verbose=0
		Bpeaktime=combbpeakmjd
	endif else begin
		Bpeaktime=dt.bb[0,where(dt.bb[1,*] eq min(dt.bb[1,*],/nan))]
	endelse
	if finite(Bpeaktime) eq 0 then Bpeaktime=dt.vv[0,where(dt.vv[1,*] eq min(dt.vv[1,*],/nan))]

	if keyword_set(peakdates) then Bpeaktime=peakdates[n]

	index=where(host.snname_array eq SNname)



	mu_best=host.dm_best_array[index]
	mu_best_err=host.dm_best_err_array[index]
	mu_best_ref=host.dm_best_ref_array[index]

w1good=where( finite(dt.w1[1,*]) eq 1 )

oploterror, dt.w1[0,w1good]-Bpeaktime[0], dt.w1[1,w1good]-mu_best[0], dt.w1[2,w1good], linestyle=2, psym=-symcat(symbols[n]), symsize=sizes[0], color=colors[n], errcolor=cgcolor(colors[n])

	cgoplot, dt.w1[0,w1good]-Bpeaktime[0], dt.w1[1,w1good]-mu_best[0], psym=symbols[n], color=colors[n], symsize=1

endfor
print, 'moving on to m2'

cgplot, xdata, ydata, /nodata, /noerase, position=[x1,2.0*(y2+y2-y1-y1)/3.0+y1,x2,3.0*(y2+y2-y1-y1)/3.0+y1],  $      
ytitle='uvm2 - $\mu$', charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=m2absrange, ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1)





for n=0,n_elements(shortnormallist)-1 do begin
	SNname=shortnormallist[n]
	print, SNname

	pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat',   dt=dt


	fitsfilename='~/Desktop/Dropbox/SN/SOUSA/fitting/'+SNname+'_6fits16.sav'

	if file_test(fitsfilename) eq 1 then  begin
		restore, fitsfilename, verbose=0
		Bpeaktime=combbpeakmjd
	endif else begin
		Bpeaktime=dt.bb[0,where(dt.bb[1,*] eq min(dt.bb[1,*],/nan))]
	endelse
	Bpeaktime=Bpeaktime[0]
;;;;;;


		if SNname eq 'SN2005cf' or SNname eq 'SN2005df' or  SNname eq 'SN2011fe'  or  SNname eq 'SN2013aa'   or  SNname eq 'SN2011iv' then begin

			groundbbfile='$SOUSA/grounddata/'+SNname+'_bb_ground.dat'

			readcol,groundbbfile,bbmjd, bbmag, bbmagerr,/silent

			dt.mag_array[4,*]=interpol(bbmag,bbmjd,dt.time_array)
			dt.magerr_array[4,*]=interpol(bbmagerr,bbmjd,dt.time_array)

			groundvvfile='$SOUSA/grounddata/'+SNname+'_vv_ground.dat'

			readcol,groundvvfile,vvmjd, vvmag, vvmagerr,/silent

			dt.mag_array[5,*]=interpol(vvmag,vvmjd,dt.time_array)
			dt.magerr_array[5,*]=interpol(vvmagerr,vvmjd,dt.time_array)

		endif

;;;;;;;





	if finite(Bpeaktime) eq 0 then Bpeaktime=dt.vv[0,where(dt.vv[1,*] eq min(dt.vv[1,*],/nan))]


	if keyword_set(peakdates) then Bpeaktime=peakdates[n]

	index=where(host.snname_array eq SNname)


	mu_best=host.dm_best_array[index]
	mu_best_err=host.dm_best_err_array[index]
	mu_best_ref=host.dm_best_ref_array[index]

m2good=where( finite(dt.m2[1,*]) eq 1 )

	cgoplot, dt.m2[0,m2good]-Bpeaktime[0], dt.m2[1,m2good]-mu_best[0], psym=0, color=shortcolors[n], symsize=0, thick=10

endfor




for n=0,n_elements(SNnames)-1 do begin
	SNname=SNnames[n]
	print, SNname

	pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat',   dt=dt


	fitsfilename='~/Desktop/Dropbox/SN/SOUSA/fitting/'+SNname+'_6fits16.sav'

	if file_test(fitsfilename) eq 1 then  begin
		restore, fitsfilename, verbose=0
		Bpeaktime=combbpeakmjd
	endif else begin
		Bpeaktime=dt.bb[0,where(dt.bb[1,*] eq min(dt.bb[1,*],/nan))]
	endelse

;;;;;;


		if SNname eq 'SN2005cf' or SNname eq 'SN2005df' or  SNname eq 'SN2011fe'  or  SNname eq 'SN2013aa'   or  SNname eq 'SN2011iv' then begin

			groundbbfile='$SOUSA/grounddata/'+SNname+'_bb_ground.dat'

			readcol,groundbbfile,bbmjd, bbmag, bbmagerr,/silent

			dt.mag_array[4,*]=interpol(bbmag,bbmjd,dt.time_array)
			dt.magerr_array[4,*]=interpol(bbmagerr,bbmjd,dt.time_array)

			groundvvfile='$SOUSA/grounddata/'+SNname+'_vv_ground.dat'

			readcol,groundvvfile,vvmjd, vvmag, vvmagerr,/silent

			dt.mag_array[5,*]=interpol(vvmag,vvmjd,dt.time_array)
			dt.magerr_array[5,*]=interpol(vvmagerr,vvmjd,dt.time_array)

		endif

;;;;;;;





	if finite(Bpeaktime) eq 0 then Bpeaktime=dt.vv[0,where(dt.vv[1,*] eq min(dt.vv[1,*],/nan))]


	if keyword_set(peakdates) then Bpeaktime=peakdates[n]

	index=where(host.snname_array eq SNname)


	mu_best=host.dm_best_array[index]
	mu_best_err=host.dm_best_err_array[index]
	mu_best_ref=host.dm_best_ref_array[index]

m2good=where( finite(dt.m2[1,*]) eq 1 )

oploterror, dt.m2[0,m2good]-Bpeaktime[0], dt.m2[1,m2good]-mu_best[0], dt.m2[2,m2good], linestyle=2, psym=-symcat(symbols[n]), symsize=sizes[0], color=colors[n], errcolor=cgcolor(colors[n])

	cgoplot, dt.m2[0,m2good]-Bpeaktime[0], dt.m2[1,m2good]-mu_best[0], psym=symbols[n], color=colors[n], symsize=1

endfor

;al_legend, SNnames,   psym=symbols, symsize=1.0, color=colors, $
;pos=[0.68,0.6], /norm, charsize=0.8, box=1, background='white'

device, /close
SET_PLOT, 'X'



print, 'final stop'

stop
end

