pro redfilterplots

EBV=1.2

SN2011fe='SN2011fe_110913_full.dat'
readcol,SN2011fe,wave_SN2011fe,flux_SN2011fe,/silent

pjb_uvotspec_all, SN2011fe, counts_array=counts_array, lambda=lambda, filter_array=filter_array, countspec_array=countspec_array,  spectrum_uvotres=spectrum_uvotres


SN2011fe_spectrum_uvotres=spectrum_uvotres
SN2011fe_countspec_array=countspec_array

Amw=sne_mw_reddening(lambda,Ebv)
specredmw=SN2011fe_spectrum_uvotres*10.0^(-Amw/2.5)


pjb_uvotspec_all, [transpose(lambda),transpose(specredmw)], mag_array=mag_array, counts_array=counts_array, lambda=lambda, filter_array=filter_array, countspec_array=countspec_array,  spectrum_uvotres=spectrum_uvotres
SN2011fe_mwred_countspec_array=countspec_array
SN2011fe_mwred_spectrum_uvotres=spectrum_uvotres
SN2011fe_mwred_mag_array=mag_array


Acslmc=sne_goobarlmc_reddening(lambda,Ebv)
specredcslmc=SN2011fe_spectrum_uvotres*10.0^(-Acslmc/2.5)


pjb_uvotspec_all, [transpose(lambda),transpose(specredcslmc)], mag_array=mag_array, counts_array=counts_array, lambda=lambda, filter_array=filter_array, countspec_array=countspec_array,  spectrum_uvotres=spectrum_uvotres
SN2011fe_cslmcred_countspec_array=countspec_array
SN2011fe_cslmcred_mag_array=mag_array
SN2011fe_cslmcred_spectrum_uvotres=spectrum_uvotres



;

nxplots=1
nyplots=9
nplots=nxplots*nyplots

; from http://www.iluvatar.org/~dwijn/idlfigures
!p.font = 1
!p.thick = 2
!x.thick = 2
!y.thick = 2
!z.thick = 2
; the default size is given in centimeters
; 8.8 is made to match a single journal column width
xsize = 18.0
;xsize = 8.8

; putting these in real size (centimeters) 
; so they are the same regardless of the number of plots
wall = 0.03*8.8 ; these values were originally made for the 8.8 width
margin=0.16*8.8 ; a little big if the y-axis label is single digits but accomodates more digits
; width and height of single plots
a = (xsize - margin - wall)/nxplots
; for square plots
; b = a 
; for golden ratio plots
b = a * 2d / (1 + sqrt(5))
b = 0.5 * a * 2d / (1 + sqrt(5))

ysize= b * nyplots + margin + wall

ysize = (margin + nyplots*(b + wall ) )
ysize=30.0
b = (ysize-margin)/nyplots - wall 

ticklen = 0.01
xticklen = ticklen/b
yticklen = ticklen/a

;ymax=40
;nyticks=7


ytitle0='Normalized Counts'
ytitle1='Normalized Counts'
ytitle2='Normalized Counts'
ytitle3='Normalized Counts'
ytitle4='Normalized Counts'
ytitle5='Normalized Counts'
ytitle6='log (Effective Area)'
ytitle7='A_lambda/E(B-V)'
ytitle8=' Flux Density *10^16'


xrange0=[1500,8000]
nxticks=13


xtitle0='Wavelength [Angstroms]'


fontsize=12

xdata=[0,1,2,3,4]
ydata=[2,4,6,8,10]


figurename='redfilterplots.eps'

SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize, bits_per_pixel=8, /color


yrange8=[0.0000001,100000.0]


nx=0
ny=8
plot, xdata,ydata,  /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle=' ',   ytitle=ytitle8, charsize=1.0,  $
xminor=1, yminor=1, /ylog, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange8, ystyle=1, xrange=xrange0, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)


;cgoplot,  wave_SN2011fe,flux_SN2011fe*10.0^(16.0), color='black'
cgoplot,  lambda, SN2011fe_spectrum_uvotres*10.0^(16.0), color='black'

cgoplot,  lambda, SN2011fe_mwred_spectrum_uvotres*10.0^(16.0), color='red', linestyle=2
cgoplot,  lambda, SN2011fe_cslmcred_spectrum_uvotres*10.0^(16.0), color='dark green', linestyle=3

al_legend, ['unreddened SN 2011fe spectrum','MW reddened SN 2011fe spectrum','CSLMC reddened SN 2011fe spectrum'],   linestyle=[0,2,3], color=[cgcolor('black'),cgcolor('red'),cgcolor('dark green')], background_color='white', $
pos=[0.5,0.88], /norm, charsize=1, box=0



yrange7=[0,24.0]


nx=0
ny=7
plot, xdata,ydata,  /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle=' ',   ytitle=ytitle7, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange7, ystyle=1, xrange=xrange0, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

EBV=1.2
rv=3.1
	lambda_av_mw31=sne_mw_reddening(lambda,Ebv,rv=rv)
rv=1.7	
	lambda_av_cslmc=sne_goobarlmc_reddening(lambda,Ebv)

cgoplot, lambda, 	lambda_av_mw31, color='red', linestyle=2

cgoplot, lambda, 	lambda_av_cslmc, color='dark green', linestyle=3

al_legend, ['CSLMC R_V=1.7','MW R_V=3.1'],   linestyle=[3,2],  color=[cgcolor('dark green'),cgcolor('red')], background_color='white', $
pos=[0.5,0.81], /norm, charsize=1, box=0

nx=0
ny=6

yrange6=[-5,1]

plot, xdata,ydata,  /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle=' ',   ytitle=ytitle6, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange6, ystyle=1, xrange=xrange0, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

filters=['uvw2','uvm2','uvw1','u','b','v']
colors=['purple','red','black','violet','blue','dark green']
for f=0,5 do cgoplot,  lambda, alog10(filter_array[f,*]/max(filter_array[f,*])), color=colors[f]



al_legend, filters,   linestyle=0,  color=cgcolor(colors), background_color='white', $
pos=[0.8,0.72], /norm, charsize=1, box=0



nx=0
ny=5
f=0
yrange5=[0,max(filter_array[f,*]/total(filter_array[f,*]))*1.2]
yrange5=[0,0.02]

plot, xdata,ydata,  /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle=' ',   ytitle=ytitle5, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange5, ystyle=1, xrange=xrange0, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=2, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

f=0

cgoplot,  lambda, filter_array[f,*]/total(filter_array[f,*]), color='black', linestyle=0
cgoplot,  lambda, SN2011fe_countspec_array[f,*]/total(SN2011fe_countspec_array[f,*]), color='blue', linestyle=1
cgoplot,  lambda, SN2011fe_mwred_countspec_array[f,*]/total(SN2011fe_mwred_countspec_array[f,*]), color='red', linestyle=2
cgoplot,  lambda, SN2011fe_cslmcred_countspec_array[f,*]/total(SN2011fe_cslmcred_countspec_array[f,*]), color='dark green', linestyle=3


al_legend, ['UVOT uvw2 filter','SN2011fe unred','SN2011fe MW', 'SN2011fe CSLMC'],   linestyle=[0,1,2,3],  color=[cgcolor('black'),cgcolor('blue'),cgcolor('red'),cgcolor('dark green')], background_color='white', $
pos=[0.6,0.62], /norm, charsize=1, box=0



nx=0
ny=4
f=1
yrange4=[0,max(filter_array[f,*]/total(filter_array[f,*]))*1.2]
yrange4=[0,0.02]
plot, xdata,ydata,  /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle=' ',   ytitle=ytitle4, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange4, ystyle=1, xrange=xrange0, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=2, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

f=1

cgoplot,  lambda, filter_array[f,*]/total(filter_array[f,*]), color='black', linestyle=0
cgoplot,  lambda, SN2011fe_countspec_array[f,*]/total(SN2011fe_countspec_array[f,*]), color='blue', linestyle=1
cgoplot,  lambda, SN2011fe_mwred_countspec_array[f,*]/total(SN2011fe_mwred_countspec_array[f,*]), color='red', linestyle=2
cgoplot,  lambda, SN2011fe_cslmcred_countspec_array[f,*]/total(SN2011fe_cslmcred_countspec_array[f,*]), color='dark green', linestyle=3


al_legend, ['UVOT uvm2 filter','SN2011fe unred','SN2011fe MW', 'SN2011fe CSLMC'],   linestyle=[0,1,2,3],  color=[cgcolor('black'),cgcolor('blue'),cgcolor('red'),cgcolor('dark green')], background_color='white', $
pos=[0.6,0.52], /norm, charsize=1, box=0



nx=0
ny=3
f=2
yrange3=[0,max(SN2011fe_countspec_array[f,*]/total(SN2011fe_countspec_array[f,*]))*1.2]
yrange3=[0,0.02]
plot, xdata,ydata,  /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle=' ',   ytitle=ytitle3, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange3, ystyle=1, xrange=xrange0, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=2, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

f=2

cgoplot,  lambda, filter_array[f,*]/total(filter_array[f,*]), color='black', linestyle=0
cgoplot,  lambda, SN2011fe_countspec_array[f,*]/total(SN2011fe_countspec_array[f,*]), color='blue', linestyle=1
cgoplot,  lambda, SN2011fe_mwred_countspec_array[f,*]/total(SN2011fe_mwred_countspec_array[f,*]), color='red', linestyle=2
cgoplot,  lambda, SN2011fe_cslmcred_countspec_array[f,*]/total(SN2011fe_cslmcred_countspec_array[f,*]), color='dark green', linestyle=3


al_legend, ['UVOT uvw1 filter','SN2011fe unred','SN2011fe MW', 'SN2011fe CSLMC'],   linestyle=[0,1,2,3],  color=[cgcolor('black'),cgcolor('blue'),cgcolor('red'),cgcolor('dark green')], background_color='white', $
pos=[0.6,0.42], /norm, charsize=1, box=0



nx=0
ny=2
f=3
yrange2=[0,max(SN2011fe_cslmcred_countspec_array[f,*]/total(SN2011fe_cslmcred_countspec_array[f,*]))*1.2]
plot, xdata,ydata,  /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle=' ',   ytitle=ytitle2, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange2, ystyle=1, xrange=xrange0, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

f=3

cgoplot,  lambda, filter_array[f,*]/total(filter_array[f,*]), color='black', linestyle=0
cgoplot,  lambda, SN2011fe_countspec_array[f,*]/total(SN2011fe_countspec_array[f,*]), color='blue', linestyle=1
cgoplot,  lambda, SN2011fe_mwred_countspec_array[f,*]/total(SN2011fe_mwred_countspec_array[f,*]), color='red', linestyle=2
cgoplot,  lambda, SN2011fe_cslmcred_countspec_array[f,*]/total(SN2011fe_cslmcred_countspec_array[f,*]), color='dark green', linestyle=3


al_legend, ['UVOT u filter','SN2011fe unred','SN2011fe MW', 'SN2011fe CSLMC'],   linestyle=[0,1,2,3],  color=[cgcolor('black'),cgcolor('blue'),cgcolor('red'),cgcolor('dark green')], background_color='white', $
pos=[0.6,0.31], /norm, charsize=1, box=0



nx=0
ny=1
f=4
yrange1=[0,max(SN2011fe_cslmcred_countspec_array[f,*]/total(SN2011fe_cslmcred_countspec_array[f,*]))*1.2]
yrange1=[0,0.02]
plot, xdata,ydata,  /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle=' ',   ytitle=ytitle1, charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yrange1, ystyle=1, xrange=xrange0, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=2, ytickv=ytickvalues , xtickname=replicate(' ',nxticks+1)

f=4

cgoplot,  lambda, filter_array[f,*]/total(filter_array[f,*]), color='black', linestyle=0
cgoplot,  lambda, SN2011fe_countspec_array[f,*]/total(SN2011fe_countspec_array[f,*]), color='blue', linestyle=1
cgoplot,  lambda, SN2011fe_mwred_countspec_array[f,*]/total(SN2011fe_mwred_countspec_array[f,*]), color='red', linestyle=2
cgoplot,  lambda, SN2011fe_cslmcred_countspec_array[f,*]/total(SN2011fe_cslmcred_countspec_array[f,*]), color='dark green', linestyle=3


al_legend, ['UVOT b filter','SN2011fe unred','SN2011fe MW', 'SN2011fe CSLMC'],   linestyle=[0,1,2,3],  color=[cgcolor('black'),cgcolor('blue'),cgcolor('red'),cgcolor('dark green')], background_color='white', $
pos=[0.1,0.21], /norm, charsize=1, box=0



nx=0
ny=0
f=5
yrange0=[0,max(SN2011fe_countspec_array[f,*]/total(SN2011fe_countspec_array[f,*]))*1.2]


plot, xdata,ydata,  /nodata, /noerase, $
position=[(margin+nx*a)/xsize,(margin+b*ny)/ysize,(margin+(nx+1)*a)/xsize,(margin+b*(ny+1))/ysize], $
xtitle=xtitle0, ytitle=ytitle0, charsize=1.0,  $
 yrange=[0,0.02], ystyle=1, xrange=xrange0, xstyle=1, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xticks=xticks, xtickv=xtickvalues, yticks=2, ytickv=ytickvalues


cgoplot,  lambda, filter_array[f,*]/total(filter_array[f,*]), color='black', linestyle=0
cgoplot,  lambda, SN2011fe_countspec_array[f,*]/total(SN2011fe_countspec_array[f,*]), color='blue', linestyle=1
cgoplot,  lambda, SN2011fe_mwred_countspec_array[f,*]/total(SN2011fe_mwred_countspec_array[f,*]), color='red', linestyle=2
cgoplot,  lambda, SN2011fe_cslmcred_countspec_array[f,*]/total(SN2011fe_cslmcred_countspec_array[f,*]), color='dark green', linestyle=3



al_legend, ['UVOT v filter','SN2011fe unred','SN2011fe MW', 'SN2011fe CSLMC'],   linestyle=[0,1,2,3],  color=[cgcolor('black'),cgcolor('blue'),cgcolor('red'),cgcolor('dark green')], background_color='white', $
pos=[0.1,0.12], /norm, charsize=1, box=0




device, /close
SET_PLOT, 'X'

$open redfilterplots.eps



stop
end

