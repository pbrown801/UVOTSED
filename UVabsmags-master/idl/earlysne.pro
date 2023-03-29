pro earlysne



color=['royal blue', 'dark green', 'red', 'orange', 'pink', 'maroon', 'black', 'dodger blue', 'aqua', 'brown', 'yellow', 'royal blue', 'green', 'red', 'orange', 'pink', 'maroon', 'black', 'dodger blue', 'aqua', 'brown', 'yellow','royal blue', 'dark green', 'red', 'orange', 'pink', 'maroon', 'black', 'dodger blue', 'aqua', 'brown', 'yellow', 'royal blue', 'green', 'red', 'orange', 'pink', 'maroon', 'black', 'dodger blue', 'aqua', 'brown', 'yellow','royal blue', 'dark green', 'red', 'orange', 'pink', 'maroon', 'black', 'dodger blue', 'aqua', 'brown', 'yellow', 'royal blue', 'green', 'red', 'orange', 'pink', 'maroon', 'black', 'dodger blue', 'aqua', 'brown', 'yellow','royal blue', 'dark green', 'red', 'orange', 'pink', 'maroon', 'black', 'dodger blue', 'aqua', 'brown', 'yellow', 'royal blue', 'green', 'red', 'orange', 'pink', 'maroon', 'black', 'dodger blue', 'aqua', 'brown', 'yellow','royal blue', 'dark green', 'red', 'orange', 'pink', 'maroon', 'black', 'dodger blue', 'aqua', 'brown', 'yellow', 'royal blue', 'green', 'red', 'orange', 'pink', 'maroon', 'black', 'dodger blue', 'aqua', 'brown', 'yellow','royal blue', 'dark green', 'red', 'orange', 'pink', 'maroon', 'black', 'dodger blue', 'aqua', 'brown', 'yellow', 'royal blue', 'green', 'red', 'orange', 'pink', 'maroon', 'black', 'dodger blue', 'aqua', 'brown', 'yellow','royal blue', 'dark green', 'red', 'orange', 'pink', 'maroon', 'black', 'dodger blue', 'aqua', 'brown', 'yellow', 'royal blue', 'green', 'red', 'orange', 'pink', 'maroon', 'black', 'dodger blue', 'aqua', 'brown', 'yellow','royal blue', 'dark green', 'red', 'orange', 'pink', 'maroon', 'black', 'dodger blue', 'aqua', 'brown', 'yellow', 'royal blue', 'green', 'red', 'orange', 'pink', 'maroon', 'black', 'dodger blue', 'aqua', 'brown', 'yellow']

; SN2011fe removed because it doesn't have v band data

SNlist=['iPTF14atg',  'SN2012cg', 'SN2012ij','iPTF14bdn', 'iPTF13dge', 'SN2009ig', 'SN2013gy', 'SN2012cp', 'SN2015F', 'SN2007af', 'iPTF13ebh', 'SN2008hv', 'SN2011by', 'SN2012ht', 'SN2011ao', 'SN2011B', 'SN2013dy', 'SN2013gs', 'SN2010Y', 'SN2008ec', 'SN2005ke']



readcol, '../snnamesquality.txt', longlistnames, longlistquality, format='A,A'

SNlist=longlistnames[where(longlistquality eq 'good' or longlistquality eq 'risepeak')]
nprimesne=6

SNlist=['SN2019yvq','SN2020hvf','iPTF14atg', 'SN2012cg', 'SN2012ij','iPTF14bdn',SNlist]

symbol=[15,17,14,15,17,14,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6,4,5,6]

symbol=symbol[0:n_elements(SNlist)-1]
color=color[0:n_elements(SNlist)-1]

restore, '/Users/pbrown/Desktop/SN/github/SOUSA/www/host.sav'

nplots=1
nyplots=1
; from http://www.iluvatar.org/~dwijn/idlfigures
!p.font = 1
!p.thick = 2
!x.thick = 2
!y.thick = 2
!z.thick = 2
; the default size is given in centimeters
; 8.8 is made to match a journal column width
xsize = 8.8
wall = 0.03
margin=0.16
a = xsize/8.8 - (margin + wall)
b = a * 2d / (1 + sqrt(5))

ysize = (margin + nplots*(b + wall ) )*xsize
ticklen = 0.01
xticklen = ticklen/b
yticklen = ticklen/a

x1 = margin*8.8/xsize
x2 = x1 + a*8.8/xsize
xc = x2 + wall*8.8/xsize
y1 = margin*8.8/ysize
y2 = y1 + b*8.8/ysize
xrange=[-5,45]


nxticks=10


xdata=[2,4,6]
ydata=[6,8,9]


fontsize=12

SET_PLOT, 'PS'

device, filename='earlysne_w1.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize, bits_per_pixel=8, /color

x=0
plot, xdata,ydata, /nodata, /noerase, position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle=' Days from Maximum Light',   ytitle='uvw1 mag - uvw1 max', charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xrange=[-15,15],yrange=[6,-1], ystyle=1, xstyle=1

for n=0,n_elements(snlist)-1 do begin

SNname=SNlist[n]

pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat', dt=dt

SNindex=where(host.SNname_array eq SNname)
pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat', dt=dt

good=where(finite(dt.w1[1,*]) eq 1)
maxindex=where(dt.w1[1,good] eq min(dt.w1[1,good]))
w1maxtime=dt.w1[0,good[maxindex]]
w1maxmag=dt.w1[1,good[maxindex]]
maxindex=where(dt.w1[1,good] eq min(dt.w1[1,good]))
w1maxtime=host.peaktime_array[SNindex,2]
w1maxmag=host.appmag_array[SNindex,2]


if SNname eq 'SN2012ij' then w1maxtime=dt.w1[0,good[maxindex+1]]
if SNname eq 'SN2012ij' then w1maxmag=dt.w1[1,good[maxindex+1]]



oploterror, dt.w1[0,where(finite(dt.w1[1,*]) eq 1)]-w1maxtime[0],dt.w1[1,where(finite(dt.w1[1,*]) eq 1)]-w1maxmag[0],dt.w1[2,where(finite(dt.w1[1,*]) eq 1)], psym=-symbol[n] , color=color[n], errcolor=color[n], linestyle=0, symsize=0.5


endfor

al_legend, snlist, color=color[0:n_elements(SNlist)-1], psym=symbol[0:n_elements(SNlist)-1], position=[10,-1], box=1, charsize=0.6, linestyle=1, background='white', symsize=0.5



device, /close
SET_PLOT, 'X'
$open earlysne_w1.eps




SET_PLOT, 'PS'

device, filename='earlysne_v.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize, bits_per_pixel=8, /color

x=0
plot, xdata,ydata, /nodata, /noerase, position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle=' Days from Maximum Light',   ytitle='v mag - v max', charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xrange=[-20,10],yrange=[6,-1], ystyle=1, xstyle=1

for n=0,n_elements(snlist)-1 do begin

SNname=SNlist[n]

pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat', dt=dt

good=where(finite(dt.vv[1,*]) eq 1)
maxindex=where(dt.vv[1,good] eq min(dt.vv[1,good]))
vvmaxtime=dt.vv[0,good[maxindex]]
vvmaxmag=dt.vv[1,good[maxindex]]

oploterror, dt.vv[0,where(finite(dt.vv[1,*]) eq 1)]-vvmaxtime[0],dt.vv[1,where(finite(dt.vv[1,*]) eq 1)]-vvmaxmag[0],dt.vv[2,where(finite(dt.vv[1,*]) eq 1)], psym=-symbol[n] , color=color[n], errcolor=color[n], linestyle=0, symsize=0.5


endfor

al_legend, snlist, color=color[0:n_elements(SNlist)-1], psym=symbol[0:n_elements(SNlist)-1], position=[5,-1], box=1, charsize=0.6, linestyle=1, background='white', symsize=0.5



device, /close
SET_PLOT, 'X'
$open earlysne_v.eps

;;;;;;;; switch others to grey


color=['royal blue','black', 'dark green', 'red', 'purple', 'maroon','grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey','grey','grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey','grey','grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey','grey','grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey','grey','grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey','grey','grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey','grey','grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey','grey','grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey','grey','grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey','grey','grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey','grey','grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey', 'grey','grey']

color=color[0:n_elements(SNlist)-1]

SET_PLOT, 'PS'

device, filename='earlysne_w1_grey.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize, bits_per_pixel=8, /color

x=0
plot, xdata,ydata, /nodata, /noerase, position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle=' Days from Maximum Light',   ytitle='uvw1 mag - uvw1 max', charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xrange=[-15,15],yrange=[6,-1], ystyle=1, xstyle=1

for n=0,n_elements(snlist)-1 do begin

SNname=SNlist[n]

SNindex=where(host.SNname_array eq SNname)
pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat', dt=dt

good=where(finite(dt.w1[1,*]) eq 1)
maxindex=where(dt.w1[1,good] eq min(dt.w1[1,good]))
w1maxtime=dt.w1[0,good[maxindex]]
w1maxmag=dt.w1[1,good[maxindex]]
maxindex=where(dt.w1[1,good] eq min(dt.w1[1,good]))
w1maxtime=host.peaktime_array[SNindex,2]
w1maxmag=host.appmag_array[SNindex,2]



if SNname eq 'SN2012ij' then w1maxtime=dt.w1[0,good[maxindex+1]]
if SNname eq 'SN2012ij' then w1maxmag=dt.w1[1,good[maxindex+1]]

if SNname eq 'SN2019yvq' then w1maxtime=dt.w1[0,8]
if SNname eq 'SN2019yvq' then w1maxmag=dt.w1[1,8]


oploterror, dt.w1[0,where(finite(dt.w1[1,*]) eq 1)]-w1maxtime[0],dt.w1[1,where(finite(dt.w1[1,*]) eq 1)]-w1maxmag[0],dt.w1[2,where(finite(dt.w1[1,*]) eq 1)], psym=-symbol[n] , color=color[n], errcolor=color[n], linestyle=0, symsize=0.5


endfor


for n=0,3 do begin

SNname=SNlist[n]

pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat', dt=dt

SNindex=where(host.SNname_array eq SNname)
pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat', dt=dt

good=where(finite(dt.w1[1,*]) eq 1)
maxindex=where(dt.w1[1,good] eq min(dt.w1[1,good]))
w1maxtime=dt.w1[0,good[maxindex]]
w1maxmag=dt.w1[1,good[maxindex]]
maxindex=where(dt.w1[1,good] eq min(dt.w1[1,good]))
w1maxtime=host.peaktime_array[SNindex,2]
w1maxmag=host.appmag_array[SNindex,2]



if SNname eq 'SN2012ij' then w1maxtime=dt.w1[0,good[maxindex+1]]
if SNname eq 'SN2012ij' then w1maxmag=dt.w1[1,good[maxindex+1]]


oploterror, dt.w1[0,where(finite(dt.w1[1,*]) eq 1)]-w1maxtime[0],dt.w1[1,where(finite(dt.w1[1,*]) eq 1)]-w1maxmag[0],dt.w1[2,where(finite(dt.w1[1,*]) eq 1)], psym=-symbol[n] , color=color[n], errcolor=color[n], linestyle=0, symsize=0.8


endfor


al_legend, snlist, color=color[0:n_elements(SNlist)-1], psym=symbol[0:n_elements(SNlist)-1], position=[10,-1], box=1, charsize=0.6, linestyle=1, background='white', symsize=0.5



device, /close
SET_PLOT, 'X'
$open earlysne_w1_grey.eps





SET_PLOT, 'PS'

device, filename='earlysne_v_grey.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize, bits_per_pixel=8, /color

x=0
plot, xdata,ydata, /nodata, /noerase, position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle=' Days from Maximum Light',   ytitle='v mag - v max', charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xrange=[-20,10],yrange=[6,-1], ystyle=1, xstyle=1

for n=0,n_elements(snlist)-1 do begin

	SNname=SNlist[n]

	pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat', dt=dt

	good=where(finite(dt.vv[1,*]) eq 1)
	maxindex=where(dt.vv[1,good] eq min(dt.vv[1,good]))
	vvmaxtime=dt.vv[0,good[maxindex]]
	vvmaxmag=dt.vv[1,good[maxindex]]

	oploterror, dt.vv[0,where(finite(dt.vv[1,*]) eq 1)]-vvmaxtime[0],dt.vv[1,where(finite(dt.vv[1,*]) eq 1)]-	vvmaxmag[0],dt.vv[2,where(finite(dt.vv[1,*]) eq 1)], psym=-symbol[n] , 	color=color[n], errcolor=color[n], linestyle=0, symsize=0.5


endfor

al_legend, snlist, color=color[0:n_elements(SNlist)-1], psym=symbol[0:n_elements(SNlist)-1], position=[5,-1], box=1, charsize=0.6, linestyle=1, background='white', symsize=0.5



device, /close
SET_PLOT, 'X'
$open earlysne_v_grey.eps

restore, 'kasenshockmodelcountsnew.sav'


SET_PLOT, 'PS'

device, filename='earlysne_w1_shock.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize, bits_per_pixel=8, /color

x=0
plot, xdata,ydata, /nodata, /noerase, position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle=' Days from Explosion',   ytitle='Absolute UVW1 mag', charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xrange=[0,30],yrange=[-12,-20], ystyle=1, xstyle=1

for n=0,n_elements(snlist)-1 do begin

SNname=SNlist[n]

index=where(host.snname_array eq SNname)
dm=host.dm_best_array[index]

pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat', dt=dt

SNindex=where(host.SNname_array eq SNname)
pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat', dt=dt

good=where(finite(dt.w1[1,*]) eq 1)
maxindex=where(dt.w1[1,good] eq min(dt.w1[1,good]))
w1maxtime=dt.w1[0,good[maxindex]]
w1maxmag=dt.w1[1,good[maxindex]]
maxindex=where(dt.w1[1,good] eq min(dt.w1[1,good]))
w1maxtime=host.peaktime_array[SNindex,2]
w1maxmag=host.appmag_array[SNindex,2]


if SNname eq 'SN2012ij' then w1maxtime=dt.w1[0,good[maxindex+1]]
if SNname eq 'SN2012ij' then w1maxmag=dt.w1[1,good[maxindex+1]]


oploterror, dt.w1[0,where(finite(dt.w1[1,*]) eq 1)]-w1maxtime[0]+16.0,dt.w1[1,where(finite(dt.w1[1,*]) eq 1)]-dm[0],dt.w1[2,where(finite(dt.w1[1,*]) eq 1)], psym=-symbol[n] , color=color[n], errcolor=color[n], linestyle=0, symsize=0.5


endfor

print, a13_array[25]
cgoplot, d_exp,a13mag_array[3,*,25], linestyle=4, color=color[3]
print, a13_array[29]
cgoplot, d_exp,a13mag_array[3,*,29], linestyle=3, color=color[2]
print, a13_array[53]
cgoplot, d_exp,a13mag_array[3,*,53], linestyle=2, color=color[1]
print, a13_array[86]
cgoplot, d_exp,a13mag_array[3,*,86], linestyle=0, color=color[0]

al_legend, ['1 Msun RG', '6 Msun MS', '2 Msun MS', '1 Msun MS'],  color=color[0:3], position=[15,-14.2], box=1, charsize=0.6, linestyle=[0,2,3,4], background='white', symsize=0.5


al_legend, snlist, color=color[0:n_elements(SNlist)-1], psym=symbol[0:n_elements(SNlist)-1], position=[25, -20], box=1, charsize=0.6, linestyle=1, background='white', symsize=0.5



device, /close
SET_PLOT, 'X'
$open earlysne_w1_shock.eps





SET_PLOT, 'PS'

device, filename='earlysne_v_shock.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize, bits_per_pixel=8, /color

x=0
plot, xdata,ydata, /nodata, /noerase, position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle=' Days from Explosion',   ytitle='Absolute V mag', charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
xrange=[0,30],yrange=[-14,-22], ystyle=1, xstyle=1

for n=0,n_elements(snlist)-1 do begin

	SNname=SNlist[n]

	index=where(host.snname_array eq SNname)
	dm=host.dm_best_array[index]


	pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat', dt=dt

	good=where(finite(dt.vv[1,*]) eq 1)
	maxindex=where(dt.vv[1,good] eq min(dt.vv[1,good]))
	vvmaxtime=dt.vv[0,good[maxindex]]
	vvmaxmag=dt.vv[1,good[maxindex]]

	oploterror, dt.vv[0,where(finite(dt.vv[1,*]) eq 1)]-vvmaxtime[0]+19.0,dt.vv[1,where(finite(dt.vv[1,*]) eq 1)]-dm[0],dt.vv[2,where(finite(dt.vv[1,*]) eq 1)], psym=-symbol[n] , color=color[n], errcolor=color[n], linestyle=0, symsize=0.5


endfor

print, a13_array[25]
cgoplot, d_exp,a13mag_array[0,*,25], linestyle=4, color=color[3]
print, a13_array[29]
cgoplot, d_exp,a13mag_array[0,*,29], linestyle=3, color=color[2]
print, a13_array[53]
cgoplot, d_exp,a13mag_array[0,*,53], linestyle=2, color=color[1]
print, a13_array[86]
cgoplot, d_exp,a13mag_array[0,*,86], linestyle=0, color=color[0]

al_legend, ['1 Msun RG', '6 Msun MS', '2 Msun MS', '1 Msun MS'],  color=color[0:3], position=[15,-16], box=1, charsize=0.6, linestyle=[0,2,3,4], background='white', symsize=0.5


al_legend, snlist, color=color[0:n_elements(SNlist)-1], psym=symbol[0:n_elements(SNlist)-1], position=[25,-22], box=1, charsize=0.6, linestyle=1, background='white', symsize=0.5



device, /close
SET_PLOT, 'X'
$open earlysne_v_shock.eps

;;;;;;;;;;;; double plot




; from http://www.iluvatar.org/~dwijn/idlfigures
!p.font = 1
!p.thick = 2
!x.thick = 2
!y.thick = 2
!z.thick = 2
xsize = 8.8
wall = 0.03
margin=0.16
a = xsize/8.8 - (margin + wall)
b = a * 2d / (1 + sqrt(5))

ysize = (margin + b + wall + b + wall)*8.8
ticklen = 0.01
xticklen = ticklen/b
yticklen = ticklen/a

nxticks=20

x1 = margin*8.8/xsize
x2 = x1 + a*8.8/xsize
xc = x2 + wall*8.8/xsize
y1 = margin*8.8/ysize
y2 = y1 + b*8.8/ysize
y3 = y2 + wall*8.8/ysize
y4 = y3 + b*8.8/ysize
yc = y4 + wall*8.8/ysize
xrange=[0,1,2,3]
yrange1=[1,2,3,4]
fontsize=14

SET_PLOT, 'PS'

device, filename='earlysne_w1v_shock.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

plot, xrange, yrange1, /nodata, /noerase, position=[x1,y2,x2,y2+y2-y1], $
xtitle=' ',   ytitle='uvw1 abs mag', charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[-10,-20], ystyle=1, xrange=[-2,38.0], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1)


for n=0,n_elements(snlist)-1 do begin

	SNname=SNlist[n]

	index=where(host.snname_array eq SNname)
	dm=host.dm_best_array[index]

	pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat', dt=dt

SNindex=where(host.SNname_array eq SNname)
pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat', dt=dt

good=where(finite(dt.w1[1,*]) eq 1)
maxindex=where(dt.w1[1,good] eq min(dt.w1[1,good]))
w1maxtime=dt.w1[0,good[maxindex]]
w1maxmag=dt.w1[1,good[maxindex]]
maxindex=where(dt.w1[1,good] eq min(dt.w1[1,good]))
w1maxtime=host.peaktime_array[SNindex,2]
w1maxmag=host.appmag_array[SNindex,2]


	if SNname eq 'SN2012ij' then w1maxtime=dt.w1[0,good[maxindex+1]]
	if SNname eq 'SN2012ij' then w1maxmag=dt.w1[1,good[maxindex+1]]


	oploterror, dt.w1[0,where(finite(dt.w1[1,*]) eq 1)]-w1maxtime[0]+16.0,dt.w1[1,where(finite(dt.w1[1,*]) eq 1)]-dm[0],dt.w1[2,where(finite(dt.w1[1,*]) eq 1)], psym=-symbol[n] , color=color[n], errcolor=color[n], linestyle=0, symsize=0.5


endfor

for n=0,nprimesne-1 do begin

	SNname=SNlist[n]

	index=where(host.snname_array eq SNname)
	dm=host.dm_best_array[index]

	pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat', dt=dt

SNindex=where(host.SNname_array eq SNname)
pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat', dt=dt

good=where(finite(dt.w1[1,*]) eq 1)
maxindex=where(dt.w1[1,good] eq min(dt.w1[1,good]))
w1maxtime=dt.w1[0,good[maxindex]]
w1maxmag=dt.w1[1,good[maxindex]]
maxindex=where(dt.w1[1,good] eq min(dt.w1[1,good]))
w1maxtime=host.peaktime_array[SNindex,2]
w1maxmag=host.appmag_array[SNindex,2]


	if SNname eq 'SN2012ij' then w1maxtime=dt.w1[0,good[maxindex+1]]
	if SNname eq 'SN2012ij' then w1maxmag=dt.w1[1,good[maxindex+1]]


	oploterror, dt.w1[0,where(finite(dt.w1[1,*]) eq 1)]-w1maxtime[0]+16.0,dt.w1[1,where(finite(dt.w1[1,*]) eq 1)]-dm[0],dt.w1[2,where(finite(dt.w1[1,*]) eq 1)], psym=-symbol[n] , color=color[n], errcolor=color[n], linestyle=0, symsize=0.8


endfor


print, a13_array[25]
cgoplot, d_exp,a13mag_array[3,*,25], linestyle=4, color=color[3]
print, a13_array[29]
cgoplot, d_exp,a13mag_array[3,*,29], linestyle=3, color=color[2]
print, a13_array[53]
cgoplot, d_exp,a13mag_array[3,*,53], linestyle=2, color=color[1]
print, a13_array[86]
cgoplot, d_exp,a13mag_array[3,*,86], linestyle=0, color=color[0]



print, 'v'

plot, xrange, yrange1, /nodata, /noerase, position=[x1,y1,x2,y2], $
xtitle='Days from Explosion',   ytitle='v abs mag', charsize=1.0, $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[-13,-19], ystyle=1, xrange=[-2,38.0], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, xtickname=[' ','0',' ','4',' ','8',' ','12',' ','16',' ','20',' ','24',' ','28',' ','32',' ','36',' ']


for n=0,n_elements(snlist)-1 do begin

	SNname=SNlist[n]

	index=where(host.snname_array eq SNname)
	dm=host.dm_best_array[index]

	pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat', dt=dt

	good=where(finite(dt.vv[1,*]) eq 1)
	maxindex=where(dt.vv[1,good] eq min(dt.vv[1,good]))
	vvmaxtime=dt.vv[0,good[maxindex]]
	vvmaxmag=dt.vv[1,good[maxindex]]


	if SNname eq 'iPTF14atg' then vvmaxtime=dt.vv[0,good[maxindex+2]]
	if SNname eq 'iPTF14atg' then vvmaxmag=dt.vv[1,good[maxindex+2]]


	oploterror, dt.vv[0,where(finite(dt.vv[1,*]) eq 1)]-vvmaxtime[0]+19.0,dt.vv[1,where(finite(dt.vv[1,*]) eq 1)]-dm[0],dt.vv[2,where(finite(dt.vv[1,*]) eq 1)], psym=-symbol[n] , color=color[n], errcolor=color[n], linestyle=0, symsize=0.5


endfor

for n=0,nprimesne-1 do begin

	SNname=SNlist[n]

	index=where(host.snname_array eq SNname)
	dm=host.dm_best_array[index]

	pjb_phot_array_B141, '$SOUSA/data/'+SNname+'_uvotB15.1.dat', dt=dt

	good=where(finite(dt.vv[1,*]) eq 1)
	maxindex=where(dt.vv[1,good] eq min(dt.vv[1,good]))
	vvmaxtime=dt.vv[0,good[maxindex]]
	vvmaxmag=dt.vv[1,good[maxindex]]

	if SNname eq 'iPTF14atg' then vvmaxtime=dt.vv[0,good[maxindex+2]]
	if SNname eq 'iPTF14atg' then vvmaxmag=dt.vv[1,good[maxindex+2]]


	oploterror, dt.vv[0,where(finite(dt.vv[1,*]) eq 1)]-vvmaxtime[0]+19.0,dt.vv[1,where(finite(dt.vv[1,*]) eq 1)]-dm[0],dt.vv[2,where(finite(dt.vv[1,*]) eq 1)], psym=-symbol[n] , color=color[n], errcolor=color[n], linestyle=0, symsize=0.8


endfor

print, a13_array[25]
cgoplot, d_exp,a13mag_array[0,*,25], linestyle=4, color=color[3]
print, a13_array[29]
cgoplot, d_exp,a13mag_array[0,*,29], linestyle=3, color=color[2]
print, a13_array[53]
cgoplot, d_exp,a13mag_array[0,*,53], linestyle=2, color=color[1]
print, a13_array[86]
cgoplot, d_exp,a13mag_array[0,*,86], linestyle=0, color=color[0]


;[0:n_elements(SNlist)-1]
al_legend, snlist[0:nprimesne-1], color=color[0:nprimesne-1], psym=symbol[0:nprimesne-1], position=[29,-25], box=1, charsize=0.6, linestyle=1, background='white', symsize=0.5


al_legend, ['1 Msun RG', '6 Msun MS', '2 Msun MS', '1 Msun MS'],  color=color[0:3], position=[16,-15], box=1, charsize=0.6, linestyle=[0,2,3,4], background='white', symsize=0.5


device, /close
SET_PLOT, 'X'

$open *w1v*eps



print, 'final stop'
stop
end