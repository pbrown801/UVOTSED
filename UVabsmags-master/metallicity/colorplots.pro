pro colorplots, m2b_range=m2b_range, w1b_range=w1b_range, ub_range=ub_range, xrange=xrange, xvalues=xvalues, err_xhigh=err_xhigh, err_xlow=err_xlow, xtitle=xtitle, figurename=figurename, bpeak_mag_array=bpeak_mag_array, bpeak_magerr_array=bpeak_magerr_array

nplots=3
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
margin=0.14
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




SET_PLOT, 'PS'

device, filename=figurename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=12, bits_per_pixel=8, /color

x=0
cgplot, /noerase, xvalues, bpeak_mag_array[*,3]-bpeak_mag_array[*,4],  err_xhigh=err_xhigh, err_xlow=err_xlow, err_ylow=sqrt(bpeak_magerr_array[*,3]-bpeak_magerr_array[*,4]), err_yhigh=sqrt(bpeak_magerr_array[*,3]-bpeak_magerr_array[*,4]), psym=16, symsize=1, color='violet', $
position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle=xtitle,   ytitle='u - b', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
yrange=ub_range, ystyle=1, xrange=xrange, xstyle=1

;y=5-x
x=1
cgplot, /noerase, xvalues, bpeak_mag_array[*,2]-bpeak_mag_array[*,4],  err_xhigh=err_xhigh, err_xlow=err_xlow, err_ylow=sqrt(bpeak_magerr_array[*,2]-bpeak_magerr_array[*,4]), err_yhigh=sqrt(bpeak_magerr_array[*,2]-bpeak_magerr_array[*,4]), psym=16, symsize=1, color='purple', $
position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle=' ',   ytitle='uvw1 - b', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
yrange=w1b_range, ystyle=1, xrange=xrange, xstyle=1
;,  xtickname=replicate(' ',nxticks+1)

x=2

cgplot, /noerase, xvalues, bpeak_mag_array[*,1]-bpeak_mag_array[*,4],  err_xhigh=err_xhigh, err_xlow=err_xlow, err_ylow=sqrt(bpeak_magerr_array[*,1]-bpeak_magerr_array[*,4]), err_yhigh=sqrt(bpeak_magerr_array[*,1]-bpeak_magerr_array[*,4]), psym=16, symsize=1, color='maroon', $
position=[x1,y1+(x)*b*8.8/ysize,x2,y1+(x+1)*b*8.8/ysize], $
xtitle=' ',   ytitle='uvm2 - b', charsize=1.0,  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
yrange=m2b_range, ystyle=1, xrange=xrange, xstyle=1
;, xtickname=replicate(' ',nxticks+1)

device, /close
SET_PLOT, 'X'
; spawn, 'open '+figurename+ ' &'


end


