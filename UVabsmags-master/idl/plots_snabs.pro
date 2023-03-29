pro plots_snabs

restore, 'host.sav'

; these are the first order coefficients for the SN1992A spectrum
rlambda_array=[6.20, 8.01,5.43, 4.92,4.16,3.16]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

nplots=6
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
b = a * 2d / (1 + sqrt(5)) *2.0/3.0

ysize = (margin + nplots*(b + wall ) )*xsize
ticklen = 0.01
xticklen = ticklen/b
yticklen = ticklen/a

x1 = margin*8.8/xsize
x2 = x1 + a*8.8/xsize
xc = x2 + wall*8.8/xsize
y1 = margin*8.8/ysize
y2 = y1 + b*8.8/ysize

xdata=[1,2,3,4]
ydata=[2,3,4,5]

xrange=[0.4,2.4]
nxticks=10

yranges=[[-12,-22],[-12,-22],[-12,-22],[-12,-22],[-12,-22],[-12,-22] ]
yranges=[[-10,-20],[-10,-20],[-10,-20],[-10,-20],[-10,-20],[-10,-20] ]
yranges=[[-11,-21],[-11,-21],[-11,-21],[-11,-21],[-11,-21],[-11,-21] ]
nyticks=5

fontsize=12


filtertitles=['uvw2','uvm2','uvw1','u','b','v']
filter=['uvw2','uvm2','uvw1','u','b','v']
colors=['maroon','red','black','violet','blue','dark green']

ytitles=filtertitles+' - $\mu$ -A!Dmw'
;xtitles=filtertitles+' - $\mu$ -A!Dmw'

; xtitles=[tex2idl(\Delta M_{15} (B))]
dmb15='$\Delta$m!D15!N(B)'

xtitles=[' ', ' ', ' ', ' ', ' ', dmb15]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;; set values of constants
H_0=72.0     ;
H_0_err=5.0;
vel_thermal=200.0     ; uncertainty in velocity of galaxies due to thermal/gravitational 
                      ; motion after large scale structure has been removed

SET_PLOT, 'PS'

device, filename='snabs6_loop.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

for f=0,5 do begin

cgplot, xdata, ydata, /nodata, /noerase, position=[x1,(5-f)*b*8.8/ysize+y1,x2,(6-f)*b*8.8/ysize+y1],  $      
ytitle=ytitles[f], charsize=1.0, xtitle=xtitles[f],$
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yranges[*,f], ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1)


oploterror, host.dm15_array[*,4], host.appmag_array[*,f]-host.dm_best_array, host.dm15err_array[*,4], sqrt(host.appmagerr_array[*,f]^2.0+host.dm_best_err_array^2.0),  linestyle=2, psym=6, symsize=symsize, color=cgcolor(colors[f]), errcolor=cgcolor(colors[f]), charsize=1.0, /NoErase, XTickformat='(A1)'


;for n=0, n_elements(host.SNname_array)-1 do print, host.SNname_array[n], host.appmag_array[n,1]-host.dm_best_array[n]

dm15=[0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7]
dm15abs=0.672*(dm15-1.1)+0.633*(dm15-1.1)^2.0-19.0

dm15back=[1.7,1.6,1.5,1.4,1.3,1.1,1,0.9]
dm15backabs=0.672*(dm15back-1.1)+0.633*(dm15back-1.1)^2.0-19.0

xpoly=[dm15,dm15back]
ypoly=[dm15abs-0.09,dm15backabs+0.09]
if f eq 5 then cgcolorfill, xpoly, ypoly, color='grey'

if f eq 0 then axis, xaxis=1, xminor=1, xticklen=xticklen,  xrange=xrange, xstyle=1, xticks=nxticks, xtickv=xtickvalues 

if f eq 5 then axis, xaxis=0, xminor=1, xticklen=xticklen,  xrange=xrange, xstyle=1, xticks=nxticks, xtickv=xtickvalues

endfor


;al_legend, allSNlegend,   psym=allsymnum, symsize=allsizes, color=allcolors, background_color='white', $
;pos=[0.75,0.94], /norm, charsize=0.8, box=1




device, /close
SET_PLOT, 'X'


$open snabs6_loop.eps

stop


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

yranges=[[-14,-20],[-14,-20],[-14,-20],[-14,-20],[-14,-20],[-14,-20] ]
nyticks=nyticks-4

;;;;;;;;;;;;;;;;;;  abs mag plot with forced standard candle ;;;;;;;

SET_PLOT, 'PS'

device, filename='SN_abs_all6_forcebv.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

for f=0,5 do begin

cgplot, xdata, ydata, /nodata, /noerase, position=[x1,(5-f)*b*8.8/ysize+y1,x2,(6-f)*b*8.8/ysize+y1],  $      
ytitle=ytitles[f], charsize=1.0, xtitle=xtitles[f],$
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=yranges[*,f], ystyle=1, xrange=xrange, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1)


oploterror, host.dm15_array[*,4], host.appmag_array[*,f]-host.dm_best_array-rlambda_array[f]*(host.appmag_array[*,4]-host.appmag_array[*,5])  -(host.appmag_array[*,5]-host.dm_best_array-rlambda_array[5]*(host.appmag_array[*,4]-host.appmag_array[*,5])+19.0), host.dm15err_array[*,4], sqrt(host.appmagerr_array[*,f]^2.0+host.dm_best_err_array^2.0),  linestyle=2, psym=6, symsize=symsize, color=cgcolor(colors[f]), errcolor=cgcolor(colors[f]), charsize=1.0, /NoErase, XTickformat='(A1)'


dm15=[0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7]
dm15abs=0.672*(dm15-1.1)+0.633*(dm15-1.1)^2.0-19.0

dm15back=[1.7,1.6,1.5,1.4,1.3,1.1,1,0.9]
dm15backabs=0.672*(dm15back-1.1)+0.633*(dm15back-1.1)^2.0-19.0

xpoly=[dm15,dm15back]
ypoly=[dm15abs-0.09,dm15backabs+0.09]
if f eq 5 then cgcolorfill, xpoly, ypoly, color='grey'

if f eq 0 then axis, xaxis=1, xminor=1, xticklen=xticklen,  xrange=xrange, xstyle=1, xticks=nxticks, xtickv=xtickvalues 

if f eq 5 then axis, xaxis=0, xminor=1, xticklen=xticklen,  xrange=xrange, xstyle=1, xticks=nxticks, xtickv=xtickvalues

endfor


device, /close
SET_PLOT, 'X'
$open SN_abs_all6_forcebv.eps


stop


print, 'final stop'
stop
end
