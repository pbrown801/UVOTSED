pro w1bv

;, dmb15low, dmb15high, bvhigh

dmb15low=0.5
dmb15high=2.5
bvhigh=0.35
bvlow=-0.2

fontsize=16

restore, 'host.sav'  

;  this is the version from Brown et al. 2017 paper with Nancy Landez
; restore, 'nancyhost.sav'  

indexcbv=where(host.snname_array eq 'SN2017cbv')
host.dm15_array[indexcbv,4]=1.06
host.dm15err_array[indexcbv,4]=0.1

index11fe=where(host.snname_array eq 'SN2011fe')
print, host.dm15_array[index11fe,4] 

restore, 'SN2011fe_redbolmags161.sav'
fedm=29.04


indexgi=where(host.snname_array eq 'SN2007gi')
host.dm15_array[indexgi,4]=!Values.F_NAN

; filters, ebv, epochs, laws
;  extinction laws model=0 rv=3.1, 1 1.7  , model 2 smc, model 3 CSLMC
; plot, epoch, feredmags[1,0,*,0]



;normal=where(host.sntype2_array eq 'Ia' and host.dm15_array[*,4] gt dmb15low and host.dm15_array[*,4] lt dmb15high and (host.bpeakappmag_array[*,4]-host.bpeakappmag_array[*,5]) lt bvhigh)
 normal=where(host.sntype2_array eq 'Ia' and (host.bpeakappmag_array[*,4]-host.bpeakappmag_array[*,5]) lt bvhigh and (host.bpeakappmag_array[*,4]-host.bpeakappmag_array[*,5]) gt bvlow and host.DM15ERR_ARRAY[*,4] lt 0.2 and host.dm15_array[*,4] gt dmb15low and host.dm15_array[*,4] lt dmb15high)
 morenormal=where(host.sntype2_array eq 'Ia' and (host.bpeakappmag_array[*,4]-host.bpeakappmag_array[*,5]) lt 0.25 and (host.bpeakappmag_array[*,4]-host.bpeakappmag_array[*,5]) gt -0.15 and host.DM15ERR_ARRAY[*,4] lt 0.2 and host.dm15_array[*,4] gt 1.0 and host.dm15_array[*,4] lt 1.5)
 
ia=where(host.sntype2_array eq 'Ia')

readcol,"../snialist.txt",swiftsnia,quality,format='A,A',/silent
;normal=swiftsnia

firstepoch_array=make_array(n_elements(host.snname_array),value=!Values.F_NAN)
hstspectra_array=make_array(n_elements(host.snname_array),value=' ')
swiftspectra_array=make_array(n_elements(host.snname_array),value=' ')
nuvb_array=make_array(n_elements(host.snname_array),value=' ')
nuvr_array=make_array(n_elements(host.snname_array),value=' ')
milne_array=make_array(n_elements(host.snname_array),value=' ')

stritcolor_array=make_array(n_elements(host.snname_array),value=' ')

readcol,"swiftspectralist.txt",swiftspectra,format='A',/silent
readcol,"hstspectralist.txt",hstspectra,format='A',/silent

;readcol, 'Stritzinger_2018_redblue.txt', SN, Host, Redshift, EBV_MW, EBV_host, t_first, t_rise, DeltamB15, pm, dmerr, M_B, pm2, mberr, SpectralType, Color, References, format='(A, A, F, F, A, F, F, F, F, F, A, F, F, A, F, A, A, A, A)'

readcol, 'Stritzinger_2018_redblue.txt', strit_SNname, strit_Host, strit_Redshift, strit_EBV_MW, strit_EBV_host, pm3, strit_ebvhosterr, strit_t_first, strit_t_rise, strit_DeltamB15, pm, dmerr, strit_M_B, pm2, strit_mberr, strit_Spectralcode, strit_SpectralType, strit_Color, References, format='(A, A, A, A, A,A,A,A,A,A, A, A, A, A, A, A, A, A, A)', comment='#'

milnenuvb=['SN2011ia','SN2011fe','SN2011by','SN2008hv','SN2008Q']
milnenuvr=['SN2015F','SN2013gy','SN2013gs','SN2013ex','SN2013cs','SN2012hr', 'SN2011im','SN2008ec','SN2007co','SN2007af','SN2005df','SN2005cf']
;;; adding in those too red
milnenuvr=['SN2015F','SN2013gy','SN2013gs','SN2013ex','SN2013cs','SN2012hr', 'SN2011im','SN2008ec','SN2007co','SN2007af','SN2005df','SN2005cf', 'SN2013eu','SN2010kg','SN2010gp','SN2010ev']

milnemuvb=['SN2012ht', 'SN2006ej','SN2006dm','SN2010gn']



for n=0, n_elements(host.snname_array) -1 do begin

	SNname= host.snname_array[n]
	print, SNname

	Bpeaktime=host.bpeakmjd_array[n]


	swiftindex=where(swiftspectra eq SNname)
	hstindex=where(hstspectra eq SNname)
	if swiftindex[0] ne -1 then swiftspectra_array[n]='swift'
	if hstindex[0] ne -1  and (host.bpeakappmag_array[n,4]-host.bpeakappmag_array[n,5]) lt 0.35 then hstspectra_array[n]='hst'
	stritindex=where(strit_snname eq SNname)	
	if stritindex[0] ne -1 then stritcolor_array[n]=strit_color[stritindex]

	muvbindex  =where(milnemuvb eq SNname)
	nuvbindex  =where(milnenuvb eq SNname)
	nuvrindex  =where(milnenuvr eq SNname)
	if nuvbindex[0] ne -1 then milne_array[n]='nuvb'
	if muvbindex[0] ne -1 then milne_array[n]='muvb'
	if nuvrindex[0] ne -1 then milne_array[n]='nuvr'
endfor

morenormalhst=where(hstspectra_array eq 'hst' and host.sntype2_array eq 'Ia' and (host.bpeakappmag_array[*,4]-host.bpeakappmag_array[*,5]) lt 0.25 and (host.bpeakappmag_array[*,4]-host.bpeakappmag_array[*,5]) gt -0.15 and host.DM15ERR_ARRAY[*,4] lt 0.2 and host.dm15_array[*,4] gt 1.0 and host.dm15_array[*,4] lt 1.5)
hst=where(hstspectra_array eq 'hst')
swift=where(swiftspectra_array eq 'swift')
stritred=where(stritcolor_array eq 'red')
stritblue=where(stritcolor_array eq 'blue')


muvb=where(milne_array eq 'muvb')
nuvb=where(milne_array eq 'nuvb')
nuvr=where(milne_array eq 'nuvr')





loadct, 33

colortable=intarr(n_elements(normal))
for i=0,n_elements(normal)-1 do colortable[i]=floor( ((host.bpeakappmag_array[normal[i],4]-host.bpeakappmag_array[normal[i],5])+0.1)/0.4*255)
superblue=where(colortable lt 0.0,superbluecount)
superred=where(colortable gt 255,superredcount)
if superbluecount gt 0 then colortable[superblue]=0
if superredcount gt 0 then colortable[superred]=255
 

morecolortable=intarr(n_elements(morenormal))
for i=0,n_elements(morenormal)-1 do morecolortable[i]=floor( ((host.bpeakappmag_array[morenormal[i],4]-host.bpeakappmag_array[morenormal[i],5])+0.1)/0.4*255)
superblue=where(morecolortable lt 0.0,superbluecount)
superred=where(morecolortable gt 255,superredcount)
if superbluecount gt 0 then morecolortable[superblue]=0
if superredcount gt 0 then morecolortable[superred]=255
 

; plot, host.dm15_array[normal[n],4],w1vbluest_array[2,n]-w1vbluest_array[5,n]

;for n=0, n_elements(host.snname_array[normal]) -1 do print, host.snname_array[normal[n]], w1vbluest_array[2,n]-w1vbluest_array[5,n]
; save, 

cgWindow_SetDefs, PS_Decomposed=1



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


plotfilename = 'bpeak_w1vbv_spectra.eps'

SET_PLOT, 'PS'

device, filename= plotfilename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

restore, 'SN2011fe_redbolmags161.sav'
; filters, ebv, epochs, laws
;  extinction laws model=0 rv=3.1, 1 1.7  , model 2 smc, model 3 CSLMC
; plot, epoch, feredmags[1,0,*,0]
febpeak=where(min(feredmags[4,0,*,0]) eq feredmags[4,0,*,0],count)

; print, ((feredmags[2,10,febpeak,0]-feredmags[5,10,febpeak,0])-(feredmags[2,0,febpeak,0]-feredmags[5,0,febpeak,0]))/((feredmags[4,10,febpeak,0]-feredmags[5,10,febpeak,0])-(feredmags[4,0,febpeak,0]-feredmags[5,0,febpeak,0]))

;print, ((feredmags[2,10,febpeak,1]-feredmags[5,10,febpeak,1])-(feredmags[2,0,febpeak,1]-feredmags[5,0,febpeak,1]))/((feredmags[4,10,febpeak,1]-feredmags[5,10,febpeak,1])-(feredmags[4,0,febpeak,1]-feredmags[5,0,febpeak,1]))
;print, ((feredmags[2,10,febpeak,2]-feredmags[5,10,febpeak,2])-(feredmags[2,0,febpeak,2]-feredmags[5,0,febpeak,2]))/((feredmags[4,10,febpeak,2]-feredmags[5,10,febpeak,2])-(feredmags[4,0,febpeak,2]-feredmags[5,0,febpeak,2]))
;print, ((feredmags[2,10,febpeak,3]-feredmags[5,10,febpeak,3])-(feredmags[2,0,febpeak,3]-feredmags[5,0,febpeak,3]))/((feredmags[4,10,febpeak,3]-feredmags[5,10,febpeak,3])-(feredmags[4,0,febpeak,3]-feredmags[5,0,febpeak,3]))

cgplot, charsize=1, feredmags[4,*,febpeak,3]-feredmags[5,*,febpeak,3], feredmags[2,*,febpeak,3]-feredmags[5,*,febpeak,3], xrange=[-0.15,0.3], yrange=[0.8,2.2], ystyle=1, xstyle=1, ytitle='(w1-v)!BB!Lpeak', $ 
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
oploterror, host.bpeakappmag_array[*,4]-host.bpeakappmag_array[*,5], host.bpeakappmag_array[*,2]-host.bpeakappmag_array[*,5], sqrt(host.bpeakappmagerr_array[*,4]^2.0+host.bpeakappmagerr_array[*,5]^2.0), sqrt(host.bpeakappmagerr_array[*,2]^2.0+host.bpeakappmagerr_array[*,5]^2.0), symsize=0.3, psym=15

if hst[0] ne -1 then cgoplot, host.bpeakappmag_array[4,hst]-host.bpeakappmag_array[5,hst], host.bpeakappmag_array[2,hst]-host.bpeakappmag_array[5,hst], psym=16, symsize=1, color='red'

if swift[0] ne -1 then cgoplot, host.bpeakappmag_array[4,swift]-host.bpeakappmag_array[5,swift], host.bpeakappmag_array[2,swift]-host.bpeakappmag_array[5,swift], psym=46, symsize=1.2, color='blue'

xyouts, host.bpeakappmag_array[4,hst]-host.bpeakappmag_array[5,hst] - 0.06, host.bpeakappmag_array[2,hst]-host.bpeakappmag_array[5,hst]+0.02, host.snname_array[hst], charsize=0.5
xyouts, host.bpeakappmag_array[4,swift]-host.bpeakappmag_array[5,swift] - 0.06, host.bpeakappmag_array[2,swift]-host.bpeakappmag_array[5,swift]+0.02, host.snname_array[normal[swift]], charsize=0.5

al_legend, ['HST','Swift'], psym=[16,46], color=['red','blue'], symsize=[1,1.2], $
pos=[0.7,0.45], /norm, charsize=0.8, box=1

device, /close
SET_PLOT, 'X'
;spawn, 'open bpeak_w1vbv_spectra.eps'

;;;;;;;;;;;;;;;;;;;;;;;;

plotfilename = 'bpeak_m2w1vw1v_spectra.eps'

SET_PLOT, 'PS'

device, filename= plotfilename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

restore, 'SN2011fe_redbolmags161.sav'
; filters, ebv, epochs, laws
;  extinction laws model=0 rv=3.1, 1 1.7  , model 2 smc, model 3 CSLMC
; plot, epoch, feredmags[1,0,*,0]

cgplot, charsize=1, feredmags[1,*,febpeak,3]-feredmags[2,*,febpeak,3], feredmags[2,*,febpeak,3]-feredmags[5,*,febpeak,3], xrange=[1,3.5], yrange=[0.8,2.2], $ 
xtitle=' !A (m2-w1)!NBpeak   ', $
ystyle=1, xstyle=1, ytitle='(w1-v)!BB!Lpeak', position=[x1,y1,x2,y2], linestyle=0, color='black'
oplot, feredmags[1,*,febpeak,2]-feredmags[2,*,febpeak,2], feredmags[2,*,febpeak,2]-feredmags[5,*,febpeak,2], linestyle=1
oplot, feredmags[1,*,febpeak,1]-feredmags[2,*,febpeak,1], feredmags[2,*,febpeak,1]-feredmags[5,*,febpeak,1], linestyle=2
oplot, feredmags[1,*,febpeak,0]-feredmags[2,*,febpeak,0], feredmags[2,*,febpeak,0]-feredmags[5,*,febpeak,0], linestyle=3

;;;;;;;;;;;;;;;;;
oploterror, host.bpeakappmag_array[*,1]-host.bpeakappmag_array[*,2], host.bpeakappmag_array[*,2]-host.bpeakappmag_array[*,5], sqrt(host.bpeakappmagerr_array[*,1]^2.0+host.bpeakappmagerr_array[*,2]^2.0), sqrt(host.bpeakappmagerr_array[*,2]^2.0+host.bpeakappmagerr_array[*,5]^2.0), symsize=0.3, psym=15

if hst[0] ne -1 then cgoplot, host.bpeakappmag_array[1,hst]-host.bpeakappmag_array[2,hst], host.bpeakappmag_array[2,hst]-host.bpeakappmag_array[5,hst], psym=16, symsize=1, color='red'

if swift[0] ne -1 then cgoplot, host.bpeakappmag_array[1,swift]-host.bpeakappmag_array[2,swift], host.bpeakappmag_array[2,swift]-host.bpeakappmag_array[5,swift], psym=46, symsize=1.2, color='blue'

xyouts, host.bpeakappmag_array[1,hst]-host.bpeakappmag_array[2,hst] - 0.06, host.bpeakappmag_array[2,hst]-host.bpeakappmag_array[5,hst]+0.02, host.snname_array[hst], charsize=0.5
xyouts, host.bpeakappmag_array[1,swift]-host.bpeakappmag_array[2,swift] - 0.06, host.bpeakappmag_array[2,swift]-host.bpeakappmag_array[5,swift]+0.02, host.snname_array[normal[swift]], charsize=0.5

al_legend, ['HST','Swift'], psym=[16,46], color=['red','blue'], symsize=[1,1.2], $
pos=[0.7,0.45], /norm, charsize=0.8, box=1

device, /close
SET_PLOT, 'X'
;spawn, 'open bpeak_m2w1vw1v_spectra.eps'

;;;;;;;;;;;;;;;;;;;;;;;;

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
oploterror, host.bpeakappmag_array[*,4]-host.bpeakappmag_array[*,5], host.bpeakappmag_array[*,2]-host.bpeakappmag_array[*,5], sqrt(host.bpeakappmagerr_array[*,4]^2.0+host.bpeakappmagerr_array[*,5]^2.0), sqrt(host.bpeakappmagerr_array[*,2]^2.0+host.bpeakappmagerr_array[*,5]^2.0), symsize=0.3, psym=15

cgoplot, host.bpeakappmag_array[4,stritred]-host.bpeakappmag_array[5,stritred], host.bpeakappmag_array[2,stritred]-host.bpeakappmag_array[5,stritred], psym=16, symsize=1, color='red'

cgoplot, host.bpeakappmag_array[4,stritblue]-host.bpeakappmag_array[5,stritblue], host.bpeakappmag_array[2,stritblue]-host.bpeakappmag_array[5,stritblue], psym=46, symsize=1.2, color='blue'

;xyouts, host.bpeakappmag_array[4,hst]-host.bpeakappmag_array[5,hst] - 0.06, host.bpeakappmag_array[2,hst]-host.bpeakappmag_array[5,hst]+0.02, host.snname_array[hst], charsize=0.5
;xyouts, host.bpeakappmag_array[4,swift]-host.bpeakappmag_array[5,swift] - 0.06, host.bpeakappmag_array[2,swift]-host.bpeakappmag_array[5,swift]+0.02, host.snname_array[normal[swift]], charsize=0.5

al_legend, ['Strit+18 red','Strit+18 blue'], psym=[16,46], color=['red','blue'], symsize=[1,1.2], $
pos=[0.5,0.45], /norm, charsize=0.8, box=1

device, /close
SET_PLOT, 'X'
;spawn, 'open bpeak_w1vbv_strit.eps'

;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

SET_PLOT, 'PS'

device, filename='w1v_bv_bpeak_nuvbr.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

xrange=[-0.2,0.2]
nxticks=4

cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y1,x2,y2], $
xtitle='b-v',   ytitle='(w1-v)!BB!Lpeak', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[0.5,2.5], ystyle=1, xrange=[-0.2,0.3], xstyle=1, $
xticks=5, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, color='black'


;;;;;;;;;;   at b band maximum light


cgoplot, color='black', feredmags[4,*,febpeak,3]-feredmags[5,*,febpeak,3], feredmags[2,*,febpeak,3]-feredmags[5,*,febpeak,3], linestyle=0
cgoplot, color='black', feredmags[4,*,febpeak,2]-feredmags[5,*,febpeak,2], feredmags[2,*,febpeak,2]-feredmags[5,*,febpeak,2], linestyle=1
cgoplot, color='black', feredmags[4,*,febpeak,1]-feredmags[5,*,febpeak,1], feredmags[2,*,febpeak,1]-feredmags[5,*,febpeak,1], linestyle=2
cgoplot, color='black', feredmags[4,*,febpeak,0]-feredmags[5,*,febpeak,0], feredmags[2,*,febpeak,0]-feredmags[5,*,febpeak,0], linestyle=3




w1vsne=where(finite(host.bpeakappmag_array[*,2]) eq 1 and finite(host.bpeakappmag_array[*,4]) eq 1 and finite(host.bpeakappmag_array[*,5]) eq 1 )

nuvrw1vsne=where(finite(host.bpeakappmag_array[*,2]) eq 1 and finite(host.bpeakappmag_array[*,4]) eq 1 and finite(host.bpeakappmag_array[*,5]) eq 1  and milne_array eq 'nuvr')

nuvbw1vsne=where(finite(host.bpeakappmag_array[*,2]) eq 1 and finite(host.bpeakappmag_array[*,4]) eq 1 and finite(host.bpeakappmag_array[*,5]) eq 1  and milne_array eq 'nuvb')

muvbw1vsne=where(finite(host.bpeakappmag_array[*,2]) eq 1 and finite(host.bpeakappmag_array[*,4]) eq 1 and finite(host.bpeakappmag_array[*,5]) eq 1  and milne_array eq 'muvb')

oploterror, host.bpeakappmag_array[4,w1vsne]-host.bpeakappmag_array[5,w1vsne], host.bpeakappmag_array[2,w1vsne]-host.bpeakappmag_array[5,w1vsne], sqrt(host.bpeakappmagerr_array[4,w1vsne]^2.0+host.bpeakappmagerr_array[5,w1vsne]^2.0), sqrt(host.bpeakappmagerr_array[2,w1vsne]^2.0+host.bpeakappmagerr_array[5,w1vsne]^2.0), psym=15, symsize=0.5, color='black'

oploterror, host.bpeakappmag_array[4,nuvrw1vsne]-host.bpeakappmag_array[5,nuvrw1vsne], host.bpeakappmag_array[2,nuvrw1vsne]-host.bpeakappmag_array[5,nuvrw1vsne], sqrt(host.bpeakappmagerr_array[4,nuvrw1vsne]^2.0+host.bpeakappmagerr_array[5,nuvrw1vsne]^2.0), sqrt(host.bpeakappmagerr_array[2,nuvrw1vsne]^2.0+host.bpeakappmagerr_array[5,nuvrw1vsne]^2.0), psym=16, color='red'

oploterror, host.bpeakappmag_array[4,nuvbw1vsne]-host.bpeakappmag_array[5,nuvbw1vsne], host.bpeakappmag_array[2,nuvbw1vsne]-host.bpeakappmag_array[5,nuvbw1vsne], sqrt(host.bpeakappmagerr_array[4,nuvbw1vsne]^2.0+host.bpeakappmagerr_array[5,nuvbw1vsne]^2.0), sqrt(host.bpeakappmagerr_array[2,nuvbw1vsne]^2.0+host.bpeakappmagerr_array[5,nuvbw1vsne]^2.0), psym=15, color='royal blue'


oploterror, host.bpeakappmag_array[4,muvbw1vsne]-host.bpeakappmag_array[5,muvbw1vsne], host.bpeakappmag_array[2,muvbw1vsne]-host.bpeakappmag_array[5,muvbw1vsne], sqrt(host.bpeakappmagerr_array[4,muvbw1vsne]^2.0+host.bpeakappmagerr_array[5,muvbw1vsne]^2.0), sqrt(host.bpeakappmagerr_array[2,muvbw1vsne]^2.0+host.bpeakappmagerr_array[5,muvbw1vsne]^2.0), psym=46, color='dark green', symsize=1.2

al_legend, ['NUV-red','MUV-blue','NUV-blue'], psym=[16,46,15], color=['red', 'dark green','royal blue'], position=[-0.2,2.5], box=0, charsize=0.9

;xyouts, host.bpeakappmag_array[4,w1vsne]-host.bpeakappmag_array[5,w1vsne] + 0.01, host.bpeakappmag_array[2,w1vsne]-host.bpeakappmag_array[5,w1vsne] - 0.02, host.snname_array[normal[w1vsne]], charsize=0.5

device, /close
SET_PLOT, 'X'


;;spawn, 'open w1v_bv_bpeak_nuvbr.eps'

;;;;;;;;;

SET_PLOT, 'PS'

device, filename='w1m2_bv_bpeak_nuvbr.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

xrange=[-0.2,0.2]
nxticks=4

cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y1,x2,y2], $
xtitle='b-v',   ytitle='(m2-w1)!BB!Lpeak', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[0.5,4.5], ystyle=1, xrange=[-0.2,0.3], xstyle=1, $
xticks=5, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, color='black'


;;;;;;;;;;   at b band maximum light


cgoplot, color='black', feredmags[4,*,febpeak,3]-feredmags[5,*,febpeak,3], feredmags[1,*,febpeak,3]-feredmags[2,*,febpeak,3], linestyle=0
cgoplot, color='black', feredmags[4,*,febpeak,2]-feredmags[5,*,febpeak,2], feredmags[1,*,febpeak,2]-feredmags[2,*,febpeak,2], linestyle=1
cgoplot, color='black', feredmags[4,*,febpeak,1]-feredmags[5,*,febpeak,1], feredmags[1,*,febpeak,1]-feredmags[2,*,febpeak,1], linestyle=2
cgoplot, color='black', feredmags[4,*,febpeak,0]-feredmags[5,*,febpeak,0], feredmags[1,*,febpeak,0]-feredmags[2,*,febpeak,0], linestyle=3




w1bvsne=where(finite(host.bpeakappmag_array[*,1]) eq 1 and finite(host.bpeakappmag_array[*,2]) eq 1 and finite(host.bpeakappmag_array[*,4]) eq 1 and finite(host.bpeakappmag_array[*,5]) eq 1 )

nuvrw1vsne=where(finite(host.bpeakappmag_array[*,1]) eq 1 and finite(host.bpeakappmag_array[*,2]) eq 1 and finite(host.bpeakappmag_array[*,4]) eq 1 and finite(host.bpeakappmag_array[*,5]) eq 1  and milne_array eq 'nuvr')

nuvbw1vsne=where(finite(host.bpeakappmag_array[*,1]) eq 1 and finite(host.bpeakappmag_array[*,2]) eq 1 and finite(host.bpeakappmag_array[*,4]) eq 1 and finite(host.bpeakappmag_array[*,5]) eq 1  and milne_array eq 'nuvb')

muvbw1vsne=where(finite(host.bpeakappmag_array[*,1]) eq 1 and finite(host.bpeakappmag_array[*,2]) eq 1 and finite(host.bpeakappmag_array[*,4]) eq 1 and finite(host.bpeakappmag_array[*,5]) eq 1  and milne_array eq 'muvb')

oploterror, host.bpeakappmag_array[4,w1vsne]-host.bpeakappmag_array[5,w1vsne], host.bpeakappmag_array[1,w1vsne]-host.bpeakappmag_array[2,w1vsne], sqrt(host.bpeakappmagerr_array[4,w1vsne]^2.0+host.bpeakappmagerr_array[5,w1vsne]^2.0), sqrt(host.bpeakappmagerr_array[2,w1vsne]^2.0+host.bpeakappmagerr_array[1,w1vsne]^2.0), psym=15, symsize=0.5, color='black'

oploterror, host.bpeakappmag_array[4,nuvrw1vsne]-host.bpeakappmag_array[5,nuvrw1vsne], host.bpeakappmag_array[1,nuvrw1vsne]-host.bpeakappmag_array[2,nuvrw1vsne], sqrt(host.bpeakappmagerr_array[4,nuvrw1vsne]^2.0+host.bpeakappmagerr_array[5,nuvrw1vsne]^2.0), sqrt(host.bpeakappmagerr_array[2,nuvrw1vsne]^2.0+host.bpeakappmagerr_array[1,nuvrw1vsne]^2.0), psym=16, color='red'

oploterror, host.bpeakappmag_array[4,nuvbw1vsne]-host.bpeakappmag_array[5,nuvbw1vsne], host.bpeakappmag_array[1,nuvbw1vsne]-host.bpeakappmag_array[2,nuvbw1vsne], sqrt(host.bpeakappmagerr_array[4,nuvbw1vsne]^2.0+host.bpeakappmagerr_array[5,nuvbw1vsne]^2.0), sqrt(host.bpeakappmagerr_array[2,nuvbw1vsne]^2.0+host.bpeakappmagerr_array[5,nuvbw1vsne]^2.0), psym=15, color='royal blue'


oploterror, host.bpeakappmag_array[4,muvbw1vsne]-host.bpeakappmag_array[5,muvbw1vsne], host.bpeakappmag_array[1,muvbw1vsne]-host.bpeakappmag_array[2,muvbw1vsne], sqrt(host.bpeakappmagerr_array[4,muvbw1vsne]^2.0+host.bpeakappmagerr_array[5,muvbw1vsne]^2.0), sqrt(host.bpeakappmagerr_array[2,muvbw1vsne]^2.0+host.bpeakappmagerr_array[1,muvbw1vsne]^2.0), psym=46, color='dark green', symsize=1.2

al_legend, ['NUV-red','MUV-blue','NUV-blue'], psym=[16,46,15], color=['red', 'dark green','royal blue'], position=[-0.2,4.5], box=0, charsize=0.9

;xyouts, host.bpeakappmag_array[4,w1vsne]-host.bpeakappmag_array[5,w1vsne] + 0.01, host.bpeakappmag_array[2,w1vsne]-host.bpeakappmag_array[5,w1vsne] - 0.02, host.snname_array[normal[w1vsne]], charsize=0.5

device, /close
SET_PLOT, 'X'


;;spawn, 'open w1m2_bv_bpeak_nuvbr.eps'



;;;;;;;;;

plotfilename = 'bpeak_m2w1vw1v_hst.eps'

SET_PLOT, 'PS'

device, filename= plotfilename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

restore, 'SN2011fe_redbolmags161.sav'
; filters, ebv, epochs, laws
;  extinction laws model=0 rv=3.1, 1 1.7  , model 2 smc, model 3 CSLMC
; plot, epoch, feredmags[1,0,*,0]

cgplot, charsize=1, feredmags[1,*,febpeak,3]-feredmags[2,*,febpeak,3], feredmags[2,*,febpeak,3]-feredmags[5,*,febpeak,3], xrange=[0.5,4], yrange=[0,3], $ 
xtitle=' !A (m2-w1)!NBpeak   ', $
ystyle=1, xstyle=1, ytitle='(w1-v)!BB!Lpeak', position=[x1,y1,x2,y2], linestyle=0, color='black'
oplot, feredmags[1,*,febpeak,2]-feredmags[2,*,febpeak,2], feredmags[2,*,febpeak,2]-feredmags[5,*,febpeak,2], linestyle=1
oplot, feredmags[1,*,febpeak,1]-feredmags[2,*,febpeak,1], feredmags[2,*,febpeak,1]-feredmags[5,*,febpeak,1], linestyle=2
oplot, feredmags[1,*,febpeak,0]-feredmags[2,*,febpeak,0], feredmags[2,*,febpeak,0]-feredmags[5,*,febpeak,0], linestyle=3

;;;;;;;;;;;;;;;;;
oploterror, host.bpeakappmag_array[*,1]-host.bpeakappmag_array[*,2], host.bpeakappmag_array[*,2]-host.bpeakappmag_array[*,5], sqrt(host.bpeakappmagerr_array[*,1]^2.0+host.bpeakappmagerr_array[*,2]^2.0), sqrt(host.bpeakappmagerr_array[*,2]^2.0+host.bpeakappmagerr_array[*,5]^2.0), symsize=0.3, psym=15

for i=0,n_elements(normal)-1 do  cgplots, host.bpeakappmag_array[normal[i],1]-host.bpeakappmag_array[normal[i],2], host.bpeakappmag_array[normal[i],2]-host.bpeakappmag_array[normal[i],5], psym=16, color=colortable[i], symsize=1.0



oploterror, host.bpeakappmag_array[*,1]-host.bpeakappmag_array[*,2], host.bpeakappmag_array[*,2]-host.bpeakappmag_array[*,5], sqrt(host.bpeakappmagerr_array[*,1]^2.0+host.bpeakappmagerr_array[*,2]^2.0),  sqrt(host.bpeakappmagerr_array[*,2]^2.0+host.bpeakappmagerr_array[*,5]^2.0), psym=16, color='black', symsize=0.1


if hst[0] ne -1 then cgoplot, host.bpeakappmag_array[1,hst]-host.bpeakappmag_array[2,hst], host.bpeakappmag_array[2,hst]-host.bpeakappmag_array[5,hst], psym=9, symsize=1.5, color='black'

device, /close
SET_PLOT, 'X'
;spawn, 'open bpeak_m2w1vw1v_hst.eps'


;;;;;;;;;

plotfilename = 'bpeak_m2vvw1v.eps'

SET_PLOT, 'PS'

device, filename= plotfilename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

restore, 'SN2011fe_redbolmags161.sav'
; filters, ebv, epochs, laws
;  extinction laws model=0 rv=3.1, 1 1.7  , model 2 smc, model 3 CSLMC
; plot, epoch, feredmags[1,0,*,0]
xrange=[0.5,5.5]
yrange=[0,3]
;cgplot, charsize=1, feredmags[1,*,febpeak,3]-feredmags[5,*,febpeak,3], feredmags[2,*,febpeak,3]-feredmags[5,*,febpeak,3], $
cgplot, xdata, ydata, /nodata, /noerase, $
xrange=xrange, yrange=yrange, $ 
xtitle=' !A (m2-v)!NBpeak   ', charsize=1,  $
ystyle=1, xstyle=1, ytitle='(w1-v)!BB!Lpeak', position=[x1,y1,x2,y2], linestyle=0, color='black'
;oplot, feredmags[1,*,febpeak,2]-feredmags[5,*,febpeak,2], feredmags[2,*,febpeak,2]-feredmags[5,*,febpeak,2], linestyle=1
;oplot, feredmags[1,*,febpeak,1]-feredmags[5,*,febpeak,1], feredmags[2,*,febpeak,1]-feredmags[5,*,febpeak,1], linestyle=2
;oplot, feredmags[1,*,febpeak,0]-feredmags[5,*,febpeak,0], feredmags[2,*,febpeak,0]-feredmags[5,*,febpeak,0], linestyle=3

;;;;;;;;;;;;;;;;;
oploterror, host.bpeakappmag_array[*,1]-host.bpeakappmag_array[*,5], host.bpeakappmag_array[*,2]-host.bpeakappmag_array[*,5], sqrt(host.bpeakappmagerr_array[*,1]^2.0+host.bpeakappmagerr_array[*,5]^2.0), sqrt(host.bpeakappmagerr_array[*,2]^2.0+host.bpeakappmagerr_array[*,5]^2.0), symsize=0.3, psym=15

for i=0,n_elements(colortable)-1 do  cgplots, host.bpeakappmag_array[i,1]-host.bpeakappmag_array[i,5], host.bpeakappmag_array[i,2]-host.bpeakappmag_array[i,5], psym=16, color=colortable[i], symsize=1.0



oploterror, host.bpeakappmag_array[*,1]-host.bpeakappmag_array[*,5], host.bpeakappmag_array[*,2]-host.bpeakappmag_array[*,5], sqrt(host.bpeakappmagerr_array[*,1]^2.0+host.bpeakappmagerr_array[*,5]^2.0),  sqrt(host.bpeakappmagerr_array[*,2]^2.0+host.bpeakappmagerr_array[*,5]^2.0), psym=16, color='black', symsize=0.1

oplot, xrange, [1.0,1.0]
oplot, [2.7,2.7], yrange

xyouts, host.bpeakappmag_array[where(host.snname_array eq 'iPTF14atg'),1]-host.bpeakappmag_array[where(host.snname_array eq 'iPTF14atg'),5]-0.3, host.bpeakappmag_array[where(host.snname_array eq 'iPTF14atg'),2]-host.bpeakappmag_array[where(host.snname_array eq 'iPTF14atg'),5] + 0.15, 'iPTF14atg', charsize=0.5

xyouts, host.bpeakappmag_array[where(host.snname_array eq 'SN2016ccj'),1]-host.bpeakappmag_array[where(host.snname_array eq 'SN2016ccj'),5]-0.3, host.bpeakappmag_array[where(host.snname_array eq 'SN2016ccj'),2]-host.bpeakappmag_array[where(host.snname_array eq 'SN2016ccj'),5] - 0.3, 'SN2016ccj', charsize=0.5


xyouts, host.bpeakappmag_array[where(host.snname_array eq 'ASASSN-15pz'),1]-host.bpeakappmag_array[where(host.snname_array eq 'ASASSN-15pz'),5]-0.0, host.bpeakappmag_array[where(host.snname_array eq 'ASASSN-15pz'),2]-host.bpeakappmag_array[where(host.snname_array eq 'ASASSN-15pz'),5] - 0.15, 'A15pz', charsize=0.5

xyouts, host.bpeakappmag_array[where(host.snname_array eq 'SN2011aa'),1]-host.bpeakappmag_array[where(host.snname_array eq 'SN2011aa'),5]-0.0, host.bpeakappmag_array[where(host.snname_array eq 'SN2011aa'),2]-host.bpeakappmag_array[where(host.snname_array eq 'SN2011aa'),5] + 0.05, 'SN2011aa', charsize=0.5

xyouts, host.bpeakappmag_array[where(host.snname_array eq 'SN2012dn'),1]-host.bpeakappmag_array[where(host.snname_array eq 'SN2012dn'),5]-0.5, host.bpeakappmag_array[where(host.snname_array eq 'SN2012dn'),2]-host.bpeakappmag_array[where(host.snname_array eq 'SN2012dn'),5] - 0.2, 'SN2012dn', charsize=0.5

xyouts, host.bpeakappmag_array[where(host.snname_array eq 'LSQ12gdj'),1]-host.bpeakappmag_array[where(host.snname_array eq 'LSQ12gdj'),5]+0.15, host.bpeakappmag_array[where(host.snname_array eq 'LSQ12gdj'),2]-host.bpeakappmag_array[where(host.snname_array eq 'LSQ12gdj'),5] - 0.3, 'LSQ12gdj', charsize=0.5
xyouts, host.bpeakappmag_array[where(host.snname_array eq 'LSQ12gdj'),1]-host.bpeakappmag_array[where(host.snname_array eq 'LSQ12gdj'),5]+0.05, host.bpeakappmag_array[where(host.snname_array eq 'LSQ12gdj'),2]-host.bpeakappmag_array[where(host.snname_array eq 'LSQ12gdj'),5] - 0.2, '\', charsize=0.5



xyouts, host.bpeakappmag_array[where(host.snname_array eq 'ASASSN-15rq'),1]-host.bpeakappmag_array[where(host.snname_array eq 'ASASSN-15rq'),5]-0.2, host.bpeakappmag_array[where(host.snname_array eq 'ASASSN-15rq'),2]-host.bpeakappmag_array[where(host.snname_array eq 'ASASSN-15rq'),5] - 0.25, 'A15rq', charsize=0.5

xyouts, host.bpeakappmag_array[where(host.snname_array eq 'SN2016bln'),1]-host.bpeakappmag_array[where(host.snname_array eq 'SN2016bln'),5]+0.05, host.bpeakappmag_array[where(host.snname_array eq 'SN2016bln'),2]-host.bpeakappmag_array[where(host.snname_array eq 'SN2016bln'),5] - 0.15, 'SN2016bln', charsize=0.5

xyouts, host.bpeakappmag_array[where(host.snname_array eq 'SN2011ay'),1]-host.bpeakappmag_array[where(host.snname_array eq 'SN2011ay'),5]-0.9, host.bpeakappmag_array[where(host.snname_array eq 'SN2011ay'),2]-host.bpeakappmag_array[where(host.snname_array eq 'SN2011ay'),5] + 0.05, 'SN2011ay', charsize=0.5

xyouts, host.bpeakappmag_array[where(host.snname_array eq 'SN2011ia'),1]-host.bpeakappmag_array[where(host.snname_array eq 'SN2011ia'),5]+0.05, host.bpeakappmag_array[where(host.snname_array eq 'SN2011ia'),2]-host.bpeakappmag_array[where(host.snname_array eq 'SN2011ia'),5] - 0.15, 'SN2011ia', charsize=0.5


device, /close
SET_PLOT, 'X'
spawn, 'open bpeak_m2vvw1v.eps'

;;;;;;
;;;;;;;;;

plotfilename = 'bpeak_uvvw1v.eps'

SET_PLOT, 'PS'

device, filename= plotfilename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

restore, 'SN2011fe_redbolmags161.sav'
; filters, ebv, epochs, laws
;  extinction laws model=0 rv=3.1, 1 1.7  , model 2 smc, model 3 CSLMC
; plot, epoch, feredmags[1,0,*,0]

cgplot, charsize=1, feredmags[3,*,febpeak,3]-feredmags[5,*,febpeak,3], feredmags[2,*,febpeak,3]-feredmags[5,*,febpeak,3], xrange=[-1.5,3], yrange=[-1,5], $ 
xtitle=' !A (u-v)!NBpeak   ', $
ystyle=1, xstyle=1, ytitle='(w1-v)!BB!Lpeak', position=[x1,y1,x2,y2], linestyle=0, color='black'
oplot, feredmags[3,*,febpeak,2]-feredmags[5,*,febpeak,2], feredmags[2,*,febpeak,2]-feredmags[5,*,febpeak,2], linestyle=1
oplot, feredmags[3,*,febpeak,1]-feredmags[5,*,febpeak,1], feredmags[2,*,febpeak,1]-feredmags[5,*,febpeak,1], linestyle=2
oplot, feredmags[3,*,febpeak,0]-feredmags[5,*,febpeak,0], feredmags[2,*,febpeak,0]-feredmags[5,*,febpeak,0], linestyle=3

;;;;;;;;;;;;;;;;;
oploterror, host.bpeakappmag_array[*,3]-host.bpeakappmag_array[*,5], host.bpeakappmag_array[*,2]-host.bpeakappmag_array[*,5], sqrt(host.bpeakappmagerr_array[*,3]^2.0+host.bpeakappmagerr_array[*,5]^2.0), sqrt(host.bpeakappmagerr_array[*,2]^2.0+host.bpeakappmagerr_array[*,5]^2.0), symsize=0.3, psym=15

for i=0,n_elements(normal)-1 do  cgplots, host.bpeakappmag_array[i,3]-host.bpeakappmag_array[i,5], host.bpeakappmag_array[i,2]-host.bpeakappmag_array[i,5], psym=16, color=colortable[i], symsize=1.0



oploterror, host.bpeakappmag_array[*,3]-host.bpeakappmag_array[*,5], host.bpeakappmag_array[*,2]-host.bpeakappmag_array[*,5], sqrt(host.bpeakappmagerr_array[*,3]^2.0+host.bpeakappmagerr_array[*,5]^2.0),  sqrt(host.bpeakappmagerr_array[*,2]^2.0+host.bpeakappmagerr_array[*,5]^2.0), psym=16, color='black', symsize=0.1


device, /close
SET_PLOT, 'X'
;spawn, 'open bpeak_uvvw1v.eps'

;;;;;;

;;;;;;;;;;;;;;;;;;;;
plotfilename = 'bpeak_bvw1v_hst.eps'

SET_PLOT, 'PS'

device, filename= plotfilename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

restore, 'SN2011fe_redbolmags161.sav'
; filters, ebv, epochs, laws
;  extinction laws model=0 rv=3.1, 1 1.7  , model 2 smc, model 3 CSLMC
; plot, epoch, feredmags[1,0,*,0]

cgplot, charsize=1, feredmags[4,*,febpeak,3]-feredmags[5,*,febpeak,3], feredmags[2,*,febpeak,3]-feredmags[5,*,febpeak,3], xrange=[-0.15,0.25], yrange=[0.5,2.5], $ 
xtitle=' !A (B-V)!NBpeak   ', $
ystyle=1, xstyle=1, ytitle='(w1-v)!BB!Lpeak', position=[x1,y1,x2,y2], linestyle=0, color='black'
oplot, feredmags[4,*,febpeak,2]-feredmags[5,*,febpeak,2], feredmags[2,*,febpeak,2]-feredmags[5,*,febpeak,2], linestyle=1
oplot, feredmags[4,*,febpeak,1]-feredmags[5,*,febpeak,1], feredmags[2,*,febpeak,1]-feredmags[5,*,febpeak,1], linestyle=2
oplot, feredmags[4,*,febpeak,0]-feredmags[5,*,febpeak,0], feredmags[2,*,febpeak,0]-feredmags[5,*,febpeak,0], linestyle=3

;;;;;;;;;;;;;;;;;
oploterror, host.bpeakappmag_array[morenormal,4]-host.bpeakappmag_array[morenormal,5], host.bpeakappmag_array[morenormal,2]-host.bpeakappmag_array[morenormal,5], sqrt(host.bpeakappmagerr_array[morenormal,4]^2.0+host.bpeakappmagerr_array[morenormal,5]^2.0), sqrt(host.bpeakappmagerr_array[2,morenormal]^2.0+host.bpeakappmagerr_array[5,morenormal]^2.0), symsize=0.3, psym=15

for i=0,n_elements(morenormal)-1 do  cgplots, host.bpeakappmag_array[morenormal[i],4]-host.bpeakappmag_array[morenormal[i],5], host.bpeakappmag_array[morenormal[i],2]-host.bpeakappmag_array[morenormal[i],5], psym=16, color=morecolortable[i], symsize=1.0



oploterror, host.bpeakappmag_array[morenormal,4]-host.bpeakappmag_array[morenormal,5], host.bpeakappmag_array[morenormal,2]-host.bpeakappmag_array[morenormal,5], sqrt(host.bpeakappmagerr_array[morenormal,1]^2.0+host.bpeakappmagerr_array[morenormal,2]^2.0),  sqrt(host.bpeakappmagerr_array[morenormal,2]^2.0+host.bpeakappmagerr_array[morenormal,5]^2.0), psym=16, color='black', symsize=0.1


cgoplot, host.bpeakappmag_array[morenormalhst,4]-host.bpeakappmag_array[morenormalhst,5], host.bpeakappmag_array[morenormalhst,2]-host.bpeakappmag_array[morenormalhst,5], psym=9, symsize=1.5, color='black'

device, /close
SET_PLOT, 'X'
;spawn , 'open bpeak_bvw1v_hst.eps'

;;;;;;;;;;;;;;;;;;;;


xsize = 8.8
wall = 0.04
margin=0.16
a = xsize/8.8 - (margin + wall)
b = a * 2d / (1 + sqrt(5))

ysize = (margin + b + wall + b + wall + b + wall)*8.8
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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

SET_PLOT, 'PS'

device, filename='w1v_bpeak_colors_abs_color.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

xrange=[-0.2,0.2]
nxticks=4

                                              
cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y2+y2-y1,x2,y2+y2+y2-y1-y1], $
xrange=[-0.1,0.3], yrange=[3.0,0.5], ytitle='(w1-v)!BB!Lpeak',  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
ystyle=1, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1), color='black'


cgcolorbar, ncolors=256 , POSITION=[x1,y2+y2+y2-y1-y1, x2, 0.95], range=[-0.1,0.3], /top, charsize=1

print, 'test 3'



cgoplot, color='black', feredmags[4,*,febpeak,3]-feredmags[5,*,febpeak,3], feredmags[2,*,febpeak,3]-feredmags[5,*,febpeak,3], linestyle=0
cgoplot, color='black', feredmags[4,*,febpeak,2]-feredmags[5,*,febpeak,2], feredmags[2,*,febpeak,2]-feredmags[5,*,febpeak,2], linestyle=1
cgoplot, color='black', feredmags[4,*,febpeak,1]-feredmags[5,*,febpeak,1], feredmags[2,*,febpeak,1]-feredmags[5,*,febpeak,1], linestyle=2
cgoplot, color='black', feredmags[4,*,febpeak,0]-feredmags[5,*,febpeak,0], feredmags[2,*,febpeak,0]-feredmags[5,*,febpeak,0], linestyle=3

; could maybe replace with SN2017erp
;
;restore, 'SN2005cf_redbolmags161.sav'
; filters, ebv, epochs, laws
;  extinction laws model=0 rv=3.1, 1 1.7  , model 2 smc, model 3 CSLMC
; plot, epoch, feredmags[1,0,*,0]

;oplot, cfredmags[4,*,0,3]-cfredmags[5,*,0,3], cfredmags[2,*,0,3]-cfredmags[5,*,0,3], linestyle=0
;oplot, cfredmags[4,*,0,2]-cfredmags[5,*,0,2], cfredmags[2,*,0,2]-cfredmags[5,*,0,2], linestyle=1
;oplot, cfredmags[4,*,0,1]-cfredmags[5,*,0,1], cfredmags[2,*,0,1]-cfredmags[5,*,0,1], linestyle=2
;oplot, cfredmags[4,*,0,0]-cfredmags[5,*,0,0], cfredmags[2,*,0,0]-cfredmags[5,*,0,0], linestyle=3



for i=0,n_elements(normal)-1 do if host.snname_array[normal[i]] ne 'SN2008Q' and  host.snname_array[normal[i]] ne 'SN2008hv'  and  host.snname_array[normal[i]] ne  'SN2011by'  and  host.snname_array[normal[i]] ne  'SN2011fe'  and  host.snname_array[normal[i]] ne  'SN2011ia' then cgplots, host.bpeakappmag_array[normal[i],4]-host.bpeakappmag_array[normal[i],5], host.bpeakappmag_array[normal[i],2]-host.bpeakappmag_array[normal[i],5], psym=16, color=colortable[i], symsize=1.0


index11fe=where(host.snname_array eq 'SN2011fe')
index11ia=where(host.snname_array eq 'SN2011ia')
index11by=where(host.snname_array eq 'SN2011by')
index08Q=where(host.snname_array eq 'SN2008Q')
index08hv=where(host.snname_array eq 'SN2008hv')

cgplots, host.bpeakappmag_array[4,index11ia]-host.bpeakappmag_array[5,index11ia], host.bpeakappmag_array[2,index11ia]-host.bpeakappmag_array[5,index11ia], psym=34, color=colortable[index11ia], symsize=1.0
cgplots, host.bpeakappmag_array[4,index11fe]-host.bpeakappmag_array[5,index11fe], host.bpeakappmag_array[2,index11fe]-host.bpeakappmag_array[5,index11fe], psym=46, color=colortable[index11fe], symsize=1.5
cgplots, host.bpeakappmag_array[4,index11by]-host.bpeakappmag_array[5,index11by], host.bpeakappmag_array[2,index11by]-host.bpeakappmag_array[5,index11by], psym=17, color=colortable[index11by], symsize=1.0
cgplots, host.bpeakappmag_array[4,index08hv]-host.bpeakappmag_array[5,index08hv], host.bpeakappmag_array[2,index08hv]-host.bpeakappmag_array[5,index08hv], psym=15, color=colortable[index08hv], symsize=1.0
cgplots, host.bpeakappmag_array[4,index08Q]-host.bpeakappmag_array[5,index08Q], host.bpeakappmag_array[2,index08Q]-host.bpeakappmag_array[5,index08Q], psym=18, color=colortable[index08Q], symsize=1.0

al_legend, ['SN2008Q','SN2008hv','SN2011by','SN2011fe','SN2011ia'], color=[colortable[index08Q],colortable[index08hv],colortable[index11by],colortable[index11fe],colortable[index11ia]], psym=[18,15,17,46,34], position=[-0.1, 1.8], charsize=0.6, box=0


 oploterror, host.bpeakappmag_array[*,4]-host.bpeakappmag_array[*,5], host.bpeakappmag_array[*,2]-host.bpeakappmag_array[*,5], sqrt(host.bpeakappmagerr_array[*,4]^2.0+host.bpeakappmagerr_array[*,5]^2.0), sqrt(host.bpeakappmagerr_array[*,2]^2.0+host.bpeakappmagerr_array[*,5]^2.0), psym=16, color='black', symsize=0.1


;;;;;;;;;;;;;;;;;;;;;;
cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y2,x2,y2+y2-y1], $
 ytitle='w1!BB!Lpeak !N- $\mu$', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[-15.0,-19], ystyle=1, xrange=[-0.1,0.3], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=4, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1), color='black'


cgoplot, color='black', feredmags[4,*,febpeak,3]-feredmags[5,*,febpeak,3], feredmags[2,*,febpeak,3]-fedm, linestyle=0
cgoplot, color='black', feredmags[4,*,febpeak,2]-feredmags[5,*,febpeak,2], feredmags[2,*,febpeak,2]-fedm, linestyle=1
cgoplot, color='black', feredmags[4,*,febpeak,1]-feredmags[5,*,febpeak,1], feredmags[2,*,febpeak,1]-fedm, linestyle=2
cgoplot, color='black', feredmags[4,*,febpeak,0]-feredmags[5,*,febpeak,0], feredmags[2,*,febpeak,0]-fedm, linestyle=3


al_legend, ['R!LV!N=3.1','R!LV!N=1.7','SMC', 'CSLMC'], linestyle=[0,1,2,3], position=[0.1,-18.9], charsize=0.8, box=0, background='white'

;oploterror, host.bpeakappmag_array[4,w1vsne]-host.bpeakappmag_array[5,w1vsne], host.bpeakappmag_array[2,w1vsne]-host.dm_best_array[normal[w1vsne]], sqrt(host.bpeakappmagerr_array[4,w1vsne]^2.0+host.bpeakappmagerr_array[5,w1vsne]^2.0), sqrt(host.dm_best_err_array[normal[w1vsne]]^2.0+host.bpeakappmagerr_array[5,w1vsne]^2.0), psym=15

for i=0,n_elements(normal)-1 do if host.snname_array[normal[i]] ne 'SN2008Q' and  host.snname_array[normal[i]] ne 'SN2008hv'  and  host.snname_array[normal[i]] ne  'SN2011by'  and  host.snname_array[normal[i]] ne  'SN2011fe'  and  host.snname_array[normal[i]] ne  'SN2011ia' then  cgplots, host.bpeakappmag_array[normal[i],4]-host.bpeakappmag_array[normal[i],5], host.bpeakappmag_array[normal[i],2]-host.dm_best_array[normal[i]], psym=16, color=colortable[i], symsize=1.0

cgplots, host.bpeakappmag_array[4,index11ia]-host.bpeakappmag_array[5,index11ia], host.bpeakappmag_array[2,index11ia]-host.dm_best_array[[index11ia]], psym=34, color=colortable[index11ia], symsize=1.0
cgplots, host.bpeakappmag_array[4,index11fe]-host.bpeakappmag_array[5,index11fe], host.bpeakappmag_array[2,index11fe]-host.dm_best_array[[index11fe]], psym=46, color=colortable[index11fe], symsize=1.5
cgplots, host.bpeakappmag_array[4,index11by]-host.bpeakappmag_array[5,index11by], host.bpeakappmag_array[2,index11by]-host.dm_best_array[[index11by]], psym=17, color=colortable[index11by], symsize=1.0
cgplots, host.bpeakappmag_array[4,index08hv]-host.bpeakappmag_array[5,index08hv], host.bpeakappmag_array[2,index08hv]-host.dm_best_array[[index08hv]], psym=15, color=colortable[index08hv], symsize=1.0
cgplots, host.bpeakappmag_array[4,index08Q]-host.bpeakappmag_array[5,index08Q], host.bpeakappmag_array[2,index08Q]-host.dm_best_array[[index08Q]], psym=18, color=colortable[index08Q], symsize=1.0

oploterror, host.bpeakappmag_array[*,4]-host.bpeakappmag_array[*,5], host.bpeakappmag_array[*,2]-host.dm_best_array[normal[*]], sqrt(host.bpeakappmagerr_array[*,4]^2.0+host.bpeakappmagerr_array[*,5]^2.0), sqrt(host.dm_best_err_array[normal[*]]^2.0+host.bpeakappmagerr_array[*,2]^2.0), psym=16, color='black', symsize=0.1



cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y1,x2,y2], $
xtitle='!A(b-v)!NBpeak       ',   ytitle='v!BB!Lpeak !N- $\mu$', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[-17.0,-20], ystyle=1, xrange=[-0.1,0.3], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=3, ytickv=ytickvalues, color='black'



cgoplot, color='black', feredmags[4,*,febpeak,3]-feredmags[5,*,febpeak,3], feredmags[5,*,febpeak,3]-fedm, linestyle=0
cgoplot, color='black', feredmags[4,*,febpeak,2]-feredmags[5,*,febpeak,2], feredmags[5,*,febpeak,2]-fedm, linestyle=1
cgoplot, color='black', feredmags[4,*,febpeak,1]-feredmags[5,*,febpeak,1], feredmags[5,*,febpeak,1]-fedm, linestyle=2
cgoplot, color='black', feredmags[4,*,febpeak,0]-feredmags[5,*,febpeak,0], feredmags[5,*,febpeak,0]-fedm, linestyle=3


;oploterror, host.bpeakappmag_array[4,w1vsne]-host.bpeakappmag_array[5,w1vsne], host.bpeakappmag_array[5,w1vsne]-host.dm_best_array[normal[w1vsne]], sqrt(host.bpeakappmagerr_array[4,w1vsne]^2.0+host.bpeakappmagerr_array[5,w1vsne]^2.0), sqrt(host.dm_best_err_array[normal[w1vsne]]^2.0+host.bpeakappmagerr_array[5,w1vsne]^2.0), psym=16

cgplots, host.bpeakappmag_array[4,index11ia]-host.bpeakappmag_array[5,index11ia], host.bpeakappmag_array[5,index11ia]-host.dm_best_array[[index11ia]], psym=34, color=colortable[index11ia], symsize=1.0
cgplots, host.bpeakappmag_array[4,index11fe]-host.bpeakappmag_array[5,index11fe], host.bpeakappmag_array[5,index11fe]-host.dm_best_array[[index11fe]], psym=46, color=colortable[index11fe], symsize=1.5
cgplots, host.bpeakappmag_array[4,index11by]-host.bpeakappmag_array[5,index11by], host.bpeakappmag_array[5,index11by]-host.dm_best_array[[index11by]], psym=17, color=colortable[index11by], symsize=1.0
cgplots, host.bpeakappmag_array[4,index08hv]-host.bpeakappmag_array[5,index08hv], host.bpeakappmag_array[5,index08hv]-host.dm_best_array[[index08hv]], psym=15, color=colortable[index08hv], symsize=1.0
cgplots, host.bpeakappmag_array[4,index08Q]-host.bpeakappmag_array[5,index08Q], host.bpeakappmag_array[5,index08Q]-host.dm_best_array[[index08Q]], psym=18, color=colortable[index08Q], symsize=1.0

for i=0,n_elements(normal)-1 do if host.snname_array[normal[i]] ne 'SN2008Q' and  host.snname_array[normal[i]] ne 'SN2008hv'  and  host.snname_array[normal[i]] ne  'SN2011by'  and  host.snname_array[normal[i]] ne  'SN2011fe'  and  host.snname_array[normal[i]] ne  'SN2011ia' then  cgplots, host.bpeakappmag_array[normal[i],4]-host.bpeakappmag_array[normal[i],5], host.bpeakappmag_array[normal[i],5]-host.dm_best_array[normal[i]], psym=16, color=colortable[i], symsize=1.0



oploterror, host.bpeakappmag_array[*,4]-host.bpeakappmag_array[*,5], host.bpeakappmag_array[*,5]-host.dm_best_array[normal[*]], sqrt(host.bpeakappmagerr_array[*,4]^2.0+host.bpeakappmagerr_array[*,5]^2.0), sqrt(host.dm_best_err_array[normal[*]]^2.0+host.bpeakappmagerr_array[*,5]^2.0), psym=16, color='black', symsize=0.1


device, /close
SET_PLOT, 'X'

zcolors=colortable
;scatter_surface, host.dm15_array[normal,4], host.bpeakappmag_array[normal,4]-host.bpeakappmag_array[normal,5], host.bpeakappmag_array[normal,5]-host.dm_best_array[normal], zcolors=colortable, xtitle='$\Delta$M$_{15}$(B)', ytitle='B-V', ztitle='m$_V$-$\mu$'

;spawn, 'open w1v_bpeak_colors_abs_color.eps'


;for n=0,n_elements(normal)-1 do print, host.snname_array[normal[n]], host.bpeakappmag_array[5,n]-host.dm_best_array[normal[n]], sqrt(host.dm_best_err_array[normal[n]]^2.0+host.bpeakappmagerr_array[5,n]^2.0)

;;;;;;;;;;;;;;;;;;;;;
SET_PLOT, 'PS'

device, filename='w1v_bpeak_colors_abs_dmb15_nocolorbar.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize, /color

xrange=[-0.2,0.2]
nxticks=4
                                                  
cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y2+y2-y1,x2,y2+y2+y2-y1-y1], $
xrange=[dmb15low,dmb15high], yrange=[3.0,0.5], ytitle='(w1 - v)!BB!Lpeak',  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
ystyle=1, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1), color='black'


for i=0,n_elements(normal)-1 do if host.snname_array[normal[i]] ne 'SN2008Q' and  host.snname_array[normal[i]] ne 'SN2008hv'  and  host.snname_array[normal[i]] ne  'SN2011by'  and  host.snname_array[normal[i]] ne  'SN2011fe'  and  host.snname_array[normal[i]] ne  'SN2011ia' then  cgplots, host.dm15_array[normal[i],4], host.bpeakappmag_array[normal[i],2]-host.bpeakappmag_array[normal[i],5], psym=16, color=colortable[i], symsize=1.0


cgplots, host.dm15_array[[index11ia],4], host.bpeakappmag_array[index11ia,2]-host.bpeakappmag_array[index11ia,5], psym=34, color=colortable[index11ia], symsize=1.0
cgplots, host.dm15_array[[index11fe],4], host.bpeakappmag_array[index11fe,2]-host.bpeakappmag_array[index11fe,5], psym=46, color=colortable[index11fe], symsize=1.5
cgplots, host.dm15_array[[index11by],4], host.bpeakappmag_array[index11by,2]-host.bpeakappmag_array[index11by,5], psym=17, color=colortable[index11by], symsize=1.0
cgplots, host.dm15_array[[index08hv],4], host.bpeakappmag_array[index08hv,2]-host.bpeakappmag_array[index08hv,5], psym=15, color=colortable[index08hv], symsize=1.0
cgplots, host.dm15_array[[index08Q],4],  host.bpeakappmag_array[index08Q,2] -host.bpeakappmag_array[index08Q,5],  psym=18, color=colortable[index08Q], symsize=1.0

oploterror, host.dm15_array[normal,4], host.bpeakappmag_array[*,2]-host.bpeakappmag_array[*,5], host.dm15err_array[normal,4], sqrt(host.bpeakappmagerr_array[*,2]^2.0+host.bpeakappmagerr_array[*,5]^2.0), psym=16, symsize=0.1, color='black'

;;;;;;;;;;;;;;;;;;;;;;
cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y2,x2,y2+y2-y1], $
 ytitle='w1!BB!Lpeak !N- $\mu$', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[-15.0,-19], ystyle=1, xrange=[dmb15low,dmb15high], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=4, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1), color='black'

for i=0,n_elements(normal)-1 do if host.snname_array[normal[i]] ne 'SN2008Q' and  host.snname_array[normal[i]] ne 'SN2008hv'  and  host.snname_array[normal[i]] ne  'SN2011by'  and  host.snname_array[normal[i]] ne  'SN2011fe'  and  host.snname_array[normal[i]] ne  'SN2011ia' then  cgplots, host.dm15_array[normal[i],4], host.bpeakappmag_array[normal[i],2]-host.dm_best_array[normal[i]], psym=16, color=colortable[i], symsize=1.0

cgplots, host.dm15_array[[index11ia],4], host.bpeakappmag_array[2,index11ia]-host.dm_best_array[[index11ia]], psym=34, color=colortable[index11ia], symsize=1.0
cgplots, host.dm15_array[[index11fe],4], host.bpeakappmag_array[2,index11fe]-host.dm_best_array[[index11fe]], psym=46, color=colortable[index11fe], symsize=1.5
cgplots, host.dm15_array[[index11by],4], host.bpeakappmag_array[2,index11by]-host.dm_best_array[[index11by]], psym=17, color=colortable[index11by], symsize=1.0
cgplots, host.dm15_array[[index08hv],4], host.bpeakappmag_array[2,index08hv]-host.dm_best_array[[index08hv]], psym=15, color=colortable[index08hv], symsize=1.0
cgplots, host.dm15_array[[index08Q],4], host.bpeakappmag_array[2,index08Q]-host.dm_best_array[[index08Q]], psym=18, color=colortable[index08Q], symsize=1.0

oploterror, host.dm15_array[normal,4], host.bpeakappmag_array[*,2]-host.dm_best_array[normal[*]], host.dm15err_array[normal,4], sqrt(host.dm_best_err_array[normal[*]]^2.0+host.bpeakappmagerr_array[*,5]^2.0), psym=16, symsize=0.1, color='black'

cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y1,x2,y2], $
xtitle='!A$\Delta$ M!N15!A(B)',   ytitle='v!BB!Lpeak !N- $\mu$',  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[-17.0,-20], ystyle=1, xrange=[dmb15low,dmb15high], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=3, ytickv=ytickvalues, color='black'

for i=0,n_elements(normal)-1 do if host.snname_array[normal[i]] ne 'SN2008Q' and  host.snname_array[normal[i]] ne 'SN2008hv'  and  host.snname_array[normal[i]] ne  'SN2011by'  and  host.snname_array[normal[i]] ne  'SN2011fe'  and  host.snname_array[normal[i]] ne  'SN2011ia' then  cgplots, host.dm15_array[normal[i],4], host.bpeakappmag_array[normal[i],5]-host.dm_best_array[normal[i]], psym=16, color=colortable[i], symsize=1.0

cgplots, host.dm15_array[[index11ia],4], host.bpeakappmag_array[5,index11ia]-host.dm_best_array[[index11ia]], psym=34, color=colortable[index11ia], symsize=1.0
cgplots, host.dm15_array[[index11fe],4], host.bpeakappmag_array[5,index11fe]-host.dm_best_array[[index11fe]], psym=46, color=colortable[index11fe], symsize=1.5
cgplots, host.dm15_array[[index11by],4], host.bpeakappmag_array[5,index11by]-host.dm_best_array[[index11by]], psym=17, color=colortable[index11by], symsize=1.0
cgplots, host.dm15_array[[index08hv],4], host.bpeakappmag_array[5,index08hv]-host.dm_best_array[[index08hv]], psym=15, color=colortable[index08hv], symsize=1.0
cgplots, host.dm15_array[[index08Q],4],  host.bpeakappmag_array[5,index08Q] -host.dm_best_array[[index08Q]],  psym=18, color=colortable[index08Q], symsize=1.0

oploterror, host.dm15_array[normal,4],  host.bpeakappmag_array[*,5]-host.dm_best_array[normal[*]], host.dm15err_array[normal,4],  sqrt(host.dm_best_err_array[normal[*]]^2.0+host.bpeakappmagerr_array[*,5]^2.0), psym=16, symsize=0.1, color='black'

device, /close
SET_PLOT, 'X'

;spawn, 'open w1v_bpeak_colors_abs_dmb15_nocolorbar.eps'




SET_PLOT, 'PS'

device, filename='w1v_bpeak_colors_abs_dmb15.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize, /color

nxticks=4
                                                  
cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y2+y2-y1,x2,y2+y2+y2-y1-y1], $
xrange=[dmb15low,dmb15high], yrange=[3.0,0.5], ytitle='(w1-v)!BB!Lpeak',$
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
ystyle=1, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1), color=black


for i=0,n_elements(normal)-1 do cgplots, host.dm15_array[normal[i],4], host.bpeakappmag_array[normal[i],2]-host.bpeakappmag_array[normal[i],5], psym=16, color=colortable[i], symsize=1.0


oploterror, host.dm15_array[normal,4], host.bpeakappmag_array[normal,2]-host.bpeakappmag_array[normal,5], host.dm15err_array[normal,4], sqrt(host.bpeakappmagerr_array[normal,2]^2.0+host.bpeakappmagerr_array[normal,5]^2.0), psym=16, symsize=0.1

;cgLOADCT, 33, NCOLORS=100
cgcolorbar, ncolors=256 , POSITION=[x1,y2+y2+y2-y1-y1, x2, 0.95], range=[-0.1,0.3], /top, charsize=1

;;;;;;;;;;;;;;;;;;;;;;
cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y2,x2,y2+y2-y1], $
 ytitle='w1!BB!Lpeak !N- $\mu$', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[-16.0,-19], ystyle=1, xrange=[dmb15low,dmb15high], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=3, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1), color=black


for i=0,n_elements(normal)-1 do cgplots, host.dm15_array[normal[i],4], host.bpeakappmag_array[normal[i],2]-host.dm_best_array[normal[i]], psym=16, color=colortable[i], symsize=1.0



oploterror, host.dm15_array[normal,4], host.bpeakappmag_array[normal,2]-host.dm_best_array[normal], host.dm15err_array[normal,4], sqrt(host.dm_best_err_array[normal]^2.0+host.bpeakappmagerr_array[normal,5]^2.0), psym=16, symsize=0.1


cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y1,x2,y2], $
xtitle='!A$\Delta$ M!N15!A(B)',   ytitle='v!BB!Lpeak !N- $\mu$', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[-17.0,-20], ystyle=1, xrange=[dmb15low,dmb15high], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=3, ytickv=ytickvalues, color=black

for i=0,n_elements(normal)-1 do cgplots, host.dm15_array[normal[i],4], host.bpeakappmag_array[normal[i],5]-host.dm_best_array[normal[i]], psym=16, color=colortable[i], symsize=1.0


oploterror, host.dm15_array[normal,4],  host.bpeakappmag_array[normal,5]-host.dm_best_array[normal], host.dm15err_array[normal,4],  sqrt(host.dm_best_err_array[normal]^2.0+host.bpeakappmagerr_array[normal,5]^2.0), psym=16, symsize=0.1


device, /close
SET_PLOT, 'X'


;spawn, 'open w1v_bpeak_colors_abs_dmb15.eps'


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



SET_PLOT, 'PS'

device, filename='w1v_peak_colors_abs_dmb15_nocolorbar.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize, /color

xrange=[-0.2,0.2]
nxticks=4
                                                  
cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y2+y2-y1,x2,y2+y2+y2-y1-y1], $
xrange=[dmb15low,dmb15high], yrange=[3.0,0.0], ytitle='w1!Bpeak !N-!Nv!Bpeak !N',  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
ystyle=1, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1), color='black'


for i=0,n_elements(normal)-1 do if host.snname_array[normal[i]] ne 'SN2008Q' and  host.snname_array[normal[i]] ne 'SN2008hv'  and  host.snname_array[normal[i]] ne  'SN2011by'  and  host.snname_array[normal[i]] ne  'SN2011fe'  and  host.snname_array[normal[i]] ne  'SN2011ia' then  cgplots, host.dm15_array[normal[i],4], host.appmag_array[normal[i],2]-host.appmag_array[normal[i],5], psym=16, color=colortable[i], symsize=1.0


cgplots, host.dm15_array[[index11ia],4], host.appmag_array[index11ia,2]-host.appmag_array[index11ia,5], psym=34, color=colortable[index11ia], symsize=1.0
cgplots, host.dm15_array[[index11fe],4], host.appmag_array[index11fe,2]-host.appmag_array[index11fe,5], psym=46, color=colortable[index11fe], symsize=1.5
cgplots, host.dm15_array[[index11by],4], host.appmag_array[index11by,2]-host.appmag_array[index11by,5], psym=17, color=colortable[index11by], symsize=1.0
cgplots, host.dm15_array[[index08hv],4], host.appmag_array[index08hv,2]-host.appmag_array[index08hv,5], psym=15, color=colortable[index08hv], symsize=1.0
cgplots, host.dm15_array[[index08Q],4],  host.appmag_array[index08Q,2] -host.appmag_array[index08Q,5],  psym=18, color=colortable[index08Q],  symsize=1.0

oploterror, host.dm15_array[normal,4], host.appmag_array[*,2]-host.appmag_array[*,5], host.dm15err_array[normal,4], sqrt(host.appmagerr_array[*,2]^2.0+host.appmagerr_array[*,5]^2.0), psym=16, symsize=0.1, color='black'

;;;;;;;;;;;;;;;;;;;;;;
cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y2,x2,y2+y2-y1], $
 ytitle='w1!Bpeak !N- $\mu$', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[-15.0,-20], ystyle=1, xrange=[dmb15low,dmb15high], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=4, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1), color='black'

for i=0,n_elements(normal)-1 do if host.snname_array[normal[i]] ne 'SN2008Q' and  host.snname_array[normal[i]] ne 'SN2008hv'  and  host.snname_array[normal[i]] ne  'SN2011by'  and  host.snname_array[normal[i]] ne  'SN2011fe'  and  host.snname_array[normal[i]] ne  'SN2011ia' then  cgplots, host.dm15_array[normal[i],4], host.appmag_array[normal[i],2]-host.dm_best_array[normal[i]], psym=16, color=colortable[i], symsize=1.0

cgplots, host.dm15_array[[index11ia],4], host.appmag_array[2,index11ia]-host.dm_best_array[[index11ia]], psym=34, color=colortable[index11ia], symsize=1.0
cgplots, host.dm15_array[[index11fe],4], host.appmag_array[2,index11fe]-host.dm_best_array[[index11fe]], psym=46, color=colortable[index11fe], symsize=1.5
cgplots, host.dm15_array[[index11by],4], host.appmag_array[2,index11by]-host.dm_best_array[[index11by]], psym=17, color=colortable[index11by], symsize=1.0
cgplots, host.dm15_array[[index08hv],4], host.appmag_array[2,index08hv]-host.dm_best_array[[index08hv]], psym=15, color=colortable[index08hv], symsize=1.0
cgplots, host.dm15_array[[index08Q],4],  host.appmag_array[2,index08Q]-host.dm_best_array[[index08Q]], psym=18, color=colortable[index08Q], symsize=1.0

oploterror, host.dm15_array[normal,4], host.appmag_array[*,2]-host.dm_best_array[normal[*]], host.dm15err_array[normal,4], sqrt(host.dm_best_err_array[normal[*]]^2.0+host.appmagerr_array[*,5]^2.0), psym=16, symsize=0.1, color='black'

cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y1,x2,y2], $
xtitle='!A$\Delta$ M!N15!A(B)',   ytitle='v!Bpeak !N- $\mu$',  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[-17.0,-20], ystyle=1, xrange=[dmb15low,dmb15high], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=3, ytickv=ytickvalues, color='black'

for i=0,n_elements(normal)-1 do if host.snname_array[normal[i]] ne 'SN2008Q' and  host.snname_array[normal[i]] ne 'SN2008hv'  and  host.snname_array[normal[i]] ne  'SN2011by'  and  host.snname_array[normal[i]] ne  'SN2011fe'  and  host.snname_array[normal[i]] ne  'SN2011ia' then  cgplots, host.dm15_array[normal[i],4], host.appmag_array[normal[i],5]-host.dm_best_array[normal[i]], psym=16, color=colortable[i], symsize=1.0

cgplots, host.dm15_array[[index11ia],4], host.appmag_array[5,index11ia]-host.dm_best_array[[index11ia]], psym=34, color=colortable[index11ia], symsize=1.0
cgplots, host.dm15_array[[index11fe],4], host.appmag_array[5,index11fe]-host.dm_best_array[[index11fe]], psym=46, color=colortable[index11fe], symsize=1.5
cgplots, host.dm15_array[[index11by],4], host.appmag_array[5,index11by]-host.dm_best_array[[index11by]], psym=17, color=colortable[index11by], symsize=1.0
cgplots, host.dm15_array[[index08hv],4], host.appmag_array[5,index08hv]-host.dm_best_array[[index08hv]], psym=15, color=colortable[index08hv], symsize=1.0
cgplots, host.dm15_array[[index08Q],4],  host.appmag_array[5,index08Q] -host.dm_best_array[[index08Q]],  psym=18, color=colortable[index08Q], symsize=1.0

oploterror, host.dm15_array[normal,4],  host.appmag_array[*,5]-host.dm_best_array[normal[*]], host.dm15err_array[normal,4],  sqrt(host.dm_best_err_array[normal[*]]^2.0+host.appmagerr_array[*,5]^2.0), psym=16, symsize=0.1, color='black'

device, /close
SET_PLOT, 'X'

;spawn, 'open w1v_peak_colors_abs_dmb15_nocolorbar.eps'




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
fedm=29.04
;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

xsize = 8.8
wall = 0.04
margin=0.18
a = xsize/8.8 - (margin + wall)
b = a * 2d / (1 + sqrt(5))

ysize = (margin + b + wall + b + wall + b + wall)*8.8
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


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

SET_PLOT, 'PS'

device, filename='colors_dm15b.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

xrange=[-0.2,0.2]
nxticks=4

                                              
cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y2+y2-y1,x2,y2+y2+y2-y1-y1], $
xrange=[dmb15low,dmb15high], yrange=[bvhigh, bvlow], ytitle='(B-V)!BB!Lpeak',  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
ystyle=1, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1), color='black'


cgcolorbar, ncolors=256 , POSITION=[0.9,y2+y2-y1, x2, y2+y2+y2-y1-y1], range=[0.3,-0.15], /vertical, /reverse, /right, charsize=1, ticknames=[' ', ' ', ' ', ' ', ' ']



for i=0,n_elements(normal)-1 do  cgplots, host.dm15_array[normal[i],4], host.bpeakappmag_array[normal[i],4]-host.bpeakappmag_array[normal[i],5], psym=16, color=colortable[i], symsize=1.0


 oploterror, host.dm15_array[normal,4], host.bpeakappmag_array[normal,4]-host.bpeakappmag_array[normal,5], host.dm15err_array[normal,4], sqrt(host.bpeakappmagerr_array[normal,4]^2.0+host.bpeakappmagerr_array[normal,5]^2.0), psym=16, color='black', symsize=0.1


;;;;;;;;;;;;;;;;;;;;;;
cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y2,x2,y2+y2-y1], $
 ytitle='w1-b!BB!Lpeak ', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[3,0], ystyle=1, xrange=[dmb15low,dmb15high], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1), color='black'



for i=0,n_elements(normal)-1 do cgplots, host.dm15_array[normal[i],4], host.bpeakappmag_array[normal[i],2]-host.bpeakappmag_array[normal[i],4], psym=16, color=colortable[i], symsize=1.0


oploterror, host.dm15_array[normal,4], host.bpeakappmag_array[normal,2]-host.bpeakappmag_array[normal,4], host.dm15err_array[normal,4],  sqrt(host.bpeakappmagerr_array[normal,2]^2.0+host.bpeakappmagerr_array[normal,4]^2.0), psym=16, color='black', symsize=0.1



cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y1,x2,y2], $
xtitle='!A$\Delta$ M!N15!A(B)',   ytitle='m2-w1!BB!Lpeak', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[4,0], ystyle=1, xrange=[dmb15low,dmb15high], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, ytickv=ytickvalues, color='black'


for i=0,n_elements(normal)-1 do  cgplots, host.dm15_array[normal[i],4], host.bpeakappmag_array[normal[i],1]-host.bpeakappmag_array[normal[i],2], psym=16, color=colortable[i], symsize=1.0



oploterror, host.dm15_array[normal,4], host.bpeakappmag_array[normal,1]-host.bpeakappmag_array[normal,2], host.dm15err_array[normal,4],  sqrt(host.bpeakappmagerr_array[normal,1]^2.0+host.bpeakappmagerr_array[normal,2]^2.0), psym=16, color='black', symsize=0.1


device, /close
SET_PLOT, 'X'


;spawn , 'open colors_dm15b.eps'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


readcol,"uvmodel.data",wave_foley, flux_foley, delta_foley,format='F,F,F',/silent

dmb15_array=[0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7]
n_dmb15=n_elements(dmb15_array)
foleyflux_array=make_array(n_dmb15,n_elements(wave_foley))
foleymag_array=make_array(n_dmb15,6)

for n=0,n_dmb15-1 do foleyflux_array[n,*]=flux_foley+delta_foley*(dmb15_array[n]-1.1)


for n=0,n_dmb15-1 do begin

	pjb_uvotspec_all, [transpose(wave_foley), foleyflux_array[n,*]], mag_array=mag_array

	foleymag_array[n,*]=mag_array[0:5]

endfor



SET_PLOT, 'PS'

device, filename='colors_dm15b_hst.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

xrange=[-0.2,0.2]
nxticks=4

                                              
cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y2+y2-y1,x2,y2+y2+y2-y1-y1], $
xrange=[dmb15low,dmb15high], yrange=[bvhigh,bvlow], ytitle='(B-V)!BB!Lpeak',  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
ystyle=1, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1), color='black'


cgcolorbar, ncolors=256 , POSITION=[0.9,y2+y2-y1, x2, y2+y2+y2-y1-y1], range=[0.3,-0.15], /vertical, /reverse, /right, charsize=1, ticknames=[' ', ' ', ' ', ' ', ' ']



for i=0,n_elements(normal)-1 do  cgplots, host.dm15_array[normal[i],4], host.bpeakappmag_array[normal[i],4]-host.bpeakappmag_array[normal[i],5], psym=16, color=colortable[i], symsize=1.0


 oploterror, host.dm15_array[normal,4], host.bpeakappmag_array[normal,4]-host.bpeakappmag_array[normal,5], host.dm15err_array[normal,4], sqrt(host.bpeakappmagerr_array[normal,4]^2.0+host.bpeakappmagerr_array[normal,5]^2.0), psym=16, color='black', symsize=0.1

if hst[0] ne -1 then cgoplot, host.dm15_array[hst,4], host.bpeakappmag_array[hst,4]-host.bpeakappmag_array[hst,5], psym=9, symsize=1.5, color='black'


;;;;;;;;;;;;;;;;;;;;;;
cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y2,x2,y2+y2-y1], $
 ytitle='uvw1-b!BB!Lpeak ', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[3,0], ystyle=1, xrange=[dmb15low,dmb15high], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1), color='black'



for i=0,n_elements(normal)-1 do cgplots, host.dm15_array[normal[i],4], host.bpeakappmag_array[normal[i],2]-host.bpeakappmag_array[normal[i],4], psym=16, color=colortable[i], symsize=1.0


oploterror, host.dm15_array[normal,4], host.bpeakappmag_array[normal,2]-host.bpeakappmag_array[normal,4], host.dm15err_array[normal,4],  sqrt(host.bpeakappmagerr_array[normal,2]^2.0+host.bpeakappmagerr_array[normal,4]^2.0), psym=16, color='black', symsize=0.1

if hst[0] ne -1 then cgoplot, host.dm15_array[hst,4], host.bpeakappmag_array[hst,2]-host.bpeakappmag_array[hst,4], psym=9, symsize=1.5, color='black'

cgoplot, dmb15_array, foleymag_array[*,2]-foleymag_array[*,4], thick=8


cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y1,x2,y2], $
xtitle='!A$\Delta$ M!N15!A(B)',   ytitle='uvm2-uvw1!BB!Lpeak', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[4,0], ystyle=1, xrange=[dmb15low,dmb15high], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, ytickv=ytickvalues, color='black'


for i=0,n_elements(normal)-1 do  cgplots, host.dm15_array[normal[i],4], host.bpeakappmag_array[normal[i],1]-host.bpeakappmag_array[normal[i],2], psym=16, color=colortable[i], symsize=1.0



oploterror, host.dm15_array[normal,4], host.bpeakappmag_array[normal,1]-host.bpeakappmag_array[normal,2], host.dm15err_array[normal,4],  sqrt(host.bpeakappmagerr_array[normal,1]^2.0+host.bpeakappmagerr_array[normal,2]^2.0), psym=16, color='black', symsize=0.1


if hst[0] ne -1 then cgoplot, host.dm15_array[hst,4], host.bpeakappmag_array[hst,1]-host.bpeakappmag_array[hst,2], psym=9, symsize=1.5, color='black'


cgoplot, dmb15_array, foleymag_array[*,1]-foleymag_array[*,2], thick=8


device, /close
SET_PLOT, 'X'


;spawn , 'open colors_dm15b_hst.eps'

;


SET_PLOT, 'PS'

device, filename='colors_dm15b_hst_bvbottom.eps', /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

xrange=[-0.2,0.2]
nxticks=4

                                              
cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y2+y2-y1,x2,y2+y2+y2-y1-y1], $
xrange=[dmb15low,dmb15high], yrange=[4,0], ytitle='uvm2-uvw1!BB!Lpeak',  $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
ystyle=1, xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, yticks=nyticks, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1), color='black'


for i=0,n_elements(normal)-1 do  cgplots, host.dm15_array[normal[i],4], host.bpeakappmag_array[normal[i],1]-host.bpeakappmag_array[normal[i],2], psym=16, color=colortable[i], symsize=1.0



oploterror, host.dm15_array[normal,4], host.bpeakappmag_array[normal,1]-host.bpeakappmag_array[normal,2], host.dm15err_array[normal,4],  sqrt(host.bpeakappmagerr_array[normal,1]^2.0+host.bpeakappmagerr_array[normal,2]^2.0), psym=16, color='black', symsize=0.1


if hst[0] ne -1 then cgoplot, host.dm15_array[hst,4], host.bpeakappmag_array[hst,1]-host.bpeakappmag_array[hst,2], psym=9, symsize=1.5, color='black'


cgoplot, dmb15_array, foleymag_array[*,1]-foleymag_array[*,2], thick=8



;;;;;;;;;;;;;;;;;;;;;;
cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y2,x2,y2+y2-y1], $
 ytitle='uvw1-B!BB!Lpeak ', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[3,0], ystyle=1, xrange=[dmb15low,dmb15high], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, ytickv=ytickvalues, xtickname=replicate(' ',nxticks+1), color='black'



for i=0,n_elements(normal)-1 do cgplots, host.dm15_array[normal[i],4], host.bpeakappmag_array[normal[i],2]-host.bpeakappmag_array[normal[i],4], psym=16, color=colortable[i], symsize=1.0


oploterror, host.dm15_array[normal,4], host.bpeakappmag_array[normal,2]-host.bpeakappmag_array[normal,4], host.dm15err_array[normal,4],  sqrt(host.bpeakappmagerr_array[normal,2]^2.0+host.bpeakappmagerr_array[normal,4]^2.0), psym=16, color='black', symsize=0.1

if hst[0] ne -1 then cgoplot, host.dm15_array[hst,4], host.bpeakappmag_array[hst,2]-host.bpeakappmag_array[hst,4], psym=9, symsize=1.5, color='black'

cgoplot, dmb15_array, foleymag_array[*,2]-foleymag_array[*,4], thick=8


cgplot, charsize=1, xdata, ydata, /nodata, /noerase, position=[x1,y1,x2,y2], $
xtitle='!A$\Delta$ M!N15!A(B)',   ytitle='(B-V)!BB!Lpeak', $
xminor=1, yminor=1, xticklen=xticklen, yticklen=yticklen, $
 yrange=[bvhigh,bvlow], ystyle=1, xrange=[dmb15low,dmb15high], xstyle=1, $
xticks=nxticks, xtickv=xtickvalues, ytickv=ytickvalues, color='black'


; when on top panel  cgcolorbar, ncolors=256 , POSITION=[0.9,y2+y2-y1, x2, y2+y2+y2-y1-y1], range=[0.3,-0.15], /vertical, /reverse, /right, charsize=1, ticknames=[' ', ' ', ' ', ' ', ' ']

cgcolorbar, ncolors=256 , POSITION=[0.9,y1,x2,y2], range=[0.3,-0.15], /vertical, /reverse, /right, charsize=1, ticknames=[' ', ' ', ' ', ' ', ' ']



for i=0,n_elements(normal)-1 do  cgplots, host.dm15_array[normal[i],4], host.bpeakappmag_array[normal[i],4]-host.bpeakappmag_array[normal[i],5], psym=16, color=colortable[i], symsize=1.0


 oploterror, host.dm15_array[normal,4], host.bpeakappmag_array[normal,4]-host.bpeakappmag_array[normal,5], host.dm15err_array[normal,4], sqrt(host.bpeakappmagerr_array[normal,4]^2.0+host.bpeakappmagerr_array[normal,5]^2.0), psym=16, color='black', symsize=0.1

if hst[0] ne -1 then cgoplot, host.dm15_array[hst,4], host.bpeakappmag_array[hst,4]-host.bpeakappmag_array[hst,5], psym=9, symsize=1.5, color='black'



device, /close
SET_PLOT, 'X'


spawn, 'open colors_dm15b_hst_bvbottom.eps'

order=sort(host.bpeakappmag_array[*,2]-host.bpeakappmag_array[*,5])
; for i=0,30 do print, host.snname_array[order[i]], host.bpeakappmag_array[order[i],1]-host.bpeakappmag_array[order[i],2], host.bpeakappmag_array[order[i],2]-host.bpeakappmag_array[order[i],5]
; for i=0,50 do print, host.snname_array[order[i]], host.bpeakappmag_array[order[i],1]-host.bpeakappmag_array[order[i],5], host.bpeakappmag_array[order[i],2]-host.bpeakappmag_array[order[i],5]



print, 'final stop'
stop

end
