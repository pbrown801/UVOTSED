pro redvectors

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


plotfilename = 'bpeak_w1vbv_redvector.eps'

SET_PLOT, 'PS'

device, filename= plotfilename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

restore, 'SN2011fe_redbolmags161.sav'
; filters, ebv, epochs, laws
;  extinction laws model=0 rv=3.1, 1 1.7  , model 2 smc, model 3 CSLMC
; plot, epoch, feredmags[1,0,*,0]
febpeak=where(min(feredmags[4,0,*,0]) eq feredmags[4,0,*,0],count)


cgplot, charsize=1, feredmags[4,*,febpeak,3]-feredmags[5,*,febpeak,3], feredmags[2,*,febpeak,3]-feredmags[5,*,febpeak,3], xrange=[-0.5,1.5], yrange=[-1,5], ystyle=1, xstyle=1, ytitle='(w1-v)!BB!Lpeak', $ 
xtitle=' !A (b-v)!NBpeak   ', $
; double subscripts falling off page
; xtitle='!S!U (b-v) !N B !R!I peak', $
position=[x1,y1,x2,y2], linestyle=0, color='blue'
;  not getting this to work
;xyouts, 0.02, 0.5, '(b-v) !R!I B peak'

cgoplot, feredmags[4,*,febpeak,2]-feredmags[5,*,febpeak,2], feredmags[2,*,febpeak,2]-feredmags[5,*,febpeak,2], linestyle=1, color='blue'
cgoplot, feredmags[4,*,febpeak,1]-feredmags[5,*,febpeak,1], feredmags[2,*,febpeak,1]-feredmags[5,*,febpeak,1], linestyle=2, color='blue'
cgoplot, feredmags[4,*,febpeak,0]-feredmags[5,*,febpeak,0], feredmags[2,*,febpeak,0]-feredmags[5,*,febpeak,0], linestyle=3, color='blue'


restore, 'SN2016ccj_160514_uvopt_obswave_obsflux.dat_redbolmags161z.sav'

cgoplot, redzmags_array[4,*,0,0,0]-redzmags_array[5,*,0,0,0], redzmags_array[2,*,0,0,0]-redzmags_array[5,*,0,0,0], linestyle=1, color='violet'
cgoplot, redzmags_array[4,*,0,1,0]-redzmags_array[5,*,0,1,0], redzmags_array[2,*,0,1,0]-redzmags_array[5,*,0,1,0], linestyle=1, color='violet'
cgoplot, redzmags_array[4,*,0,2,0]-redzmags_array[5,*,0,2,0], redzmags_array[2,*,0,2,0]-redzmags_array[5,*,0,2,0], linestyle=1, color='violet'
cgoplot, redzmags_array[4,*,0,3,0]-redzmags_array[5,*,0,3,0], redzmags_array[2,*,0,3,0]-redzmags_array[5,*,0,3,0], linestyle=1, color='violet'


restore, 'SN2011erp_redbolmags161z.sav'
cgoplot, erpredmagsz[4,*,0,0,0]-erpredmagsz[5,*,0,0,0], erpredmagsz[2,*,0,0,0]-erpredmagsz[5,*,0,0,0], linestyle=1, color='red'
cgoplot, erpredmagsz[4,*,0,1,0]-erpredmagsz[5,*,0,1,0], erpredmagsz[2,*,0,1,0]-erpredmagsz[5,*,0,1,0], linestyle=1, color='red'
cgoplot, erpredmagsz[4,*,0,2,0]-erpredmagsz[5,*,0,2,0], erpredmagsz[2,*,0,2,0]-erpredmagsz[5,*,0,2,0], linestyle=1, color='red'
cgoplot, erpredmagsz[4,*,0,3,0]-erpredmagsz[5,*,0,3,0], erpredmagsz[2,*,0,3,0]-erpredmagsz[5,*,0,3,0], linestyle=1, color='red'



restore, 'SN2017erp_170702_uvopt_obswave_obsflux.dat_redbolmags161z.sav'
cgoplot, redzmags_array[4,*,0,0,0]-redzmags_array[5,*,0,0,0], redzmags_array[2,*,0,0,0]-redzmags_array[5,*,0,0,0], linestyle=1, color='red'
cgoplot, redzmags_array[4,*,0,1,0]-redzmags_array[5,*,0,1,0], redzmags_array[2,*,0,1,0]-redzmags_array[5,*,0,1,0], linestyle=1, color='red'
cgoplot, redzmags_array[4,*,0,2,0]-redzmags_array[5,*,0,2,0], redzmags_array[2,*,0,2,0]-redzmags_array[5,*,0,2,0], linestyle=1, color='red'
cgoplot, redzmags_array[4,*,0,3,0]-redzmags_array[5,*,0,3,0], redzmags_array[2,*,0,3,0]-redzmags_array[5,*,0,3,0], linestyle=1, color='red'


;;;;;;;;;;;;;;;;;
oploterror, host.bpeakappmag_array[*,4]-host.bpeakappmag_array[*,5], host.bpeakappmag_array[*,2]-host.bpeakappmag_array[*,5], sqrt(host.bpeakappmagerr_array[*,4]^2.0+host.bpeakappmagerr_array[*,5]^2.0), sqrt(host.bpeakappmagerr_array[*,2]^2.0+host.bpeakappmagerr_array[*,5]^2.0), symsize=0.3, psym=15


xyouts, host.bpeakappmag_array[where(host.snname_array eq 'SN2014J'),4]-host.bpeakappmag_array[where(host.snname_array eq 'SN2014J'),5]-0.05, host.bpeakappmag_array[where(host.snname_array eq 'SN2014J'),2]-host.bpeakappmag_array[where(host.snname_array eq 'SN2014J'),5] +0.2, 'SN2014J', charsize=0.5


;xyouts, host.bpeakappmag_array[where(host.snname_array eq 'SN2013gs'),4]-host.bpeakappmag_array[where(host.snname_array eq 'SN2013gs'),5]+0.05, host.bpeakappmag_array[where(host.snname_array eq 'SN2013gs'),2]-host.bpeakappmag_array[where(host.snname_array eq 'SN2013gs'),5] - 0.15, 'SN2013gs', charsize=0.5

cgoplot, 15.54-15.37, 17.08-15.37, psym=4, color='red'
xyouts, 15.54-15.37+0.25, 17.08-15.37, 'SN2013gs Zhang et al.', charsize=0.5

xyouts, host.bpeakappmag_array[where(host.snname_array eq 'SN2006X'),4]-host.bpeakappmag_array[where(host.snname_array eq 'SN2006X'),5]-0.15, host.bpeakappmag_array[where(host.snname_array eq 'SN2006X'),2]-host.bpeakappmag_array[where(host.snname_array eq 'SN2006X'),5] +0.28, 'SN2006X', charsize=0.5



device, /close
SET_PLOT, 'X'
spawn, 'open bpeak_w1vbv_redvector.eps'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


plotfilename = 'bpeak_m2vbv_redvector.eps'

SET_PLOT, 'PS'

device, filename= plotfilename, /encapsulated, xsize=xsize, ysize=ysize, $
/tt_font, set_font='Times', font_size=fontsize

restore, 'SN2011fe_redbolmags161.sav'
; filters, ebv, epochs, laws
;  extinction laws model=0 rv=3.1, 1 1.7  , model 2 smc, model 3 CSLMC
; plot, epoch, feredmags[1,0,*,0]
febpeak=where(min(feredmags[4,0,*,0]) eq feredmags[4,0,*,0],count)


cgplot, charsize=1, feredmags[4,*,febpeak,3]-feredmags[5,*,febpeak,3], feredmags[1,*,febpeak,3]-feredmags[5,*,febpeak,3], xrange=[-0.5,1.5], yrange=[0,9], ystyle=1, xstyle=1, ytitle='(m2-v)!BB!Lpeak', $ 
xtitle=' !A (b-v)!NBpeak   ', $
; double subscripts falling off page
; xtitle='!S!U (b-v) !N B !R!I peak', $
position=[x1,y1,x2,y2], linestyle=0, color='blue'
;  not getting this to work
;xyouts, 0.02, 0.5, '(b-v) !R!I B peak'

cgoplot, feredmags[4,*,febpeak,2]-feredmags[5,*,febpeak,2], feredmags[1,*,febpeak,2]-feredmags[5,*,febpeak,2], linestyle=1, color='blue'
cgoplot, feredmags[4,*,febpeak,1]-feredmags[5,*,febpeak,1], feredmags[1,*,febpeak,1]-feredmags[5,*,febpeak,1], linestyle=2, color='blue'
cgoplot, feredmags[4,*,febpeak,0]-feredmags[5,*,febpeak,0], feredmags[1,*,febpeak,0]-feredmags[5,*,febpeak,0], linestyle=3, color='blue'


restore, 'SN2016ccj_160514_uvopt_obswave_obsflux.dat_redbolmags161z.sav'

cgoplot, redzmags_array[4,*,0,0,0]-redzmags_array[5,*,0,0,0], redzmags_array[1,*,0,0,0]-redzmags_array[5,*,0,0,0], linestyle=1, color='violet'
cgoplot, redzmags_array[4,*,0,1,0]-redzmags_array[5,*,0,1,0], redzmags_array[1,*,0,1,0]-redzmags_array[5,*,0,1,0], linestyle=1, color='violet'
cgoplot, redzmags_array[4,*,0,2,0]-redzmags_array[5,*,0,2,0], redzmags_array[1,*,0,2,0]-redzmags_array[5,*,0,2,0], linestyle=1, color='violet'
cgoplot, redzmags_array[4,*,0,3,0]-redzmags_array[5,*,0,3,0], redzmags_array[1,*,0,3,0]-redzmags_array[5,*,0,3,0], linestyle=1, color='violet'


restore, 'SN2011erp_redbolmags161z.sav'
cgoplot, erpredmagsz[4,*,0,0,0]-erpredmagsz[5,*,0,0,0], erpredmagsz[1,*,0,0,0]-erpredmagsz[5,*,0,0,0], linestyle=1, color='red'
cgoplot, erpredmagsz[4,*,0,1,0]-erpredmagsz[5,*,0,1,0], erpredmagsz[1,*,0,1,0]-erpredmagsz[5,*,0,1,0], linestyle=1, color='red'
cgoplot, erpredmagsz[4,*,0,2,0]-erpredmagsz[5,*,0,2,0], erpredmagsz[1,*,0,2,0]-erpredmagsz[5,*,0,2,0], linestyle=1, color='red'
cgoplot, erpredmagsz[4,*,0,3,0]-erpredmagsz[5,*,0,3,0], erpredmagsz[1,*,0,3,0]-erpredmagsz[5,*,0,3,0], linestyle=1, color='red'



restore, 'SN2017erp_170702_uvopt_obswave_obsflux.dat_redbolmags161z.sav'
cgoplot, redzmags_array[4,*,0,0,0]-redzmags_array[5,*,0,0,0], redzmags_array[1,*,0,0,0]-redzmags_array[5,*,0,0,0], linestyle=1, color='red'
cgoplot, redzmags_array[4,*,0,1,0]-redzmags_array[5,*,0,1,0], redzmags_array[1,*,0,1,0]-redzmags_array[5,*,0,1,0], linestyle=1, color='red'
cgoplot, redzmags_array[4,*,0,2,0]-redzmags_array[5,*,0,2,0], redzmags_array[1,*,0,2,0]-redzmags_array[5,*,0,2,0], linestyle=1, color='red'
cgoplot, redzmags_array[4,*,0,3,0]-redzmags_array[5,*,0,3,0], redzmags_array[1,*,0,3,0]-redzmags_array[5,*,0,3,0], linestyle=1, color='red'


;;;;;;;;;;;;;;;;;
oploterror, host.bpeakappmag_array[*,4]-host.bpeakappmag_array[*,5], host.bpeakappmag_array[*,1]-host.bpeakappmag_array[*,5], sqrt(host.bpeakappmagerr_array[*,4]^2.0+host.bpeakappmagerr_array[*,5]^2.0), sqrt(host.bpeakappmagerr_array[*,1]^2.0+host.bpeakappmagerr_array[*,5]^2.0), symsize=0.3, psym=15


xyouts, host.bpeakappmag_array[where(host.snname_array eq 'SN2014J'),4]-host.bpeakappmag_array[where(host.snname_array eq 'SN2014J'),5]-0.05, host.bpeakappmag_array[where(host.snname_array eq 'SN2014J'),1]-host.bpeakappmag_array[where(host.snname_array eq 'SN2014J'),5] +0.2, 'SN2014J', charsize=0.5


;xyouts, host.bpeakappmag_array[where(host.snname_array eq 'SN2013gs'),4]-host.bpeakappmag_array[where(host.snname_array eq 'SN2013gs'),5]+0.05, host.bpeakappmag_array[where(host.snname_array eq 'SN2013gs'),1]-host.bpeakappmag_array[where(host.snname_array eq 'SN2013gs'),5] - 0.15, 'SN2013gs', charsize=0.5

;cgoplot, 15.54-15.37, 17.08-15.37, psym=4, color='red'
;xyouts, 15.54-15.37+0.25, 17.08-15.37, 'SN2013gs Zhang et al.', charsize=0.5

;xyouts, host.bpeakappmag_array[where(host.snname_array eq 'SN2006X'),4]-host.bpeakappmag_array[where(host.snname_array eq 'SN2006X'),5]-0.15, host.bpeakappmag_array[where(host.snname_array eq 'SN2006X'),1]-host.bpeakappmag_array[where(host.snname_array eq 'SN2006X'),5] +0.28, 'SN2006X', charsize=0.5



device, /close
SET_PLOT, 'X'
spawn, 'open bpeak_m2vbv_redvector.eps'


order=sort(host.bpeakappmag_array[*,5]-host.bpeakappmag_array[*,2])
; for i=0,30 do print, host.snname_array[order[i]], host.bpeakappmag_array[order[i],4]-host.bpeakappmag_array[order[i],5], host.bpeakappmag_array[order[i],2]-host.bpeakappmag_array[order[i],5]



print, 'final stop'
stop

end
