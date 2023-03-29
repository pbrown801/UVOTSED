
;;;;;;  this takes an input spectrum (ideally unreddened and at low redshift or with rest-frame wavelengths)
;;;;;;  and creates an array of its magnitudes at different redshifts and reddening values


pro makeredzarray, spectrum, spectra_labels=spectra_labels


nspectra=n_elements(spectrum)  ;; maybe later this can have multiple spectra
nlaws=4
nebv=81 ;161
nz=18; 
nfilters=6 ; just doing uvot for now

redzmags_array    =fltarr(nfilters,nebv,nspectra,nlaws,nz)
redzspectra_array =fltarr(641,     nebv,nspectra,nlaws)    ; spectra for different z's can be made by changing the lambda array

ebv=fltarr(nebv)
zarray=fltarr(nz)
for zindex=0,nz-1 do zarray[zindex]=zindex*0.002
for e=1,nebv-1    do ebv[e]=e*0.02

filename='spectrum_array_redbolmags161z.sav'
s=size(spectrum)
if (s[1] eq 7) and n_elements(spectrum) eq 1 then filename=spectrum+'_redbolmags161z.sav'

for n=0,nspectra-1 do begin


	; check to see if it is a string (ie a filename) or an array
	s=size(spectrum)
	if (s[1] eq 7) then begin
		readcol,spectrum[n], sp_wave,sp_flux,/silent
	endif else begin
		sp_wave=spectrum[0,*]
		sp_flux=spectrum[1,*]
	endelse

	sp_flux=sp_flux(sort(sp_wave))
	sp_wave=sp_wave(sort(sp_wave))

	sp_wave=sp_wave[where(finite(sp_flux) eq 1)]
	sp_flux=sp_flux[where(finite(sp_flux) eq 1)]
	;;;;;


	print, 'working on '+spectrum[n]

	pjb_uvotspec_all, spectrum[n], mag_array=mag_array, counts_array=counts_array, lambda=lambda, filter_array=filter_array, countspec_array=countspec_array,  spectrum_uvotres=spectrum_uvotres, fluxdensityfactors=fluxdensityfactors, flatflux=flatflux,  intflux=intflux, muvflux=muvflux, nuvflux=nuvflux, optflux=optflux, zeropoints=zeropoints, effwavelength=effwavelength, sp_wave=sp_wave, sp_flux=sp_flux



	for l=0,nlaws-1 do redzspectra_array[*,0,n,l]=spectrum_uvotres

	for e=0,nebv-1  do begin

		print, 'running models for E(B-V) = ', ebv[e]

		model=0
		rv=3.1
		lambda_av_mw31=sne_mw_reddening(lambda,Ebv[e], rv=rv)
		redzspectra_array[*,e,n,model]=redzspectra_array[*,0,n,model]*10^(-lambda_av_mw31/2.5)

		for zindex=0,nz-1 do begin
			z=zarray[zindex]
			pjb_uvotspec_all, [transpose(lambda*(1.0+z)), transpose(redzspectra_array[*,e,n,model])], mag_array=mag_array
	   		redzmags_array[*,e,n,model,zindex]=mag_array[0:5]
		endfor

		model=1
		rv=1.7
		lambda_av_mw17=sne_mw_reddening(lambda,Ebv[e], rv=1.7)
		redzspectra_array[*,e,n,model]=redzspectra_array[*,0,n,model]*10^(-lambda_av_mw17/2.5)

		for zindex=0,nz-1 do begin
			z=zarray[zindex]
			pjb_uvotspec_all, [transpose(lambda/(1.0+z)), transpose(redzspectra_array[*,e,n,model])], mag_array=mag_array
		   	redzmags_array[*,e,n,model,zindex]=mag_array[0:5]
		endfor

		model=2
		lambda_av_smc=sne_smc_reddening(lambda,Ebv[e])
		redzspectra_array[*,e,n,model]=redzspectra_array[*,0,n,model]*10^(-lambda_av_smc/2.5)


		for zindex=0,nz-1 do begin
			z=zarray[zindex]
			pjb_uvotspec_all, [transpose(lambda*(1.0+z)), transpose(redzspectra_array[*,e,n,model])], mag_array=mag_array
		   	redzmags_array[*,e,n,model,zindex]=mag_array[0:5]
		endfor

		model=3
		lambda_av_cslmc=sne_goobarlmc_reddening(lambda,Ebv[e])
		redzspectra_array[*,e,n,model]=redzspectra_array[*,0,n,model]*10^(-lambda_av_cslmc/2.5)

		for zindex=0,nz-1 do begin
			z=zarray[zindex]
			pjb_uvotspec_all, [transpose(lambda*(1.0+z)), transpose(redzspectra_array[*,e,n,model])], mag_array=mag_array
		   	redzmags_array[*,e,n,model,zindex]=mag_array[0:5]
		endfor

	endfor
endfor

;;;;;;;;;;;;;;;;;;;;;
save, filename=filename, redzmags_array, redzspectra_array, spectrum, spectra_labels, nlaws, ebv, zarray, lambda



stop
end