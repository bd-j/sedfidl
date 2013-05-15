FUNCTION hii_maggies, grid, filters, attenuate = attenuate, abs_neb = abs_neb, tau_neb = tau_neb, lyc_esc = lyc_esc

;;----------------------
;;Model nebular flux contributions to the continuum
;; takes as input a grid structure, which contains information on the
;; number of lyCcontinuum photons and the dust properties of each
;; model, and a filter list.  If /nodust is set then the calculation
;; is highly simplified and speeded up and attenuation of the nebular 
;; continuum is not considered.  Returns an nfilt x nmod array of
;; nebular fluxes in units of AB maggies.  optionally returns the
;; luminosity of the 
;;----------------------
nfilt = n_elements(filters)
nmod = n_elements(grid)
if n_elements(lyc_esc) EQ 0 then lyc_esc = fltarr(nmod)+0.
if n_elements(lyc_abs) EQ 0 then lyc_abs = fltarr(nmod)+0.
if n_elements(tau_neb)  EQ 0 then tau_neb=grid.tau_bc
if N_elements(neb_dust) GT 1 then dust_pars_neb=float(neb_dust[1:*])

hii_maggies=fltarr(nfilt,nmod)
wave=findgen(3E3)*18+912.
kwave=k_lambda_to_edges(wave)
HH=4.*!PI*(10.0d*!pc2cm)^2.
;;get nebular spectrum, cgs at a distance of 10pc
hii_flux=(hii_spec(wave=wave,logNly=lognly_cloudy,cloudyfile=cloudyfile))/(HH)

if keyword_set(attenuate) EQ 0 then begin
   norm = 10.^(grid.log_nly-lognly_cloudy)*(1.-lyc_esc)*(1.-lyc_abs)
   hii_maggies = reform(k_project_filters(kwave,reform(hii_flux),filterlist=filters))
   hii_maggies = hii_maggies # norm

endif else begin
  ;;split it into blocks to
  ;;avoid ememory constraints.  this could maybe be sped up through a rewrite
   memlimit=1E4
   nblock=ceil(nmod/memlimit)
   abs_neb=fltarr(nmod)
   for iblock=0,nblock-1 do begin
      if keyword_set(verbose) then $
         print,'emission line correction: block '+ts(iblock+1)+' of '+ts(nblock)
      lo=iblock*memlimit & hi=((iblock+1)*memlimit-1) < (nmod-1)
      ntb=hi-lo+1
      ;;normalize to number of Lyc photons in model, with corrections for lost photons
      norm_hii_flux=hii_flux#(10.^(grid[lo:hi].log_nly-lognly_cloudy)*(1.-lyc_esc[lo:hi])*(1.-lyc_abs[lo:hi])) 
      ;;attenuate nebular spectrum by bc dust
      ext_neb=exp(0.-call_function(neb_dust[0],wave,dust_pars_neb)#(tau_neb[lo:hi]))
      ;;project nebular spectra onto filters
   
      hii_maggies[*,lo:hi]=transpose(reform(k_project_filters(kwave,reform(norm_hii_flux*ext_neb),filterlist=filters)))
      ;;roughly estimate out how much nebular flux is absorbed by dust
      ni=where(wave GT 912,nwni)
      ww=k_lambda_to_edges(wave[ni])
      dwni=ww[1:*]-ww[0:n_elements(ww)-2]
      abs_neb[lo:hi]=total((dwni#(fltarr(nmod)+1))* hii_flux[ni,*]*(1-ext_neb[ni,*]),1)
   endfor
endelse

;stop
;;add to stellar magnitudes
;maggies=10.^(0.-mags/2.5)+hii_maggies
;mags=0.-2.5*alog10(maggies)

return, hii_maggies

end

