;+
; NAME:
;  INIT_MODELS_Zt
;
; VERSION:
;   1.0 (Apr, 2010)
;
; PURPOSE:
;   To take a hosts(+bursts) structure, perform emission line
;   corrections, and return a n_mod X n_filt array of corrected model
;   SEDs, optionally including emission line and IR flux points
;
; REFERENCE:
;   Bruzual & Charlot, 2003, MNRAS
;
; CALLING SEQUENCE:
;   
;
; CATEGORY:
;   Population Synthesis
;
; INPUTS:
;   grid - structure containing the output from add_bursts_Zt (or photo_grid_Zt)
;   
; OPTIONAL INPUTS:

FUNCTION INIT_MODELS_Zt,grid,filters=filters,$
                        grid_pars=grid_pars,sfr_avg=sfr_avg,$
                        bc_dust=bc_dust,diff_dust=diff_dust,neb_dust=neb_dust,$
                        lyc_esc=lyc_esc,lyc_abs=lyc_abs,tau_neb=tau_neb,$
                        emcorr=emcorr,redshifts=redshifts,verbose=verbose;,scale=scale

K=alog10(4*!PI)+2.0*alog10(3.086)+38
nfilt=n_elements(filters)
nmod=n_elements(grid)

mags=grid.mag

if N_elements(bc_dust) GT 1 then dust_pars_bc=float(bc_dust[1:*])
if N_elements(diff_dust) GT 1 then dust_pars_diff=float(diff_dust[1:*])
if N_elements(neb_dust) GT 1 then dust_pars_neb=float(neb_dust[1:*])

caseB=[11.87]                   ;Ha, conversion from Q to erg/s for case B, ~solar metallicity
ewave=[6563.]                   ;Ha

 ;don't absorb or allow to escape any lyman continuum photons.. 
if n_elements(lyc_esc)  EQ 0 then lyc_esc=fltarr(nmod)+0.
if n_elements(lyc_abs)  EQ 0 then lyc_abs=fltarr(nmod)+0.
fudge=1.0 ;parameter giving the ratio of the average energy per absorbed continuum photon to the average energy of a continuum photon.  In truth this is probably lyc_abs dependent, as it is determined by the absorption curve.  1.0 is grey (in some units?) 

 ;make dust the same as the young stars
if n_elements(neb_dust) EQ 0 then begin 
   neb_dust=bc_dust
   dust_pars_neb=dust_pars_bc
endif
if n_elements(tau_neb)  EQ 0 then tau_neb=grid.tau_bc


;;-------------
;;reconstruct SFH to get to the SFR averaged over different timescales
;;--------------
if keyword_set(sfr_avg) then $
   sfr_avg=reconstruct_sfh(grid,scale=1.0)

;;----------------------
;;Model nebular flux contributions to the continuum??
;;----------------------
if keyword_set(emcorr) then begin
  ;split it into blocks to
  ;avoid ememory constraints
   memlimit=1E4
   nblock=ceil(nmod/memlimit)
   hii_maggies=fltarr(nfilt,nmod)
   abs_neb=fltarr(nmod)
   wave=findgen(3E3)*18+912.
   HH=4.*!PI*(10.0d*!pc2cm)^2.
   ;get nebular spectrum, cgs at a distance of 10pc
   hii_flux=(hii_spec(wave=wave,logNly=lognly_cloudy,cloudyfile=cloudyfile))/(HH)

   for iblock=0,nblock-1 do begin
      if keyword_set(verbose) then $
         print,'emission line correction: block '+ts(iblock+1)+' of '+ts(nblock)
      lo=iblock*memlimit & hi=((iblock+1)*memlimit-1) < (nmod-1)
      ntb=hi-lo+1
     ;normalize to number of Lyc photons in model, with corrections for lost photons
      norm_hii_flux=hii_flux#(10.^(grid[lo:hi].log_nly-lognly_cloudy)*(1.-lyc_esc[lo:hi])*(1.-lyc_abs[lo:hi])) 
     ;attenuate nebular spectrum by bc dust
      ext_neb=exp(0.-call_function(neb_dust[0],wave,dust_pars_neb)#(tau_neb[lo:hi]))
     ;project nebular spectra onto filters
      kwave=k_lambda_to_edges(wave)
      hii_maggies[*,lo:hi]=transpose(reform(k_project_filters(kwave,reform(norm_hii_flux*ext_neb),filterlist=filters)))
      ;;roughly estimate out how much nebular flux is absorbed by dust
      ;;;if you dont care about this, then the code can be speeded up
      ;;;by projecting filters before normalizing.  takes a lot less memory. 
      ni=where(wave GT 912,nwni)
      ww=k_lambda_to_edges(wave[ni])
      dwni=ww[1:*]-ww[0:n_elements(ww)-2]
      abs_neb[lo:hi]=total((dwni#(fltarr(nmod)+1))* hii_flux[ni,*]*(1-ext_neb[ni,*]),1)
   endfor

   ;add to stellar magnitudes
   maggies=10.^(0.-mags/2.5)+hii_maggies
   mags=0.-2.5*alog10(maggies)
endif

;;----------------------
;;Determine model emission line flux(es)
;;----------------------
;if keyword_set(do_ha) then begin
 ;make sure observed quantities in corresponding
 ;units (i.e. -2.5*log(erg/s/cm^2) - distance modulus 
   nline=n_elements(ewave)
   mag_eline=0.-2.5*(grid.log_nly)+2.5*K-2.5*(alog10(1.-lyc_esc)+alog10(1-lyc_abs))
   mag_eline=0.-2.5*alog10((10.0^(0.-mag_eline/2.5))#(10.^(0-caseB)))
   a_line=1.086*(tau_neb#call_function(neb_dust[0],ewave,dust_pars_neb)+reform(grid.tau)#call_function(diff_dust[0],ewave,dust_pars_diff))
   mag_emline=mag_eline+a_line
   mags=[mags,reform(mag_emline,[nline,nmod])]
;endif

;;----------------------
;;Determine model IR flux
;;----------------------
;if keyword_set(do_irmag) then begin
   ;erg/s/cm^2/M_sun at a distance of 10pc. (conversion is done when building the grid) 
   mag_tir=0.-2.5*alog10(10.^grid.log_ldust+lyc_abs*fudge*10.^grid.log_lion)
   mags=[mags,reform(mag_tir,[1,nmod])]
   ;split into multiple components to match to IR colors?
   ;e.g. color=alog10(lyc_abs*10.^grid.log_lion/10.^grid.log_ldust)
;endif

;;----------------------
;;get magnitudes in the observed space, i.e. add the distance modulus
;;----------------------

if keyword_set(redshifts) then begin
   dm=5.0*alog10(lumdist(grid.band_shift,Omega_M=0.27,lambda0=0.73))+25.
   K=dm-2.5*alog10(1+grid.band_shift) 
  ;do i need the second term?  depends on the operation of k_project_filters 
  ;yes, you need the second term, but must *subtract* 2.5log(1+z) from
  ;the distance modulus. This ***has been checked***

   mags=mags+K ;need to remove the 2.5log(1+z) term for Ha and IR, but these data likely won't be present if at high-z
endif


;grid_pars=[grid.mass,grid.sfr,grid.tau,grid.tau_bc,grid.age,grid.met]  ;trivial


return, mags
end
