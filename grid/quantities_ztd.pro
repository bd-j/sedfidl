;+
; NAME:
;   QUANTITIES_ZTD
;
; VERSION:
;   1.0 (Apr, 2010)
;
; PURPOSE:
;   to create a grid of flux points for a given set of model ages and
;   in a small range of metallicities.  Also apply dust attenuation
;   and estimate IR fluxes. Only returns results for *1* named SFH
;
; REFERENCE:
;   Bruzual & Charlot, 2003, MNRAS
;
; CALLING SEQUENCE:
;
; CATEGORY:
;   Population Synthesis
;
; INPUTS:
;   met - [n_model] length vector of the metallicity of the population
;   tdiff - [n_model] length vector of the tau_v of the dust affecting
;           all stars 
;   tbc - [n_model] length vector of the tau_v of the dust
;         additionally affecting young (t<10^7 yrs) stars
;   age - [n_model] length vector of ages of model galaxy (yrs)
;
; (OPTIONAL) INPUTS:
;   sfh - string specifiying the beginning of the SFH file to use
;         (i.e. 'cb07_exp10' for the cb07 exponential decline with
;         1-Gyr e-folding time)
;   filters - list of kcorrect filters for which to compute
;                magnitudes
;   rundir - path to the directory contining the SPS ised files
;   attenuate_ionizing - fraction of ionizing flux assumed to go into
;                        the HII regions (i.e., covering factor)
;   vers_sps - if set, use bc03, otherwise use cb07
;   band_shift - [n_model] vector of redshifts
;   bc_dust - name of the function that returns tau_lambda/tau_5500
;             given lambda, extra parameters may be passed as extra
;             elements of the array.  This is for the dust that
;             affects young stars
;   diff_dust - name of the function that returns tau_lambda/tau_5500 
;               given lambda, extra parameters may be passed as extra
;               elements of the array.  This is for diffuse dust that
;               affects all stars equally
;
; OUTPUT:
;    mags - [Nfilter x n_model] array of model absoulute magnitudes
;
; OPTIONAL OUPUT:
;    bcextras  - [6 x n_model] array of 'extra' model
;                parameters. These are: Stellar mass, total stellar
;                mass formed, instantaneous SFR, number of lyman
;                continuum photons, number of lyman continuum photons
;                if the spectrum was attenuated by normal dust
;                shortward of 912AA, and absolute bolometric magnitude
;
;    log_lums  - [5 x n_model] array of calculated model
;                luminosities.  These are: transmitted nonionizing
;                luminosity, ionizing luminosity, luminosity absorbed
;                by dust, luminosity of young stars absorbed by dust,
;                luminosity absorbed by birth cloud dust.
;
;    mags_nodust - as for mags, but giving the intrinsic, dust free
;                  magnitudes.
;
;
; COMMENTS:
;   It probably creates too many arrays, which makes it even more of a
;   memory hog than it needs to be. 
;
;   Luminosties are calculated from integrals blueward of 10 micron,
;   because int_tabulated is imperfect.
;   However, the effect on the integrals of the spectrum redward of 10
;   microns should be minimal.
;
; HISTORY:
;   written, Nov. 2010, B. Johnson

FUNCTION quantities_ztd,met,age,tbc,tdiff,sfh=sfh,$
                        bc_dust=bc_dust,diff_dust=diff_dust,$
                        band_shift=band_shift,filters=filters,$
                        attenuate_ionizing=attenuate_ionizing,$
                        vers_sps=vers_sps,zdir=zdir,$
                        rundir=rundir,$
                        spectra=spectra,$
                        log_lums=log_lums,bcextras=bcextras,$ ;output
                        mags_nodust=mags_nodust               ;output

  ntb=n_elements(met)
  if keyword_set(band_shift) EQ 0 then band_shift=fltarr(ntb)
  if arg_present(bc_dust) EQ 0 then bc_dust=['power_ztd','1.3']
  if arg_present(diff_dust) EQ 0 then diff_dust=['power_ztd','0.7']  
 
;read SPS
  specfit_libpars,lib_mu,lib_tau_v,lib_alpha_bc,lib_alpha,vers_sps=vers_sps

  bc = bdj_read_sps_Zt(met, age, lib_mu,lib_tau_v, sfh, zdir=zdir,$
                       rundir=rundir,wave=wave,vers_sps=vers_sps)
  bc_nd = bdj_read_sps_Zt(met, age, lib_mu, 0., sfh, zdir=zdir,$
                          rundir=rundir,vers_sps=vers_sps) 
  ff=reform(bc.flux)
  ff_nd=reform(bc_nd.flux)



;library stuff and dust curve initialization
  lib_tau_bc=lib_tau_v*(1-lib_mu)*(wave/5500)^(0.-lib_alpha_bc)
  lib_tau_old=lib_tau_v*lib_mu*(wave/5500)^(0.-lib_alpha)
  lib_tau_young=lib_tau_bc+lib_tau_old
  factor=exp(0.-lib_tau_young)-exp(0.-lib_tau_old)
  if N_elements(bc_dust) GT 1 then dust_pars_bc=float(bc_dust[1:*])
  if N_elements(diff_dust) GT 1 then dust_pars_diff=float(diff_dust[1:*])
  tau_law_bc=call_function(bc_dust[0],wave,dust_pars_bc)            
  tau_law_diff=call_function(diff_dust[0],wave,dust_pars_diff)      

;stop

;Modify flux array -------
 ;determine young (t<10^7) and old spectra automagically from
 ;the dust free and CF00 attenuated library spectra.  I'm pretty
 ;sure this works even on the Z+t interpolated spectra, since the
 ;interpolation is linear

  lbc=ff-ff_nd*(exp(0.-lib_tau_old)#(fltarr(ntb)+1))
  fyoung=temporary(lbc)/(factor#(fltarr(ntb)+1)) ;flux of the young component
  fold=ff_nd-fyoung

 ;attenuate each by dust
  tau_bc=tau_law_bc#tbc
  tau_old=tau_law_diff#tdiff
  tau_young=tau_bc+tau_old
  ff=fyoung*exp(0.-tau_young)+fold*exp(0.-tau_old)


 ;ionizing flux that went into HII regions
  if keyword_set(attenuate_ionizing) then begin
     wvptr_ion=where(wave LE 912.)
     ff[wvptr_ion,*]=ff[wvptr_ion,*]*(1-attenuate_ionizing)
  endif

  ;Apply IGM transmission NEEDS TO BE FIXED/IMPLEMENTED
      ;SINCE BANDSHIFT MAY BE DIFFERENT FOR EACH ELEMENT OF FLUX ARRAY
      ;igm_curve=
  if max(band_shift) GT 1.5 then $
     igm_att=igm_attenuate_ztd(band_shift,vers_sps=vers_sps) else $
        igm_att=1.0

;outputs and derived quantities
 
   mags=reform(bc_vals_new(wave,ff*igm_att,filters,$
                           band_shift=band_shift))
   mags_nodust=reform(bc_vals_new(wave,ff_nd*igm_att,filters,$
                           band_shift=band_shift))

   wvptr=where(wave GT 912. and wave LT 1E5,nwn) 
   wvptr_ion=where(wave LE 912.,nwi)

   wion=k_lambda_to_edges(wave[wvptr_ion])
   wni=k_lambda_to_edges(wave[wvptr])
   dwion=wion[1:*]-wion[0:nwi]
   dwni=wni[1:*]-wni[0:nwn]


   ;do some really crude integrals
   l_intrinsic=total((dwni#(fltarr(ntb)+1))*ff_nd[wvptr,*],1)
   log_ltrans=alog10(total((dwni#(fltarr(ntb)+1))*ff[wvptr,*],1))
   log_lion=alog10(total(dwion#(fltarr(ntb)+1)*ff_nd[wvptr_ion,*],1))
   log_ldust=alog10(l_intrinsic-10.^log_ltrans)
   nodust=where(abs(1-l_intrinsic/10^log_ltrans) LT 1E-5,nnd) ;not enough dust to matter
   ;if total(finite(log_ldust)) NE ntb then stop
   log_ldust_young=alog10(total((dwni#(fltarr(ntb)+1))*fyoung[wvptr,*]*$
                                (1-exp(0.-tau_young[wvptr,*])),1 ))
   log_ldust_bc=alog10(total((dwni#(fltarr(ntb)+1))*fyoung[wvptr,*]*$
                             (1-exp(0.-tau_bc[wvptr,*])),1))
   ;log_ldust_kappa=alog10(total(((dwni*kappa[wvptr])#(fltarr(ntb)+1))*ff_nd[wvptr,*]))
   log_lums=[[log_ltrans],[log_lion],[log_ldust],[log_ldust_young],[log_ldust_bc]];,$
;             [alog10(l_intrinsic)],[log_ldust_kappa]]
   if nnd GT 0 then log_lums[nodust,2:4]=-99
   bcextras=[[bc.m_],[bc.mgalaxy],[bc.sfr_yr],[bc_nd.nly],[bc.nly],[bc.mbol]]
   
   spectra=ff
   return,mags

end
