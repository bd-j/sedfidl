;+
; NAME:
;   PHOT_FIT_COMPRESS
;
; VERSION:
;   2.0 (Mar, 2010)
;
; PURPOSE:
;   given observed magnitudes and errors, and an array of model
;   magnitudes, calculate the chi^2 of each model, optimally
;   normalized, build the cumulative distribution functions for input
;   parameters, determine percentiles of the CDF, and report the best
;   fit parameters.
;
;
; REFERENCE:
;   Smiley, G. 2008 ApJ
;
; CATEGORY:
;   Spectral fitting
;
; CALLING SEQUENCE:
;   mass_out=phot_fit_compress(obs_mag,obs_magerr,model_mag,pars=,$
;                              parnames=, parnorm= [,nmag=,magptr=, bandmask=])
;
; INPUTS:
;
;    OBS_MAG    - an [n_obj x n_filt] vector of magnitudes
;    OBS_MAGERR - an [n_obj x n_filt] vector of magnitude errors
;    MODEL_MAG  - an [n_filt x n_model] array of model magnitudes
;    PARS       - an [n_model x n_pars] array of parameters associated
;                 with each model, for which CDFs will be built
;    PARNAMES   - an [n_par] string array of parameter names
;    PARNORM    - an [n_par] vector indicating whether each model
;                 parameter is to be multplied by the optimal
;                 normalization or not.  1 for yes, 0 for no
;
; OPTIONAL INPUTS:
;    NMAG      - the number of magnitude elements to include in the
;                output structures.  included for the case where the
;                SED being fitted has fewer observed wavelengths than
;                the general sample.
;    MAGPTR    - an [nobj x n_filt] pointer array indicating where in
;                the output structure magnitude elements the input
;                magnitudes should be saved.  again, for cases of
;                missing data. 
;    BANDMASK  - an [n_filt] mask indicating which magnitudes to fit.
;                defaults to fltarr(n_filt)+1
;
;
; KEYWORD PARAMETERS:
;
;
; OUTPUT:
;    Returns a structure containing all information about the fit:
;      chi2
;      par_cdf
;      par_grid
;
; OPTIONAL OUPUT:     
;     MINI_MASS - a structure: a summary of the CDF for each parameter,
;                 (marginalized over all others) 
;                 containing the median and percentiles
;                 as well as the bestfit value and the resdiuals, and
;                 input magnitudes (and errors thereon)
;
; COMMENTS:
;    Describe useful info
;
; REVISION HISTORY:
;    Mar 2010 - written, B. Johnson
;
;--------------------------------------------------------------------
FUNCTION phot_fit_compress,obs_mag,obs_magerr,model_mag,$
                  galid=galid,mini_mass=mini_mass,$
                  pars=pars,parnames=parnames,parnorm=parnorm,$
                  nmag=nmag,magptr=magptr,bandmask=bandmask


;window,1
;stop

sz=size(obs_mag)

if sz[0] GT 1 then begin
   nlvl=sz[1] & nfilt=sz[2] 
endif  else begin
  message,'data must be in 2-D arrays'
endelse
if n_elements(bandmask) EQ 0 then bandmask=fltarr(nfilt)+1

;allow for missing data
if n_elements(magptr) GT 0 then begin
   mptr=magptr
   nfilt=nmag
endif else mptr=findgen(nfilt)

psize=size(pars)
if psize[0] GT 1 then npars=psize[2] else npars=1
ngrid=psize[1]
opars=pars
pars=transpose(pars)
if keyword_set(parnorm) EQ 0 then parnorm=fltarr(npars) ;default to unnormalized parameters
npar=npars
mags=model_mag

if keyword_set(galid) then goo=galid else goo=findgen(nlvl)

for ilvl=0,nlvl-1 do begin
   ;get chi square
   chi2=find_xi_photo(mags,reform(obs_mag[ilvl,*]),reform(obs_magerr[ilvl,*]),norm=norm,mask=bandmask)
   bchi=min(chi2,imin)
   ;restrict to only those models that will have finite likelihood
   good=where(chi2-bchi LT 1416d,ngood) 
   if keyword_set(verbose) then $
      print,bchi,ngood


;OUTPUT STRUCTURES
   ;big structure
   mass_out={$
            gal_id:0.,corr_mags:fltarr(nfilt),magerrs:fltarr(nfilt),$
            ;mass_median:0.,mass_bestfit:0.,$
            delta_mag_best:fltarr(nfilt),$
            model_mag_best:fltarr(nfilt),$
            fit_mask:bandmask,$
            chi2:fltarr(ngrid) $
            }
   start=n_TAGS(mass_out)
   ns=strcompress(string(long(ngrid)),/remove_all)
   ngs=strcompress(string(long(ngood)),/remove_all)

if ngs EQ 0 then stop

   if npars GT 0 then $
      mass_out=struct_addtags(mass_out,[parnames+'_grid',parnames+'_cdf'],$
                              [strarr(npar)+'fltarr('+ns+')',strarr(npar)+'dblarr(2,'+ngs+')']) 
   ;small structure
   mini_mass={$
             gal_id:0.,corr_mags:fltarr(nfilt),magerrs:fltarr(nfilt),$
             delta_mag_best:fltarr(nfilt),model_mag_best:fltarr(nfilt),fit_mask:bandmask,$
             chibest:0.,bestfit_index:0l}
   ministart=n_TAGS(mini_mass)
   if npars GT 0 then $
      mini_mass=struct_addtags(mini_mass,[parnames+'_median',parnames+'_bestfit',parnames+'_2p5',$
                                          parnames+'_97p5'],strarr(npar*4)+'0.')
   mass_out[ilvl].chi2=chi2 
   mini_mass[ilvl].chibest=bchi
   mass_out[ilvl].delta_mag_best[mptr]=reform(obs_mag[ilvl,*])-(mags[*,imin]-2.5*alog10(norm[imin]))
   mass_out[ilvl].model_mag_best[mptr]=(mags[*,imin]-2.5*alog10(norm[imin]))

;note that final cdfs will have to multiplied by exp(-bchi/2) to get
;the real value out
  mini_mass[ilvl].bestfit_index=imin

;GENERIC PARAMETERS CDF
  for ipar=0, npar-1 do begin
     nn=(parnorm[ipar] EQ 0)+(parnorm[ipar] GT 0)*norm 
     mass_out[ilvl].(start+ipar)=pars[ipar,*]*nn  
     mass_out[ilvl].(start+npar+ipar)=transpose(build_cdf($
                      mass_out[ilvl].chi2[good]-bchi,mass_out[ilvl].(start+ipar)[good],ord=ord))
     mini_mass[ilvl].(ministart+npar+ipar)=mass_out[ilvl].(start+ipar)[imin]
  endfor

;GENERIC PARAMETERS PERCENTILES
  for ipar=0,npar-1 do begin
     j=nearest(mass_out[ilvl].(start+npar+ipar)[1,*]/max(mass_out[ilvl].(start+npar+ipar)[1,*]),$
               [0.0250,0.50,0.9750])
     if min(j) GE 0 then tmp='' else message, 'CDF percentiles not found for '+ts(galid[ilvl])
     mini_mass[ilvl].(ministart+ipar)=mass_out[ilvl].(start+npar+ipar)[0,j[1]]
     mini_mass[ilvl].(ministart+npar*3+ipar)=mass_out[ilvl].(start+npar+ipar)[0,j[2]]
     mini_mass[ilvl].(ministart+npar*2+ipar)=mass_out[ilvl].(start+npar+ipar)[0,j[0]]
  endfor
endfor

;;--------
;fill summary structure
;;---------
;mini_mass.photo_grid_name=gridname
mini_mass.gal_id=goo
mini_mass.corr_mags[mptr]=transpose(obs_mag)
mini_mass.magerrs[mptr]=transpose(obs_magerr)
mini_mass.delta_mag_best=mass_out.delta_mag_best
mini_mass.model_mag_best=mass_out.model_mag_best

;;------
;;fill cdf structure
;;-------
;if keyword_set(gridname)then mass_out.photo_grid_name=gridname
mass_out.gal_id=goo
mass_out.corr_mags[mptr]=transpose(obs_mag)
mass_out.magerrs[mptr]=transpose(obs_magerr)

pars=opars
return,mass_out

end
