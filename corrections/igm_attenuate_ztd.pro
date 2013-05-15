;+
; NAME:
;   IGM_ATTENUATION_Ztd
;
; VERSION:
;   1.0 (Apr, 2010)
;
; PURPOSE:
;   take a restframe wavelength vector and a set of redshifts and
;   return an array giving the igm attenuation as a function of
;   wavelength for each redshift
;
; REFERENCE:
;   Bruzual & Charlot, 2003, MNRAS
;   Meiksen 200?
;
; CALLING SEQUENCE:
;   
;
;
; CATEGORY:
;   Population Synthesis
;
; INPUTS:
;   wave - restframe wavelngths of the spectrum of each object
;   z    - N_obj vector of redshifts
;
; OPTIONAL INPUTS:
;   filterlist - list of kcorrect filters for which to compute

FUNCTION igm_attenuate_ztd,z,wave=wave,vers_sps=vers_sps
  nobj=n_elements(z)
  nzt=12
  if vers_sps EQ 'bc03' then $
     igm_trans=reform((mrdfits('$SPECFIT_DIR/data/meiksin_igmtrans_rest_bc03.fits.gz'))[1:*,*])
  if vers_sps EQ 'cb07' then $ 
     igm_trans=reform((mrdfits('$SPECFIT_DIR/data/meiksin_igmtrans_rest_cb07.fits.gz'))[1:*,*])
  ztigm=findgen(nzt)*0.5+1.5
  nw=(size(igm_trans,/dims))[1]

  ddzz=((1-ztigm)#(fltarr(nobj)+1)+(fltarr(nzt)+1.)#band_shift) < 0
  s1=total(ddzz EQ 0,1)-1
  w2=(band_shift-ztigm[s1])/0.5 ;linear interpolation weights
  w1=1-w2
  att=reform((w1#(fltarr(nw)+1.))*igm_trans[s1>0,*]+(w2#(fltarr(nw)+1.))*igm_trans[(s1+1)<11,*])
  
  return,att
  
end
