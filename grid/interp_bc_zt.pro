;+
; NAME:
;   INTERP_BC_MET
;
; VERSION:
;   1.0 (Oct, 2010)
;
; PURPOSE:
;   take two BC03 structures as produced by k_im_read_bc03_v3 (with
;   identical dust and SFH propertires)
;   and interpolate them to given ages and intermediate
;   metallicities.  Return the spectra at these values, and ancillary info.
;
; REFERENCE:
;
; CATEGORY:
;   SED fitting
;
; CALLING SEQUENCE:
;   bc=routine(z,t,bc1,bc2,z1,z2,[BC1EX=BC1EX,BC2EX=BC2EX,BCEX=BCEX,/LINEAR])
;
; INPUTS:
;    z - the desired metallicity of the population synthesis model
;    t - the desired age of the model
;    bc1 - 1-element structure containing a BC03 model with Z less
;          than the desired Z
;    bc2 - 1-element structure containing a BC03 model with Z greater 
;          than the desired Z
;    z1 - metallicity of the the model contained in bc1
;    z2 - metallicity of the model contained in bc2


; OPTIONAL INPUTS:
;    bc1ex - the 'extras' structure of the 1st model
;    bc2ex - the 'extras' structure of the 2nd model
;
; KEYWORD PARAMETERS:
;    /LINEAR - make the interpolation linear in metallicity rather
;              than logarithmic
;
; OUTPUT:
;    A structure containing the BC03 model with fluxes interpolated
;    (logarithmically) in metallicity from two models which bracket
;    the metallicity
;    BCEX - an 'extra' sructure with stellar mass and number of lyman
;           continuum photons interpolated logarithmicall in Z from
;           the two bracketing models
; COMMENTS:
;    Describe useful info
;
; REVISION HISTORY:
;    Aug 2009 - written, B. Johnson
;    Sep 2010 - adapted from interp_bc_met to include age
;               interpolation and to work on vectors
;
;--------------------------------------------------------------------

FUNCTION interp_bc_Zt,ztarg,ttarg,bc1,bc2,z1,z2,bc1ex=bc1ex,bc2ex=bc2ex,wave=wave


ages=bc1.age
nage=n_elements(ages)
nt=n_elements(ttarg)

SS1=bc1.flux
SS2=bc2.flux
nw=n_elements(bc1.wave)
wave=bc1.wave

;;--------------------
;; INTERPOLATION WEIGHTS
;;--------------------
   ;;age weights are from interpolation in log t
    tt1=ages#(fltarr(nt)+1)         ;host timeline
    tt2=(fltarr(nage)+1)#ttarg         ;burst timeline
    dt=tt1-tt2 
    dd=dt>0
    s2=total(dd EQ 0,1)-1 
    g1=where(s2 GE 0 and s2 LT nage-1,numg)
    if numg ne nt then print,'interp_bc_zt: warning, ages out of bounds'
    ll=ttarg[g1]
    kk=ages[s2[g1]]
    mm=ages[s2[g1]+1]
  ;;weights a and b are from interpolation in log(t) (following the BC03
  ;     add_bursts program)
    b=alog10(ll/kk)/alog10(mm/kk)
  ;very young bursts
    vy=where(s2[g1] EQ 0,nvy)
    if nvy GT 0 then b[vy]=1  ;use the first non-zero isochrone without weighting

    a=1-b
    
    onedindex=dindgen(nt)*nage+s2 ;into an nage by nt array
    UU=fltarr(nage,nt)
    UU[onedindex]=a
    UU[onedindex+1]=b

  ;;Z weights a and b are from interpolation in log(Z) 
   b=alog10(ztarg/z1)/alog10(z2/z1)
   bad=where(finite(b) EQ 0,nbad)
   if nbad GT 0 then b[bad]=0
   if keyword_set(linear) then b=(ztarg-z1)/(z2-z1)
   a=1-b
   AA=(fltarr(nw)+1)#a
   BB=(fltarr(nw)+1)#b



  ;;Do the interpolation to get flux (and mass and whatever)
   FF1=SS1#UU
   FF2=SS2#UU
   FF=AA*FF1+BB*FF2

  ;; return a structure
   bc={flux:fltarr(nw),Z:0.,t:0.,m_:0.,nly:0.,sfr_yr:0.,mgalaxy:0.,mbol:0.}
   bc=replicate(bc,nt)
   bc.flux=FF
   bc.t=ttarg
   bc.z=ztarg

  ;; add extra stuff
   bc.m_=reform(reform(bc1ex.m_,[1,220])#UU[1:*,*])*a+$
         reform(reform(bc2ex.m_,[1,220])#UU[1:*,*])*b
   bc.mgalaxy=reform(reform(bc1ex.mgalaxy,[1,220])#UU[1:*,*])*a+$
              reform(reform(bc2ex.mgalaxy,[1,220])#UU[1:*,*])*b
   n1=reform(10.^((bc1ex.nly)-30.),[1,220])
   n2=reform(10.^((bc2ex.nly)-30.),[1,220])
   bc.nly=alog10(reform(n1#UU[1:*,*])*a+reform(n2#UU[1:*,*])*b)+30.
   bc.sfr_yr=reform(reform(bc1ex.sfr_yr,[1,220])#UU[1:*,*])*a+$
             reform(reform(bc2ex.sfr_yr,[1,220])#UU[1:*,*])*b
   bc.mbol=alog10(reform(reform((10^bc1ex.mbol),[1,220])#UU[1:*,*])*a+$
                  reform(reform((10.^bc2ex.mbol),[1,220])#UU[1:*,*])*b)

;stop

return,bc

end
