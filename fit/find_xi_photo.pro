;+
; NAME:
;   FIND_XI_PHOTO
;
; VERSION:
;   0.0 (Oct, 2008)
;
; PURPOSE:
;   find the chi2 given observed magnitudes, and errors and a grid of
;   model magnitudes
;
; REFERENCE:
;   Smiley, G. 2008 ApJ
;
; CATEGORY:
;   Spectral fitting
;
; CALLING SEQUENCE:
;   result=routine(arg1,arg2,arg3,[ARG4=,OUT1=,/SWITCH1, /SWITCH2])
;
; INPUTS:
;    arg1 - the first argument
;    arg2 - the second argument.  It might have a long description in
;       which case the indentation would be like this.  just like this
;       with so many spaces.
;
; OPTIONAL INPUTS:
;    arg3 - the optional input
;
; KEYWORD PARAMETERS:
;    ARG4 - the keyword that must have a value
;    SWITCH1 -  a switch that can be set
;
;
; OUTPUT:
;    Describe the returned result.
;
;    OUT1 - the other returned output
;
; COMMENTS:
;    Describe useful info
;
; REVISION HISTORY:
;    Mar 2008 - written, B. Johnson
;
;--------------------------------------------------------------------

FUNCTION find_xi_photo,mag_grid,mags,magerr,mask=mask,norm=norm

;stop ;for debugging

sz1=size(mag_grid)
sz2=size(mags)
nmod=sz1[2]
nfilt=sz1[1]
if keyword_set(mask) EQ 0 then mask=fltarr(nfilt,nmod)+1

mags=mags#(fltarr(nmod)+1)
magerr=magerr#(fltarr(nmod)+1)

fmod=10.^(0.-0.4*mag_grid) ;maggies
fobj=10.^(0.-0.4*mags)     ;maggies
ivar=magerr*fobj/1.086
ivar=1/ivar^2
norm=total( (fmod*fobj*ivar)[where(mask),*],1 )/total( (fmod^2*ivar)[where(mask),*],1 )  ;from Salim et al 2007

chi2=total( ((fobj-((fltarr(nfilt)+1)#norm)*fmod)^2*ivar)[where(mask),*],1)  ;from Salim et al 2007

return,chi2

end
