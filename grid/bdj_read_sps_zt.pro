;+
; NAME:
;   INTERP_BC_MET
;
; VERSION:
;   0.0 (Oct, 2008)
;
; PURPOSE:
;   Given a range of metallicities and ages, a single set of dust
;   properties, and a string specifiying the the start of a bc03 .ised
;   filename (which encodes the SFH), read the relevant *.ised files
;   and call interp_bc_zt to interpolate them in metallicity and age.
;   Returns the interpolated bc structure
;
; REFERENCE:
;   Smiley, G. 2008 ApJ
;
; CATEGORY:
;   Spectral fitting
;
; CALLING SEQUENCE:
;   bc=bdj_read_sps(met,ages,mu,tauv,sfhroot,[rundir=rundir,wave=wave,endian=endian,/ZDIR, /VERS_BC03])
;
; INPUTS:
;    met - [Nmod] vector of desired metallicities.  IMPORTANT:
;          metallicities should all be bracketed by two adjacent BC03
;          precompiled values.  Units are Z/Z_sun
;    ages - [Nmod] vector of desired model ages, in yrs
;    mu - scalar. the mu parameter for this set of models.  will generally be 0.1
;    tauv - scalar. the 5500AA optical depth for this set of models,
;           will be 0.0 or 1.0 usually
;    sfhroot - string giving the starting root of the ised filename,
;              which should encode the sfh (e.g. 'cb07salp_exp0100000' is a 10 Gyr
;              declining exponential, cb07chab_delay000010 is 10Myr delayed SFH)
;    
;
; OPTIONAL INPUTS:
;    endian - 'big' or 'little'
;
; KEYWORD PARAMETERS:
;    rundir - the directory containing the ised files.  Defaults to
;             $SPECFIT_LIB/cb07 (or bc03)/runs/
;    zdir -  append metallicity subdirectory to the path to the ised files.
;    vers_sps - sps models to use, default is 'cb07'
;
; OUTPUT:
;    bc - an [Nmod] element structure containing the flux vectors and
;         associated information for the interpolated models
;
; OPTIONAL OUPUT:
;    wave - the output wavelength vector for the bc structure (the
;           normal bc wave vector in the bc structure gets obliterated
;           in the interpolation process).  This is matched in length
;           to the bc.flux vector.
;
; COMMENTS:
;    Describe useful info
;
; REVISION HISTORY:
;    Sep 2010 - written, B. Johnson
;
;--------------------------------------------------------------------

FUNCTION bdj_read_sps_Zt,met,ages,mu,tauv,sfhroot,rundir=rundir,$
                         endian=endian,zdir=zdir,$
                         wave=wave,vers_sps=vers_sps

if keyword_set(vers_sps) EQ 0 then vers_sps='cb07'
if keyword_set(rundir) EQ 0 then rundir='$SPECFIT_LIB/'+vers_sps+'/'
if keyword_set(endian) eq 0 then endian='little'

if vers_sps EQ 'cb07' or vers_sps EQ 'bc03' then $
   mgrid=[0.02,0.2,0.4,1.0,2.5] else $
      mgrid=[0.02,0.2,0.4,1.0,2.5]
nz=n_elements(mgrid)

;should do more checking here of metallicities
m=nearest(mgrid,met)
m1=mgrid[(m-(m GT 0)*(mgrid[m] GT met))]
m2=mgrid[m+(m LT (nz-1))*(mgrid[m] LT met)] 
m1=min(m1)
m2=max(m2)
if n_elements(m1) NE 1 then stop ;print,'bdj_read_sps_Zt: metallicities not right'
if n_elements(m2) NE 1 then stop ;print,'bdj_read_sps_Zt: metallicities not right'
m1=m1[0]
m2=m2[0]

;stop

name1=sfhroot+$
      '_z'+strcompress(string(m1*100,format='(I3.3)'),/remove_all)+$
      '_mu'+strcompress(string(mu*10,format='(I2.2)'),/remove_all)+$
      '_tv'+strcompress(string(tauv*10,format='(I2.2)'),/remove_all)+'.'
if keyword_set(zdir) then $
   name1='z'+strcompress(string(m1*100,format='(I3.3)'),/remove_all)+'/'+name1

name2=sfhroot+$
      '_z'+strcompress(string(m2*100,format='(I3.3)'),/remove_all)+$
      '_mu'+strcompress(string(mu*10,format='(I2.2)'),/remove_all)+$
      '_tv'+strcompress(string(tauv*10,format='(I2.2)'),/remove_all)+'.'
if keyword_set(zdir) then $
   name2='z'+strcompress(string(m2*100,format='(I3.3)'),/remove_all)+'/'+name2



if vers_sps EQ 'bc03' then begin $
   bc1=k_im_read_bc03_v3(isedpath=rundir,isedfile=name1+'ised',$
                         extraname=name1,bc03_extras=bcex1,endian=endian,/silent) 
   bc2=k_im_read_bc03_v3(isedpath=rundir,isedfile=name2+'ised',$
                         extraname=name2,bc03_extras=bcex2,endian=endian,/silent) 
   
endif else if vers_sps EQ 'cb07' then begin $
   bc1=k_im_read_cb07_v3(isedpath=rundir,isedfile=name1+'ised',$
                         extraname=name1,cb07_extras=bcex1,endian=endian,/silent) 
   bc2=k_im_read_cb07_v3(isedpath=rundir,isedfile=name2+'ised',$
                         extraname=name2,cb07_extras=bcex2,endian=endian,/silent)
endif else if vers_sps EQ 'fsps' then begin $
   print,'not implemented'
   bc1=bdj_read_fsps(rundir+name2,extras=bcex1)
   bc2=bdj_read_fsps(rundir+name2,extras=bcex1)
endif

bc=interp_bc_Zt(met,ages,bc1,bc2,m1,m2,bc1ex=bcex1,bc2ex=bcex2,wave=wave)
   
return,bc

end

