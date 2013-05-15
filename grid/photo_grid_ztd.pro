;+
; NAME:
;   PHOTO_GRID_Ztd
;
; VERSION:
;   1.0 (Apr, 2010)
;
; PURPOSE:
;   to create a grid of flux points for a given set of model ages and
;   parameters.  Basically a wrapper on quantities_ztd to isolate Z
;   ranges and make a nice structure of the output
;
; REFERENCE:
;   Bruzual & Charlot, 2003, MNRAS
;
; CALLING SEQUENCE:
;   
;
;
; CATEGORY:
;   Population Synthesis
;
; INPUTS:
;   met - [n_model] length vector of the metallicity of the population
;   tsf - [n_model] length vector of the timescale for SFR
;   tdiff - [n_model] length vector of the tau_v of the dust affecting
;           all stars 
;   tbc - [n_model] length vector of the tau_v of the dust
;         additionally affecting young (t<10^7 yrs) stars
;   age - [n_model] length vector of ages of model galaxy
;
; OPTIONAL INPUTS:
;   filterlist - list of kcorrect filters for which to compute
;                magnitudes
;   rundir - path to the directory containing the SPS ised files
;   attenuate_ionizing - fraction of ionizing flux assumed to go into
;                        the HII regions (i.e., covering factor)
;   vers_sps - sps models to use, default is 'cb07'
;   sfhtype - functional form of the sfh.  Choices are 'delay','exp',
;             or eventually, 'ssp' or 'step'
;   imftype - the imf to use.  'salp' or 'chab'
;   band_shift - [n_model] vector of redshifts
;   bc_dust - name of the function that returns tau_lambda/tau_5500
;             given lambda, extra parameters may be passed as extra
;             elements of the array.  This is for the dust that
;             affects young stars
;   bc_diff - name of the function that returns tau_lambda/tau_5500
;             given lambda, extra parameters may be passed as extra
;             elements of the array.  This is for diffuse dust that
;             affects all stars equally


FUNCTION photo_grid_Ztd,met,tsf,tdiff,tbc,age,$
                        filters=filters,band_shift=band_shift,$
                        vers_sps=vers_sps,sfhtype=sfhtype,imftype=imftype,$
                        bc_dust=bc_dust,diff_dust=diff_dust,$
                        attenuate_ionizing=attenuate_ionizing,$
                        rundir=rundir,keep_spectra=keep_spectra
;;----------
;;SETUP AND DEFAULTS
;;------------

memlimit=5E3

;SPS version
if keyword_set(vers_sps) EQ 0 then vers_sps='cb07'
if keyword_set(imftype) EQ 0 then imftype='chab'
if keyword_set(sfhtype) EQ 0 then sfhtype='exp'
sfhroot=vers_sps+imftype+'_'+sfhtype 
if keyword_set(rundir) EQ 0 then  rundir='$SPECFIT_LIB/'+vers_sps+'/'
if vers_sps EQ 'bc03' then nw=6900 else nw=6917

;SPS parameters
zbc=[0.02,0.2,0.4,1.0,2.5]
npz=n_elements(zbc)
nmet=n_elements(met)
nsf=n_elements(tsf)
nsed= nmet

;Default filters
gg=['galex_FUV','galex_NUV']+'.par'
ss='sdss_'+['u','g','r','i','z']+'0.par'
tt='twomass_'+['J','H','Ks']+'.par'
ii='spitzer_irac_'+['ch1','ch2','ch3']+'.par'
if n_elements(filters) EQ 0 then filters=[gg,ss,tt,ii]
nfilt=n_elements(filters)

;Redshifts
if n_elements(band_shift) EQ 0 then band_shift=fltarr(nsed)+0.0 $
   else if n_elements(band_shift) EQ 1 then band_shift=fltarr(nsed)+band_shift

;;--------------
;;CREATE STRUCTURE
;;------------

grid={z:0.,tau_sf:0.,tau:0.,tau_bc:0.,age:0.,band_shift:0.,$   ;input
      dust_law:' ',dust_law_bc:' ',sfhtype:' ',$               ;input
      mag:fltarr(nfilt),mag_nodust:fltarr(nfilt),$
      mbol:0.,mstar:0.,mgal:0.,sfr:0.,$
      log_ltrans:0.,log_ldust:0.,log_ldust_young:0.,log_ldust_bc:0.,$
      log_lion:0.,log_nly:0.,log_nly_dust:0.,modelname:' '} 
if keyword_set(keep_spectra) then $
   grid=struct_addtags(grid,['spectra'],['fltarr('+string(nw,format='(I4.4)')+')'])
grid=replicate(grid,nsed)

;;--------------
;;LOOP OVER MODEL PARAMETER COMBINATIONS
;;---------------

;;drill down into all models with the same tsf, and in a range
;;of met
usf=tsf[uniq(tsf,sort(tsf))]
npsf=n_elements(usf)
for isf=0,npsf-1 do begin
   this_sf=where(tsf EQ usf[isf])
   fsubs=this_sf
   zz=met[fsubs]
   for iz=0,npz-2 do begin
;      print,'z_index=',iz
      this_zr=where(zz GE zbc[iz] and zz LT zbc[iz+1],ntz)
      if ntz EQ 0 then continue ;no z's in this range
      allsubs=fsubs[this_zr]
      ntz=n_elements(allsubs)
 
   ;;break into pieces to avoid memory constraints
      ;print,ntz/memlimit
      nblock=ceil(ntz/memlimit)
      for iblock=0,nblock-1 do begin
         ;print,ntz,' models ',(iblock*memlimit),'->',((iblock+1)*memlimit)<(ntz-1)
         ;print,iblock+1,'of',nblock
         subs=allsubs[(iblock*memlimit):((iblock+1)*memlimit)<(ntz-1)]

   ;;read SPS appropriate for the model
   ;;pars, including interpolation in z
   ;;and t
         mags=quantities_ztd(met[subs],age[subs],tbc[subs],tdiff[subs],$
                             sfh=sfhroot+strcompress(string(fix(usf[isf]*1000),format='(I6.6)'),/remove_all),$
                             bc_dust=bc_dust,diff_dust=diff_dust,$
                             band_shift=band_shift,filters=filters,$
                             attenuate_ionizing=attenuate_ionizing,$
                             vers_sps=vers_sps,/zdir,$
                             log_lums=log_lums,bcextras=bcextras,$  
                             mags_nodust=mags_nodust,spectra=spectra)              

;Calculate quantities ---------
            
         grid[subs].mag=transpose(mags)
         grid[subs].mag_nodust=transpose(mags_nodust)
         grid[subs].mstar=bcextras[*,0]
         grid[subs].mgal=bcextras[*,1]
         grid[subs].sfr=bcextras[*,2]
         grid[subs].mbol=bcextras[*,5]
         grid[subs].log_nly_dust=bcextras[*,4]
         if keyword_set(keep_spectra) then grid[subs].spectra=spectra
      
   ;get dust and LyC photons
         grid[subs].log_nly=bcextras[*,3]
         grid[subs].log_ltrans=log_lums[*,0]
         grid[subs].log_lion=log_lums[*,1]
         grid[subs].log_ldust=log_lums[*,2]
         grid[subs].log_ldust_young=log_lums[*,3]
         grid[subs].log_ldust_bc=log_lums[*,4]
         ;grid[subs].log_ldust_kappa
      endfor
   endfor  
endfor



;;----------
;;UNIT CONVERSIONS
;;----------

;convert to absolute magnitudes (given that bc spectra are in L_sol/AA
;and k_correct expects erg/s/cm^2/AA)
K=(!lsun)/(4.*!PI)/(!pc2cm)^2./100. ;distance of 10pc for absolute magnitudes
grid.mag=grid.mag-2.5*alog10(K)
grid.mag_nodust=grid.mag_nodust-2.5*alog10(K)
if keyword_set(keep_spectra) then grid.spectra=grid.spectra*K

grid.log_ldust+=alog10(K)     ;convert L_ir to F_ir in erg/s/cm^2 at 10pc
grid.log_lion+=alog10(K)    ;convert L_ion to F_ion in erg/s/cm^2 at 10pc
grid.log_ltrans+=alog10(K)
grid.log_ldust_young+=alog10(K)
grid.log_ldust_bc+=alog10(K)
;endif

;grid.modelname=sednames
grid.z=met
grid.tau_sf=tsf
grid.tau=tdiff
grid.tau_bc=tbc
grid.age=age
grid.band_shift=band_shift
grid.sfhtype=sfhtype
grid.dust_law_bc=string(bc_dust+'_', $
                        format ='('+string(n_elements(bc_dust+'_'))+'A)')

grid.dust_law=string(diff_dust+'_', $
                        format ='('+string(n_elements(diff_dust+'_'))+'A)')

return,grid

end
