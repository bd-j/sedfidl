;+
; NAME:
;   BURST_GRID_Ztd
;
; VERSION:
;   1.0 (Apr, 2010)
;
; PURPOSE:
;   to generate fluxes from a set of bursts distributed randomly in
;   the past of a given set of models with given metallicities
;
; REFERENCE:
;   Bruzual & Charlot, 2003, MNRAS
;
; CALLING SEQUENCE:
;   bursts=BURST_GRID_ZTD(met,tdiff,tbc,age_host=,tau_sf=,freq_burst=
;
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
;
; OPTIONAL INPUTS:
;   age_host - [n_model] length vector of the age of the host galaxy.
;              not really needed except for bookkeeping unless
;              /ANCIENT is set
;   tau_sf - the exponential decline rate of the host SFR. not really
;            needed except for bookkeeping.
;   freq_burst - scalar giving the average frequency (in 1/Gyr) of
;                bursts.  Defaults to 1
;   delta_t - scalar giving the age range (in Gyr) of a burst to
;             consider. Defaults to 2 Gyr.  If /ANCIENT is not set,
;             then this gives the maximum age of a burst.  If /ANCIENT
;             is set, then this gives the maximum time after model
;             start that a burst may occur. 
;   band_shift - [n_model] vector of redshifts
;   filters - array of k-correct filter names with which to convolve
;             the spectra
;   attenuate_ionizing - factor by which to reduce the transmitted ionizing flux
;                        (due to absorption by ISM)  Defaults to 0.
;                        Only importnat for magnitudes calculated
;                        blueward of restframe 912AA
;   bc_dust - name of the function that returns tau_lambda/tau_5500
;             given lambda, extra parameters may be passed as extra
;             elements of the array.  This is for the dust that
;             affects young stars
;   bc_diff - name of the function that returns tau_lambda/tau_5500
;             given lambda, extra parameters may be passed as extra
;             elements of the array.  This is for diffuse dust that
;             affects all stars equally
;   burstdir - 
;   burstroot - 
;
; KEYWORDS:
;   vers_bc03 - if set, use the bc03 SSPs.  Default is to use CB07
;   ancient - if set, use burst ages given by
;             age_host-randomu*delta_t, instead of just
;             randomu*delta_t. This can be thought of as approximating
;             bursts of star formation in the very early history of
;             the galaxy (i.e. bulges)  It should be set to something
;             small, <~1 Gyr.  Also, should likely be used in
;             combination with 'MASS' normalization.
;
;---------------

FUNCTION burst_grid_Ztd,met,tdiff,tbc,age_host=age_host,tau_sf=tau_sf,$
                        filters=filters,band_shift=band_shift,$
                        burstroot=burstroot,burstlengths=burstlengths,$
                        vers_sps=vers_sps,sfhtype=sfhtype,imftype=imftype,$
                        delta_t=delta_t,freq_burst=freq_burst,$
                        ancient=ancient,$
                        bc_dust=bc_dust,diff_dust=diff_dust,$
                        attenuate_ionizing=attenuate_ionizing


;;------------
;;SETUP AND DEFAULTS
;;------------
memlimit=5E3 ;maximum number of models to input to quantities_ztd to avoid memory constraints

;SPS version
if keyword_set(vers_sps) EQ 0 then vers_sps='cb07'
if keyword_set(imftype) EQ 0 then imftype='chab'
if keyword_set(sfhtype) EQ 0 then sfhtype='delayed'

if keyword_set(burstdir) EQ 0 then burstdir='$SPECFIT_LIB/'+vers_sps+'/' 
if keyword_set(burstroot) EQ 0 then burstroot=vers_sps+imftype+'_'+sfhtype+'000050' ;50 Myr delayed

zbc=[0.02,0.2,0.4,1.0,2.5]
npz=n_elements(zbc)
nmod=n_elements(met)

;Default filters
gg=['galex_FUV','galex_NUV']+'.par'
ss='sdss_'+['u','g','r','i','z']+'0.par'
tt='twomass_'+['J','H','Ks']+'.par'
ii='spitzer_irac_'+['ch1','ch2','ch3']+'.par'
if n_elements(filters) EQ 0 then filters=[gg,ss,tt,ii]
nfilt=n_elements(filters)


;;burst times and model assignment if not already defined
; 
;GYr  maximum time in the past to consider a burst
if keyword_set(delta_t) EQ 0 then delta_t=2.0
if keyword_set(freq_burst) EQ 0 then freq_burst=1.0
nbursts=freq_burst*delta_t*nmod 
if keyword_set(ancient) then nbursts=nmod ;one bulge per galaxy
print,nbursts,' bursts'
tbursts=randomu(seed,nbursts)*delta_t*1E9  ;age of the burst at the host timepoint
burstmod=long(randomu(seed,nbursts)*nmod)
if keyword_set(ancient) then burstmod=findgen(nmod)


;;----------------------
;;SET UP OUTPUT
;;----------------------

;;make bursts near start of model timeline, not end (i.e., bulges)
if keyword_set(ancient) then tbursts=(age_host[burstmod]-tbursts) > 0


;model assignment
burstnum=dblarr(nbursts)
hh=histogram(burstmod,min=0d,max=nmod-1d,reverse_indices=revind)
nh=n_elements(hh)
ii=revind[0:nh]-(nh+1)
oo=revind[nh+1:*]
maxh=max(hh)
for i=1,maxh do begin
   j=where(hh GE i)
   burstnum[oo[ii[j]+i-1]]=i-1
endfor
nn=max(burstnum)+1 ;or max(hh)
oned_burstindex=burstmod*nn+burstnum 

;output structure
bt={met:0.,tau_sf:0.,tau:0.,tau_bc:0.,age_host:0.,band_shift:0.,$  ;input
    dust_law:' ',dust_law_bc:' ',$                                 ;input
    mag:fltarr(nn,nfilt)+99,mag_nodust:fltarr(nn,nfilt)+99,$
    mbol:fltarr(nn),mstar:fltarr(nn),mgal:fltarr(nn),sfr:fltarr(nn),$
    log_ltrans:fltarr(nn),log_ldust:fltarr(nn),$
    log_ldust_young:fltarr(nn),log_ldust_bc:fltarr(nn),$
    log_lion:fltarr(nn),log_nly:fltarr(nn),log_nly_dust:fltarr(nn),$
    n_burst:0.,burst_weight:fltarr(nn),burst_age:fltarr(nn),$
    max_sfr_burst:fltarr(nn),burst_func:strarr(nn)}

bt=replicate(bt,nmod)
bt.met=met
bt.tau=tdiff
bt.tau_bc=tbc
if n_elements(age_host) GT 0 then bt.age_host=age_host
if n_elements(tau_sf) GT 0 then bt.tau_sf=tau_sf

bt.n_burst=hh
temp=(bt.burst_age)
temp[oned_burstindex]=tbursts
bt.burst_age=temp
temp=bt.burst_weight
temp[oned_burstindex]=1.0
bt.burst_weight=temp

;arrays to hold calculated values temporarily
log_ldust=fltarr(nn,nmod)
log_ltrans=fltarr(nn,nmod)
log_ldust_bc=fltarr(nn,nmod)
log_ldust_young=fltarr(nn,nmod)
log_lion=fltarr(nn,nmod)
sfr=fltarr(nn,nmod)
mstar=fltarr(nn,nmod)
log_nly=fltarr(nn,nmod)
log_nly_dust=fltarr(nn,nmod)
mgal=fltarr(nn,nmod)
mbol=fltarr(nn,nmod)
burst_func=strarr(nn,nmod)

;;------------------
;;LOOP OVER MODELS
;;------------------

;drill down into all models
;in a range of met with the same tsf
;-----------------------------
nsplit=1
tsf=long(randomu(seed,nmod))*nsplit

;should implement isf loop for different burst duration/shapes: but
;then the loops should be over individual bursts, not over
;hosts models...could be done by making 1D arrays before grouping and
;looping.
; really, that would make alot of things easier bookkeeping wise too
usf=tsf[uniq(tsf,sort(tsf))]
npsf=n_elements(usf)
for ipsf=0,npsf-1 do begin
   ;print,ts(ipsf+1) +'of '+ts(npsf)
   this_sf=where(tsf EQ usf[ipsf])
   zz=met[this_sf]
    for iz=0,npz-2 do begin
      this_zr=where(zz GE zbc[iz] and zz LT zbc[iz+1],ntz)
      if keyword_set(verbose) then $
         print,'zindex='+ts(iz),ntz
      if ntz EQ 0 then continue ;no z's in this range
      allsubs=this_sf[this_zr]              ;subs[this_zr]            
;break into pieces to avoid memory constraints
      nblock=ceil(ntz*(freq_burst*delta_t)/memlimit)
      step=ceil(memlimit/(freq_burst*delta_t))
      for iblock=0,nblock-1 do begin
         if keyword_set(verbose) then $
            print,iblock+1,' of ',nblock,iblock*step,' -> ',((iblock+1)*step-1) <(ntz-1)
         subs=allsubs[iblock*step :((iblock+1)*step-1) <(ntz-1)] ;is this right?
         ntblock=n_elements(subs)

;need some tricky bookeeping to get multiple bursts per model
         burstage=transpose(bt[subs].burst_age)
         thisburst=where(burstage GT 0,ntb)*1.0d
         if ntb EQ 0 then continue
         burstz=((bt[subs].met)#(fltarr(nn)+1))[thisburst]
         burst_tbc=((bt[subs].tau_bc)#(fltarr(nn)+1))[thisburst]
         burst_tdiff=((bt[subs].tau)#(fltarr(nn)+1))[thisburst]
         burst_band_shift=((bt[subs].band_shift)#(fltarr(nn)+1))[thisburst]
         burstage=burstage[thisburst] 
         thisburstnum=fix(thisburst/ntblock) ;CHECKED, sort of
         thisburstmod=subs[(thisburst MOD ntblock)]

;read SPS, determine quantities
         mags=quantities_ztd(burstz,burstage,burst_tbc,burst_tdiff,sfh=burstroot,$
                             bc_dust=bc_dust,diff_dust=diff_dust,$
                             band_shift=band_shift,filters=filters,$
                             attenuate_ionizing=attenuate_ionizing,$
                             rundir=burstdir,vers_sps=vers_sps,/zdir,$
                             log_lums=log_lums,bcextras=bcextras,$
                             mags_nodust=mags_nodust)

   
;more tricky bookkeeping:
;figure out where each element of MAGS would go in a nmod X nn X
;nfilt array, and make everything one dimensional
         onedsub2d=(thisburstnum*nmod+thisburstmod)#(fltarr(nfilt)+1)+(fltarr(ntb)+1)#(findgen(nfilt)*nn*nmod)
         mi=reform(mags,ntb*nfilt)
         onedsub_mags=reform(onedsub2d,ntb*nfilt)
         
         temp=transpose(bt.mag,[2,0,1]) ;extract temp as a [nmod,nn,nfilt] array
         temp[onedsub_mags]=mi          ;fill temp with the correct mags using one-d subscript and value arrays
         bt.mag=transpose(temp,[1,2,0]) ;put it back in with the original order of dimensions
         temp=transpose(bt.mag_nodust,[2,0,1])  
         temp[onedsub_mags]=reform(mags_nodust,ntb*nfilt)                 
         bt.mag_nodust=transpose(temp,[1,2,0])  

    ;put in the simpler quantities. 
         onedsub=thisburstmod*nn+thisburstnum
         ;print,n_elements(onedsub),minmax(onedsub)
         mstar[onedsub]=bcextras[*,0]
         mgal[onedsub]=bcextras[*,1]
         sfr[onedsub]=bcextras[*,2]
         mbol[onedsub]=bcextras[*,5]
         log_nly_dust[onedsub]=bcextras[*,4]
         burst_func[onedsub]=burstroot
     ;get dust and lyC photons
         log_nly[onedsub]=bcextras[*,3]
         log_ltrans[onedsub]=log_lums[*,0]
         log_lion[onedsub]=log_lums[*,1]
         log_ldust[onedsub]=log_lums[*,2]
         log_ldust_young[onedsub]=log_lums[*,3]
         log_ldust_bc[onedsub]=log_lums[*,4]
      endfor
   endfor
 endfor

;;----------
;;FILL STRUCTURE INCLUDING UNIT CONVERSIONS
;;----------

bt.mstar=mstar
bt.mgal=mgal
bt.sfr=sfr
bt.mbol=mbol
bt.log_nly_dust=log_nly_dust
bt.log_nly=log_nly
bt.burst_func=burst_func

;convert to absolute magnitudes (given that bc spectra are in L_sol/AA
;and k_correct expects erg/s/cm^2/AA)
pc2cm=3.086E18
K=(3.826E33)/(4.*!PI)/(pc2cm)^2./100. ;distance of 10pc for absolute magnitudes
bt.mag=bt.mag-2.5*alog10(K)
bt.mag_nodust=bt.mag_nodust-2.5*alog10(K)
bt.log_ldust=log_ldust+alog10(K) ;convert L_ir to F_ir in erg/s/cm^2 at 10pc
bt.log_ltrans=log_ltrans+alog10(K)        
bt.log_lion=log_lion+alog10(K)         
bt.log_ldust_young=log_ldust_young+alog10(K)
bt.log_ldust_bc=log_ldust_bc+alog10(K)

return,bt

end
