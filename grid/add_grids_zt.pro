;+
; NAME:
;   ADD_GRIDS_Zt
;
; VERSION:
;   1.0 (Apr, 2010)
;
; PURPOSE:
;   to add together the fluxes and other quantities from structures
;   containing information about different grids, line by line matched
;
;
; REFERENCE:
;   Bruzual & Charlot, 2003, MNRAS
;
; CALLING SEQUENCE:
;  grid=add_grids_Zt(hosts,bursts,weights=weights,weight_type)
;  
;
; CATEGORY:
;   Population Synthesis
;
; INPUTS:
;   hosts - Nmod element structure containing the output from
;           photo_grid_Zt or add_grids_Zt, consisting of the SED and
;           ancillary information for smooth exponentially declining SFRs
;   bursts - Nmod element structure containing the output from photo_grid_Zt,
;            consisting of the SEDs and ancillary information for
;            bursts of a user defined shape.  
;   weight_type - parameter to use for scaling the burst amplitude
;                 relative to the host.  Can be 'MASS' (for stellar
;                 mass) or 'SFR' for instantaneous star formation
;                 rate.  
;   weight_time  - time at which to consider the scaling parameter.
;                  can be one of:
;                       CURRENT: scale current burst parameter to
;                       current host parameter (sfr, mass)
;                       TOTAL: scale total burst mass to current host
;                       mass (mass)
;                       PEAK : scale peak burst sfr to host sfr at
;                       time of peak burst (sfr)
;   weights - the relative weights of the different components.  the
;             first component always has a weight of 1
;   burst_peak_sfr - the peak of the SFR of the models used to
;                    generate the bursts.  Used in the normalization
;                    when weight_type='SFR'. E.g., for a BC03 top hat
;                    of 100Myr duration, the burst_peak_sfr is 1E-8.
;   burst_tot_mass - as for burst_peak_sfr, this is used in the
;                    relative normalization for weight_type='MASS'
;                    Gives the total mass formaed in the burst
;                    (usually very close to 1 M_sun)
;
;   
; OPTIONAL INPUTS:

FUNCTION ADD_GRIDS_Zt,hosts,bursts,$
                      weight_type=weight_type,weight_time=weight_time,$
                      weights=weights,$
                      burst_peak_sfr=burst_peak_sfr,$
                      burst_tot_mass=burst_tot_mass

keep_spectra=0
nmod=n_elements(hosts)
nbursts=n_elements(bursts)
nfilt=n_elements(bursts[0].mag[*])
ncomp=n_elements(hosts[0].age)

if keyword_set(tau_burst) EQ 0 then tau_burst=0 
if n_elements(weights) EQ 0 then weights=fltarr(nmod)+1

;if the 'hosts' structure is previous output of add_bursts then the
;age field is named differently than if it came from photo_grid_ztd
if tag_exist(hosts,'age') then hostage=hosts.age else hostage=hosts.age_host


;;-------------
;;NORMALIZE BURSTS RELATIVE TO HOSTS
;;------------


;need total mass formed in the bursts
;need peak sfr of burst
;need current sfr of burst
;need current mass of burst
;need current stellar mass
;need current sfr
;need sfr at time of burst

burst_tot_mass=1.025
burst_peak_sfr=1/bursts.tau_sf
delayed=where(bursts.sfhtype EQ 'delay',ndelay)
if ndelay GT 0 then burst_peak_sfr[delayed]=burst_peak_sfr[delayed]*exp(0.0-1.0)



hostsfr_peak=fltarr(nmod)
for i=0,ncomp-1 do begin
;   if i GT 0 then sc=hosts.normpar[i]*hosts.weight[i] else sc=1 
   if i GT 0 then sc = hosts.scale[i] else sc = 1
   tburst=hosts.age[i]-bursts.age+bursts.tau_sf*(strpos(bursts.sfhtype,'delay') GE 0)
   if total(tburst LE 0) GT 0 and (weight_type EQ 'SFR') and (weight_time EQ 'PEAK') then $
      print,'WARNING!!!: BURSTS ARE OLDER THAN HOST COMPONENTS, DO NOT NORMALIZE BY SFR AT BURST PEAK'
   tmp=1/(abs(hosts.tau_sf[i]*1E9)^2)*tburst*exp(0.-tburst/(hosts.tau_sf[i]*1E9))
   exped=where(hosts.sfhtype[i] EQ 'exp',nexp)
   if nexp GT 0 then tmp[exped]=1/(abs(hosts[exped].tau_sf[i]*1E9))*exp(0.-tburst/(hosts[exped].tau_sf[i]*1E9))
   hostsfr_peak+=tmp*sc
endfor
if tag_exist(hosts,'MSTAR_TOT') then begin
   ;;this is NOT the first component to be added
   host_mstar_tot=hosts.mstar_tot
   host_sfr_tot=hosts.sfr_tot
   host_mgal_tot=hosts.mgal_tot
endif else begin
   ;;this IS the first component to be added
   host_mstar_tot=hosts.mstar
   host_sfr_tot=hosts.sfr
   host_mgal_tot=hosts.mgal
endelse

;weight_type ;MASS or SFR
;weight_time ;CURRENT(mass or sfr) or TOTAL(mass) or PEAK (sfr)
weight_type=strupcase(strcompress(weight_type,/remove_all))
weight_time=strupcase(strcompress(weight_time,/remove_all))

if weight_type EQ 'MASS' then begin
   hostpar=host_mstar_tot ;relative to current host mass
   if weight_time EQ 'CURRENT' then burstpar=bursts.mstar else burstpar=burst_tot_mass
endif else if weight_type EQ 'SFR' then begin
   if weight_time EQ 'CURRENT' then begin
      hostpar=host_sfr_tot ;relative to current host sfr
      burstpar=bursts.sfr[0] ;current burst sfr
   endif else if weight_time EQ 'PEAK' then begin
      tburst=hosts.age-bursts.age+bursts.tau_sf*(strpos(bursts.sfhtype,'delay') GE 0)
      hostpar=hostsfr_peak
      burstpar=burst_peak_sfr ;peak or starting SFR of the burst
   endif
endif
normpar=(hostpar/burstpar)

;stop

;;----------
;; WEIGHTS
;;---------

;final weights include renormalization, relative weights, and a flag for a burst
;to actually exist 
;note that burst frequency can be arbitrarily
;reduced by randomly blanking the age field prior to passing it to
;add_bursts_zt, leading to zero weight

scale=normpar*weights*(bursts.age GT 0);*(bursts.age LT hosts.age) ;where weights is fraction of stellar mass (or sfr) of the first component


;;---------
; ADD BURSTS TO HOSTS
;;-----------
;next, add the  bursts to the model grid, with scaling

cmags=0.-2.5*alog10(10.^(0.-hosts.mag/2.5)+((fltarr(nfilt)+1)#scale)*10.^(0.-bursts.mag/2.5))
cmags_nd=0.-2.5*alog10(10.^(0.-hosts.mag_nodust/2.5)+((fltarr(nfilt)+1)#scale)*10.^(0.-bursts.mag_nodust/2.5))
cdust=alog10(10.^(hosts.log_ldust-30.)+scale*10.^(bursts.log_ldust-30.))+30.
cdust_bc=alog10(10.^(hosts.log_ldust_bc-30.)+scale*10.^(bursts.log_ldust_bc-30.))+30.
cdust_young=alog10(10.^(hosts.log_ldust_young-30.)+scale*10.^(bursts.log_ldust_young-30.))+30.
ctrans=alog10(10.^(hosts.log_ltrans-30.)+scale*10.^(bursts.log_ltrans-30.))+30.
clion=alog10(10.^(hosts.log_lion-30.)+scale*10.^(bursts.log_lion-30.))+30.
cnly=alog10(10.^(hosts.log_nly-30.)+scale*10.^(bursts.log_nly-30.))+30.
if tag_exist(hosts,'spectra') and tag_exist(bursts,'spectra') then begin
   keep_spectra=1
   nw=n_elements(hosts[0].spectra)
   cspec=hosts.spectra+((fltarr(nw)+1)#scale)*bursts.spectra
endif

cmasses=host_mstar_tot+scale*bursts.mstar
cmtot=host_mgal_tot+scale*bursts.mgal
csfr=host_sfr_tot+scale*bursts.sfr

;output a combined structure, keeping some of the burst parameters

grid={met:0.,tau_sf:fltarr(ncomp+1),age:fltarr(ncomp+1),$
      tau:0.,tau_bc:0.,$
      band_shift:0.,$
      mstar:fltarr(ncomp+1),sfr:fltarr(ncomp+1),mgal:fltarr(ncomp+1),$
      mstar_tot:0.,mgal_tot:0.,sfr_tot:0.,$
      mag:fltarr(nfilt)+99,mag_nodust:fltarr(nfilt)+99,$
      log_ltrans:0.,log_ldust:0.,log_ldust_young:0.,log_ldust_bc:0.,$
      log_lion:0.,log_nly:0.,log_nly_dust:0.,$
      sfhtype:strarr(ncomp+1),$
      weight:fltarr(ncomp+1),normpar:fltarr(ncomp+1),scale:fltarr(ncomp+1),$
      weight_type:strarr(ncomp+1),weight_time:strarr(ncomp+1)}
if keep_spectra then $
   grid=struct_addtags(grid,['spectra'],['fltarr('+string(nw,format='(I4.4)')+')'])
grid=replicate(grid,nmod)

;stop

if tag_exist(hosts,'Z') then grid.met=hosts.Z else grid.met=hosts.met
grid.tau=hosts.tau
grid.tau_bc=hosts.tau_bc
grid.band_shift=hosts.band_shift

if ncomp GT 1 then begin
   grid.tau_sf[0:ncomp-1]=hosts.tau_sf
   grid.age[0:ncomp-1]=hosts.age
   grid.mstar[0:ncomp-1]=hosts.mstar
   grid.mgal[0:ncomp-1]=hosts.mgal
   grid.sfr[0:ncomp-1]=hosts.sfr
   grid.sfhtype[0:ncomp-1]=hosts.sfhtype
   grid.weight[0:ncomp-1] = hosts.weight
   grid.scale[0:ncomp-1] = hosts.scale
   grid.normpar[0:ncomp-1] = hosts.normpar
endif else begin
   grid.tau_sf[0:ncomp-1]=reform(hosts.tau_sf,1,nmod)
   grid.age[0:ncomp-1]=reform(hosts.age,1,nmod)
   grid.mstar[0:ncomp-1]=reform(hosts.mstar,1,nmod)
   grid.mgal[0:ncomp-1]=reform(hosts.mgal,1,nmod)
   grid.sfr[0:ncomp-1]=reform(hosts.sfr,1,nmod)
   grid.sfhtype[0:ncomp-1]=reform(hosts.sfhtype,1,nmod)
   grid.normpar[0] = 0.0
   grid.weight[0]=1.0
   grid.scale[0] = 1.0
endelse

grid.tau_sf[ncomp]=bursts.tau_sf
grid.age[ncomp]=bursts.age
grid.mstar[ncomp]=bursts.mstar*scale
grid.mgal[ncomp]=bursts.mgal*scale
grid.sfr[ncomp]=bursts.sfr*scale

grid.weight[ncomp]=weights
grid.sfhtype[ncomp]=bursts.sfhtype
grid.scale[ncomp] = scale
grid.normpar[ncomp]=normpar

grid.mag=cmags
grid.mag_nodust=cmags_nd
grid.log_ldust=cdust
grid.log_ldust_bc=cdust_bc
grid.log_ldust_young=cdust_young
grid.log_lion=clion
grid.log_ltrans=ctrans
grid.log_nly=cnly


if keep_spectra then grid.spectra=cspec

grid.mstar_tot=cmasses
grid.mgal_tot=cmtot
grid.sfr_tot=csfr

;if max(n_burst) GT 1 then begin
;   grid.burst_weight=weight*scale
;   grid.normpar=normpar
;endif else if max(n_burst) EQ 1 then begin
;   grid.burst_weight=reform(weight*scale,1,nmod)
;   grid.normpar=reform(normpar,1,nmod)
;endif


return, grid
end
