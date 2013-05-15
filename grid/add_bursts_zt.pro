;+
; NAME:
;   ADD_BURSTS_Zt
;
; VERSION:
;   1.0 (Apr, 2010)
;
; PURPOSE:
;   to add together the fluxes and other quantities from structures
;   containg information about hosts models and burst models
;
;
; REFERENCE:
;   Bruzual & Charlot, 2003, MNRAS
;
; CALLING SEQUENCE:
;  grid=add_bursts_Zt(hosts,bursts,weight_type=weight_type,weight_dist=weight_dist,tau_burst=tau_burst,
;                     [scale=scale
;  
;
; CATEGORY:
;   Population Synthesis
;
; INPUTS:
;   hosts - Nmod element structure containing the output from photo_grid_Zt,
;           consisting of the SED and ancillary information for
;           smooth exponentially declining SFRs
;   bursts - Nmod element structure containing the output from burst_grid_Zt,
;            consisting of the SEDs and ancillary information for
;            bursts of a user defined shape.  
;   weight_type - parameter to use for scaling the burst amplitude
;                 relative to the host.  Can be 'MASS' (for stellar
;                 mass) or 'SFR' for instantaneous star formation
;                 rate.  In the case of MASS the scaling is relative
;                 to the final host mass, in the case of SFR it is
;                 relative to the host SFR at the time of the burst
;   weight_dist - distribution of relative scalings for bursts.  Can be
;                 'UNIFORM' or 'LOGARITHMIC'
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
;   tau_burst - defines an offset from the true start of the burst
;               (for example to define the normalization at the *peak*
;               of a delayed burst rather than at the true beginning
;               of the burst).  Defaults to zero
;   scale - overall scaling to apply equally to all bursts, after
;           the relative scaling applied through the weight_type
;           parameters.  Defaults to 1 if not present.
;   logrange - in the case weight_type="LOGARITHMIC', this the
;              logarithmic range for the distribution of scalings,
;              defaults to 2 (i.e., relative scalings distributed
;              logarithmically from 1 to 100)

FUNCTION ADD_BURSTS_Zt,hosts,bursts,$
                       weight_type=weight_type,weight_dist=weight_dist,$
                       scale=scale,logrange=logrange,$
                       tau_burst=tau_burst,burst_peak_sfr=burst_peak_sfr,$
                       burst_tot_mass=burst_tot_mass



if keyword_set(tau_burst) EQ 0 then tau_burst=0 
if n_elements(scale) EQ 0 then scale=1.0
scale=scale[0]

nbursts=total(bursts.n_burst)
nmod=n_elements(hosts)
nn=n_elements(bursts[0].mag[*,0])
nb=nmod*nn
nfilt=n_elements(bursts[0].mag[0,*])

;if the 'hosts' structure is previous output of add_bursts then the
;age field is named differently than if it came from photo_grid_ztd
if tag_exist(hosts,'age') then hostage=hosts.age else hostage=hosts.age_host


;;--------
;;NORMALIZE BURSTS RELATIVE TO HOSTS
case strlowcase(weight_type) of
   'mass': begin
      hostmass=hosts.mgal
      bmass=burst_tot_mass
      normpar=(hostmass/bmass)#(fltarr(nn)+1.) ;this is then relative to the *current, total* host mass, not the mass at the time of the burst.....
      ;stop

   end
   'sfr': begin
      tburst=(fltarr(nn)+1)#hosts.age-bursts.burst_age+tau_burst
      sfr_host=1/(abs((fltarr(nn)+1)#hosts.tau_sf*1E9))*$
               exp(0.-tburst/((fltarr(nn)+1)#hosts.tau_sf*1E9))
      sfr_burst=burst_peak_sfr
      normpar=sfr_host/sfr_burst
   end
endcase

;;----------
;; WEIGHTS
;;---------
case strlowcase(weight_dist) of
   'logarithmic': begin
      if n_elements(logrange) EQ 0 then logrange=2
      weights=reform(10.^(randomu(seed,nb)*logrange),nn,nmod)
   end
  'uniform': weights=reform(randomu(seed,nb),nn,nmod)
endcase



;final weights include renormalization, relative weights, and a flag for a burst
;to actually exist with age less than the age of the host.  
;note that burst frequency can be arbitrarily
;reduced by randomly blanking the burst_age field prior to passing it to
;add_bursts_zt, leading to zero weight
weight=normpar*weights*(bursts.burst_age GT 0)*(bursts.burst_age LT (fltarr(nn)+1)#hostage)
n_burst=total((bursts.burst_age GT 0)*(bursts.burst_age LT (fltarr(nn)+1)#hostage),1)

;if strlowcase(weight_type) EQ 'mass' then stop

;;-------
;COMBINE BURSTS
;;--------
;including relative weighting

;N.b. there is a possible bug (in quantities_ztd?) that causes ldust_bc and ldust_young to
;go NaN for certain ranges of (large) ages.  this is most likely caused by
;floating underflows, and should not be a problem since the young star
;component is insignificant for these ages.  thus, set NaNs and (dust)
;weights to zero

hack=where(finite(bursts.log_ldust_bc)*finite(bursts.log_ldust_young) EQ 0)
junk=bursts.log_ldust_bc
junk[hack]=0.0
bursts.log_ldust_bc=temporary(junk)
junk=bursts.log_ldust_young
junk[hack]=0.0
bursts.log_ldust_young=temporary(junk)
yd_weight=weight 
yd_weight[hack]=0.0 ;make the weights for these points zero when combining young star dust

fluxes=10.^(0.-bursts.mag/2.5)
magweight=rebin(weight,[nn,nmod,nfilt])
magweight=transpose(magweight,[0,2,1])
bmags=0.-2.5*alog10(total(fluxes*magweight,1))
bmag_nd=0.-2.5*alog10(total(10.^(0.-bursts.mag_nodust/2.5)*magweight,1))
bmasses=total(bursts.mstar*weight,1)
bgalmasses=total(bursts.mgal*weight,1)
bsfr=total(bursts.sfr*weight,1)
bdust=alog10(total(10.^(bursts.log_ldust-30.)*weight,1))+30.
bdust_bc=alog10(total(10.^(bursts.log_ldust_bc-30.)*yd_weight,1))+30. ;note the hack in yd_weight
bdust_young=alog10(total(10.^(bursts.log_ldust_young-30.)*yd_weight,1))+30. ;note the hack in yd_weight
btrans=alog10(total(10.^(bursts.log_ltrans-30.)*weight,1))+30.
blion=alog10(total(10.^(bursts.log_lion-30.)*weight,1))+30.
bnly=alog10(total(10.^(bursts.log_nly-30.)*weight,1))+30.


;;---------
; ADD BURSTS TO HOSTS
;;-----------
;next, add the combined bursts to the model grid, all scaled up or
;down the same by 'scale'.  

cmags=0.-2.5*alog10(10.^(0.-hosts.mag/2.5)+scale*10.^(0.-bmags/2.5))
cmags_nd=0.-2.5*alog10(10.^(0.-hosts.mag_nodust/2.5)+scale*10.^(0.-bmags_nd/2.5))
cmasses=hosts.mstar+scale*bmasses
cmtot=hosts.mgal+scale*bgalmasses
csfr=hosts.sfr+scale*bsfr
cdust=alog10(10.^(hosts.log_ldust-30.)+scale*10.^(bdust-30.))+30.
cdust_bc=alog10(10.^(hosts.log_ldust_bc-30.)+scale*10.^(bdust_bc-30.))+30.
cdust_young=alog10(10.^(hosts.log_ldust_young-30.)+scale*10.^(bdust_young-30.))+30.
ctrans=alog10(10.^(hosts.log_ltrans-30.)+scale*10.^(btrans-30.))+30.
clion=alog10(10.^(hosts.log_lion-30.)+scale*10.^(blion-30.))+30.
cnly=alog10(10.^(hosts.log_nly-30.)+scale*10.^(bnly-30.))+30.

;output a combined structure, keeping the burst parameters

grid={met:0.,tau_sf:0.,tau:0.,tau_bc:0.,age_host:0.,band_shift:0.,$
      mag:fltarr(nfilt)+99,mag_nodust:fltarr(nfilt)+99,$
      mstar:0.,mgal:0.,sfr:0.,$
      log_ltrans:0.,log_ldust:0.,log_ldust_young:0.,log_ldust_bc:0.,$
      log_lion:0.,log_nly:0.,log_nly_dust:0.,$
      n_burst:0.,burst_age:fltarr(nn),burst_func:strarr(nn),$
      normpar:fltarr(nn),burst_weight:fltarr(nn)}
grid=replicate(grid,nmod)

if tag_exist(hosts,'Z') then grid.met=hosts.Z else grid.met=hosts.met

grid.tau_sf=hosts.tau_sf
grid.tau=hosts.tau
grid.tau_bc=hosts.tau_bc
grid.band_shift=hosts.band_shift
grid.age_host=hostage 
;grid.age=hostage 


grid.mag=cmags
grid.mag_nodust=cmags_nd
grid.mstar=cmasses
grid.mgal=cmtot
grid.sfr=csfr
grid.log_ldust=cdust
grid.log_ldust_bc=cdust_bc
grid.log_ldust_young=cdust_young
grid.log_lion=clion
grid.log_ltrans=ctrans
grid.log_nly=cnly
grid.n_burst=n_burst
grid.burst_age=bursts.burst_age
grid.burst_func=bursts.burst_func
if max(n_burst) GT 1 then begin
   grid.burst_weight=weight*scale
   grid.normpar=normpar
endif else if max(n_burst) EQ 1 then begin
   grid.burst_weight=reform(weight*scale,1,nmod)
   grid.normpar=reform(normpar,1,nmod)
endif


return, grid
end
