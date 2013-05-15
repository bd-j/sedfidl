PRO generate_grids,sp;;GENERATE GRIDS FOR FITTING

t=systime(1) ;get set, go!

nfilt=n_elements(sp.filters)
wl=k_lambda_eff(filterlist=sp.filters)
filterlist=sp.filters

;sp={nmod:3E2,vers_sps:'bc03',imf_type:'chab',$
;    nfilt:nfilt,filters:[filterlist],keep_spectra:0,$
;    bc_dust:['power_ztd','1.3'],diff_dust:['power_ztd','0.7'],neb_dust:['power_ztd','1.3'],$
;    emcorr:1.,redshifts:0.,attenuate_ionizing:0.,$
;    bulgesname:'bulgegrid_v8_0',disksname:'diskgrid_v8_0',burstsname:'burstgrid_v8_0',fitname:''}


header_arr=struct2lines_bdj(sp)
nmod=sp.nmod

;;availabale SFHs ;should find a way to pull this out automatically
tau_delay=[10.0d,findgen(10)*20+20,findgen(3)*100+300,findgen(8)*200+600,findgen(10)*500+2500]
tau_exp=[tau_delay,findgen(8)*1000+8000.0d]
tau_delay=tau_delay*1E6
tau_exp=tau_exp*1E6

;;-------------
;;component metallicities
;;-------------
met=randomu(seed,sp.nmod)+0.1   ;metallicities evenly distributed from 0.1 to 1.1 Z_sun
bulge_met=met
disk_met=met
burst_met=met

;;-------------
;;component dust
;;-------------
nodust=fltarr(nmod)

;;------tau_v-----------
y=randomu(seed,nmod) ;for tau_v
line=(findgen(1350)+1)*0.003
c=1.1*line-(alog(cosh(5.5*line-3.0)))/5.5+alog(cosh(-3.0))/5.5 ;match lvl (v0)
c=1.1*line-(alog(cosh(4.5*line-7.0)))/4.5+alog(cosh(-7.0))/4.5 ;match lvl (v1)
c=c/max(c)
j=nearest(c,y)
tv=line[j]
plothist,tv,bin=0.1,/overplot,color=255

;;----- mu -----
;want gaussian around 0.3 with 0.1 scatter
z=randomn(seed,nmod) ;for mu
bar=0.3
sigma=0.1  ;(v2) (v1 value was 0.2)
z=(z*sigma)+bar
bad=where(z LT 0.1 or z GT 1,nbad)
while nbad GT 0 do begin
   zz=randomn(seed,nbad)
   z[bad]=(zz*sigma)+bar
   bad=where(z LT 0.1 or z GT 1,nbad)
endwhile
mu=z

tau_dust_bc=tv*(1-mu)
tau_dust_diff=tv*mu

;;-------------
;;component ages
;;-------------
age_bulge=((randomn(seed,nmod)*0.75+12.25)<13.7)*1E9 ;peak sfr at 10.0-13.0 Gyr lookback
;age_bulge=((randomn(seed,nmod)*1.0+12.0)<13.7)*1E9 ;peak sfr at 10.0-13.0 Gyr lookback
frac_age_disk=1                 ;disks same age as bulge peak SFR
frac_age_disk=(randomu(seed,nmod)*0.5+0.5) ;disks with age from t(bulge_peaksfr) to 0.5*t(bulge_peaksfr)

;;-------------
;;component timescales
;;-------------
tau_bulge=(randomu(seed,nmod)*1.4+0.1d)*1E9    ;tau_bulge from 0.1-1.5 Gyr
tau_disk=(randomu(seed,nmod)*17+3.0d)*1E9       ;tau_disk from 3 to 20 Gyr
tau_burst=(randomu(seed,nmod)*1.9+0.1d)*1E8    ;tau_burst from 10-200 Myr
;tau_burst=fltarr(nmod)+1E8 ;tau_burst=100 Myr ;tau_burst=100Myr

;;snap these to the grid of available values
tau_bulge=tau_delay[nearest(tau_delay,tau_bulge)]/1E9
tau_burst=tau_delay[nearest(tau_delay,tau_burst)]/1E9
tau_disk=tau_exp[nearest(tau_exp,tau_disk)]/1E9
print,minmax(tau_bulge),minmax(tau_disk),minmax(tau_burst)

;;-------------
;;GENERATE THE MODEL GRIDS
;;-------------

;;make bulges with peak sfr at 10.0-13.0 Gyr lookback
;;delta_t_burst=2 Gyr
;;use delayed exponentials
bulges=PHOTO_GRID_ZTD(bulge_met,tau_bulge,tau_dust_diff,nodust,age_bulge+tau_bulge,rundir=rundir,$
                      filters=sp.filters,band_shift=band_shift,keep_spectra=sp.keep_spectra,$
                      vers_sps=sp.vers_sps,sfhtype='delay',imftype=sp.imf_type,$
                      bc_dust=sp.bc_dust,diff_dust=sp.diff_dust,$
                      attenuate_ionizing=sp.attenuate_ionizing)
mwrfits,bulges,sp.bulgesname+'.fits',header_arr,/create
spawn,'gzip -f '+sp.bulgesname+'.fits'



;;make slow declines with ages given by the peak bulge SFR (times 0.5-1.0)
disks=PHOTO_GRID_ZTD(disk_met,tau_disk,tau_dust_diff,tau_dust_bc,(age_bulge)*frac_age_disk,rundir=rundir,$
                     filters=sp.filters,band_shift=band_shift,keep_spectra=sp.keep_spectra,$
                     vers_sps=sp.vers_sps,sfhtype='exp',imftype=sp.imf_type,$
                     bc_dust=sp.bc_dust,diff_dust=sp.diff_dust,$
                     attenuate_ionizing=sp.attenuate_ionizing)
mwrfits,disks,sp.disksname+'.fits',header_arr,/create
spawn,'gzip -f '+sp.disksname+'.fits'

;should add 'tau_burst' to age burst, so that *peaks* are uniformly distributed
bursts=PHOTO_GRID_ZTD(burst_met,tau_burst,tau_dust_diff,tau_dust_bc,randomu(seed,nmod)*1E9+tau_burst,rundir=rundir,$
                      filters=sp.filters,band_shift=band_shift,keep_spectra=sp.keep_spectra,$
                      vers_sps=sp.vers_sps,sfhtype='delay',imftype=sp.imf_type,$
                      bc_dust=sp.bc_dust,diff_dust=sp.diff_dust,$
                      attenuate_ionizing=sp.attenuate_ionizing)


mwrfits,bursts,sp.burstsname+'.fits',header_arr,/create
spawn,'gzip -f '+sp.burstsname+'.fits'

print,systime(1)-t

end


