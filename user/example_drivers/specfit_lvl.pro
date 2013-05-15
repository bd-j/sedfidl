
;;------------------
;; SETUP  - move to a parameter file?
;;------------------
;readcol,paramfile,params,f='A'
;readcol,params[1],filters,f='A'

nmod=1E5
gg=['galex_FUV','galex_NUV']+'.par'
ss='sdss_'+['u','g','r','i','z']+'0.par'
tt='twomass_'+['J','H','Ks']+'.par'
ii='spitzer_irac_'+['ch1','ch2','ch3']+'.par'
filters=[gg,ss,tt,ii]
attenuate_ionizing=1.0

fc=[0,1,2,3,4,5]
do_lir=1.0
do_ha=1.0
redshifts=0

;----------
vers_bc03=0.
bc_dust=['power_ztd','1.3'] ;['power',1.0]
diff_dust=['power_ztd','0.7'] ;['conroy','3.1','0.8'] ;['calzetti',4.05]
grid_sfh=[1,3,5,7,9,12,-5,-1]
;rundir=
;burstdir=
burstroot='cb07_burst40MyrDelay'
delta_t=2.0
freq_burst=1.0
tau_burst=40E6
burst_peak_sfr=1/tau_burst*exp(0.-1)
burst_tot_mass=1.025

;-------------
weight_type='sfr';'mass'
weight_dist='logarithmic' ;'uniform'
logrange=2.
scale=1.0

;------------
hostsname='hostgrid_v6_0'
burstsname='burstgrid_v6_0'
newhosts=0.0
newbursts=1.0
;fitname=

;;------------------
;;BUILD MODELS
;;------------------

t=systime(1)
;--------
;;draw from prior distributions and generate host SEDs
;---------
if newhosts then begin
   params=draw_params_Ztd(nmod,tausf=tausf,tauv=tauv,mu=mu,met=met,$
                         grid_sfh=grid_sfh,grid_tauv=grid_tauv,grid_mu=grid_mu,$
                         redshifts=redshifts)
   tsf=params[*,0] & tv=params[*,1] & mu=params[*,2] & met=params[*,3] & age=params[*,5] ;&& band_shift=params[*,4]


   tbc=tv*(1-mu)
   tdiff=tv*mu
   hosts=photo_grid_Ztd(met,tsf,tdiff,tbc,age,rundir=rundir,$
                       filters=filters,band_shift=band_shift,$
                       vers_bc03=vers_bc03,$
                       bc_dust=bc_dust,diff_dust=diff_dust,$
                       attenuate_ionizing=attenuate_ionizing)

   mwrfits,hosts,hostsname+'.fits',/create
   spawn,'gzip -f '+hostsname+'.fits'
endif else  hosts=mrdfits(hostsname+'.fits.gz',1)

;-----
;;generate burst SEDs
;------
if newbursts then begin
   tdiff=hosts.tau & tbc=hosts.tau_bc & met=hosts.Z & tsf=hosts.tau_sf & age=hosts.age ;& band_shift=hosts.band_shift  ; note that burst paroperties can be made differnt than the hosts in this step
   bursts=burst_grid_Ztd(met,tdiff,tbc,age_host=age_host,tau_sf=tsf,$
                        filters=filters,band_shift=band_shift,$
                        burstdir=burstdir,burstroot=burstroot,$
                        delta_t=delta_t,freq_burst=freq_burst,$
                        bc_dust=bc_dust,diff_dust=diff_dust,$
                        attenuate_ionizing=attenuate_ionizing,$
                        vers_bc03=vers_bc03)

   mwrfits,bursts,burstsname+'.fits',/create
   spawn,'gzip -f '+burstsname+'.fits'
endif else bursts=mrdfits(burstsname+'.fits.gz',1)


;generate bulges????? special case of old bursts.....
;this would be kind of cool. could have differnt attenuation properties

;-----
;; add bursts to hosts and set up for fitting
;------
grid=add_bursts_Zt(hosts,bursts,weight_type=weight_type,weight_dist=weight_dist,scale=scale,tau_burst=tau_burst,burst_peak_sfr=burst_peak_sfr,burst_tot_mass=burst_tot_mass)

print,systime(1)-t
sfr_avg=1 ;this will get filled with an array of SFR averaged over different timescales
model_mags=init_models_Zt(grid,do_lir=do_lir,do_ha=do_ha,filters=filters,$
                          grid_sfr=grid_sfr,grid_mass=grid_mass,sfr_avg=sfr_avg,$
                          grid_tauv=grid_tauv,grid_mu=grid_mu,grid_age=grid_age,$
                          bc_dust=bc_dust,diff_dust=diff_dust,$
                          burst_function=burst_function,emcorr=emcorr,scale=scale)

model_mags=model_mags[fc,*]
sfr7=sfr_avg[*,0]
sfr8=sfr_avg[*,1]
rsfr=sfr8/sfr7

;;------------------
;;REGULARIZE DATA
;;------------------


init_lvl_data,data_mag,data_magerr,fc=fc,do_ha=do_ha,good=good,do_irmag=do_irmag

;;-----------------
;;DO THE FITTING
;;----------------

for ig=0,ng-1 do begin
   mass_out=phot_fit_compress(data_mag[ig,*],data_magerr[ig,*],model_mags,$
                grid_mass,sfr8,$
                pars=[[grid_sfr],[sfr7],[grid_tauv],[grid_mu],[grid_z],[rsfr]],$
                parnames=['sfr_inst','sfr7','tauv','mu','met','rat_sfr'],$
                parnorm=[1.,1.,0.,0.,0.,0.],$
                galid=good[ig],mini_mass=mini_mass)

   out='cdfs/'+repstr(fitname,'lvl_','lvl_'+ts(good[ig])+'_')
   mwrfits,mass_out,out,/create
   spawn,'gzip -f '+out
   if ig EQ 0 then all_mini=mini_mass else $
      all_mini=[all_mini,mini_mass]
endfor



end
