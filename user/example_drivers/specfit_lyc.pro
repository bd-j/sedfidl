
;;------------------
;; SETUP  - move to a parameter file?
;;------------------
;readcol,paramfile,params,f='A'
;readcol,params[1],filters,f='A'

nmod=3E4
gg=['galex_FUV','galex_NUV']+'.par'
ss='sdss_'+['u','g','r','i','z']+'0.par'
tt='twomass_'+['J','H','Ks']+'.par'
ii='spitzer_irac_'+['ch1','ch2','ch3']+'.par'
filters=[gg,ss,tt,ii]
;attenuate_ionizing=0.0

fc=[0,1,2,3,4,5]
do_lir=1.0
do_ha=1.0
emcorr=1
redshifts=0

;----------
vers_bc03=0. ;use cb07 models
;['power_ztd','1.3'] 
;['power_ztd',1.0]
;['power_ztd','0.7'] 
;['conroy_ztd','3.1','0.8'] 
;['calzetti_ztd','4.05']
;['pei_ztd','2']  ;SMC
;['lmc_ztd','3.1']
neb_dust=['cardelli_ztd','3.1']
star_dust=['cardelli_ztd','3.1']
grid_sfh=[0]

;------------
hostsname='hosts_lyc_v0_0_mw'
newhosts=1.0
;fitname=

;;------------------
;;BUILD MODELS
;;------------------
t=systime(1)

;--------
;;draw from prior distributions and generate host SEDs
;---------
if newhosts then begin
   params=draw_params_Ztdl(nmod)
   tau_star=params[*,1] & tau_neb=params[*,2] & met=params[*,3] & age=params[*,4] 
   tsf=params[*,0] & lyc_abs=params[*,5] & lyc_esc=params[*,6]

   ;for young models include only one kind of dust
   tdiff=tau_star
   tbc=fltarr(nmod)

   hosts=photo_grid_Ztd(met,tsf,tdiff,tbc,age,rundir=rundir,$
                       filters=filters,band_shift=band_shift,$
                       vers_bc03=vers_bc03,$
                       bc_dust=star_dust,diff_dust=star_dust,$
                       attenuate_ionizing=attenuate_ionizing)

   
   eta=tau_neb/tau_star
   mwrfits,hosts,hostsname+'.fits',/create
   spawn,'gzip -f '+hostsname+'.fits'
endif else begin
   hosts=mrdfits(hostsname+'.fits.gz',1)
   nmod=n_elements(hosts)

   ;re-generate distributions for Lyc escape
   ;and absorbtion and for the ratio of tau_neb to tau_stars
   params=draw_params_Ztdl(nmod)
   tau_star=params[*,1] & tau_neb=params[*,2]
   lyc_abs=params[*,5] & lyc_esc=params[*,6]
   eta=tau_neb/tau_star
endelse

;-----
;; add bursts to hosts and set up for fitting
;------
grid=hosts
print,systime(1)-t

model_mags=init_models_Ztdl(grid,lyc_abs,lyc_esc,eta,$
                           grid_mass=grid_mass,grid_pars=grid_pars,$
                           neb_dust=neb_dust,diff_dust=star_dust,$
                           filters=filters,$
                           emcorr=emcorr)

model_mags=model_mags[fc,*]

;;------------------
;;REGULARIZE DATA
;;------------------

init_lyc_data,data_mag,data_magerr,fc=fc,do_ha=do_ha,good=good,do_irmag=do_irmag

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

   out='cdfs/'+repstr(fitname,'lyc_','lyc_'+ts(good[ig])+'_')
   mwrfits,mass_out,out,/create
   spawn,'gzip -f '+out
   if ig EQ 0 then all_mini=mini_mass else $
      all_mini=[all_mini,mini_mass]
endfor


end
