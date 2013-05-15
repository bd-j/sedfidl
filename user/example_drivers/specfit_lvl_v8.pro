;;------------------
;; SETUP  - move to a parameter file?
;;------------------
;readcol,paramfile,params,f='A'
;readcol,params[1],filters,f='A'

gg=['galex_FUV','galex_NUV']+'.par'
ss='sdss_'+['u','g','r','i','z']+'0.par'
bb='bessell_'+['U','B','V','R','I']+'.par'
tt='twomass_'+['J','H','Ks']+'.par'
ii='spitzer_irac_'+['ch1','ch2','ch3']+'.par'
filters=[gg,ss,bb,tt,ii]
wl=k_lambda_eff(filterlist=filters)
filterlist=filters
nfilt=n_elements(filterlist)

sp={nmod:3E5,vers_sps:'bc03',imf_type:'chab',$
    nfilt:nfilt,filters:[filterlist],keep_spectra:0,$
    bc_dust:['power_ztd','1.3'],diff_dust:['power_ztd','0.7'],neb_dust:['power_ztd','1.3'],$
    emcorr:1.,redshifts:0.,attenuate_ionizing:0.,$
    bulgesname:'bulgegrid_v8_1',disksname:'diskgrid_v8_1',burstsname:'burstgrid_v8_1',fitname:'lvl_v8_1_v3phot'}

;;-----
;; read the model libraries
;;------
make_library=0
if make_library then generate_grids,sp
;if make_extra_library then extra_grids,sp

bulges=mrdfits(sp.bulgesname+'.fits.gz',1,hdr)
disks=mrdfits(sp.disksname+'.fits.gz',1)
bursts=mrdfits(sp.burstsname+'.fits.gz',1)   

;;-----
;; Add bursts+disks+bulges with weighting and set up for fitting
;;------
;disk_weights=(((1/randomu(seed,sp.nmod)-1.1)<4.0)>0.0)*0.5
bulge_weights = ((randomn(seed,sp.nmod)*0.25+0.3)>0.01) <0.99
disk_weights = 1/bulge_weights-1
bulge_disk = ADD_GRIDS_ZT(bulges,disks,weight_type='MASS',weight_time='CURRENT',weights=disk_weights)
burst_weights = (randomu(seed,sp.nmod)*0.30-0.2)>0
grid = ADD_GRIDS_ZT(bulge_disk,bursts,weight_type='MASS',weight_time='TOTAL',weights=burst_weights)
;burst_weights=(randomu(seed,sp.nmod)*0.10)
;burst_weights=(randomu(seed,sp.nmod)*100)>1
;composite=ADD_GRIDS_ZT(bulge_disk,bursts,weight_type='SFR',weight_time='PEAK',weights=burst_weights)

mwrfits,grid,'combined_v8_1_current.fits',/create

;;-----------
;; Add Nebular Corrections, IR and Ha Data
;;-----------
model_mags = grid.mag   
hii_maggies = hii_maggies(grid,filters);,/attenuate,neb_dust = sp.bc_dust)
model_mags = 0-2.5*alog10(10^(0-model_mags/2.5) + temporary(hii_maggies))
;model_mags=[[model_mags],[0-grid.log_ldust*2.5],[0-(grid.log_nly)*2.5+tau_ha*1.086]]

;;-----------
;;Determine some extra model parameters
;;------------
sfr_avg = 1 
;model_mags=INIT_MODELS_ZT(final,filters=sp.filters,$
;                          bc_dust=sp.bc_dust,diff_dust=sp.diff_dust,$
;                          emcorr=sp.emcorr,sfr_avg=sfr_avg,/verbose)

K=(!lsun)/(4.*!PI)/(!pc2cm)^2./100. ;from Lsun to erg/s/cm^2 at a distance of 10pc for absolute magnitudes
;void=calc_model_pars(grid,sp,bulge_tot=bulge_tot,mass_to_light=m2l_grid)
pars=[[grid.mstar_tot],[grid.sfr_tot],[grid.sfr_tot/grid.mstar_tot],$
      [grid.mag[0]-grid.mag_nodust[0]],[10^(grid.log_ldust-alog10(K))],[grid.met],$
      [grid.tau_bc],[grid.tau],[grid.tau/disks.tau_bc],$
      [grid.age[2]/1e9/grid.tau_sf[2]],[grid.age[2]/1e9],[grid.tau_sf[2]],[1.025*(grid.mstar[2]/bursts.mstar)/grid.mstar_tot],$
      [grid.mstar[0]/grid.mstar_tot]]
parnames=['mass','sfr','ssfr',$
          'Afuv','Ldust','met',$
          'tau_bc','tau_diff','mu',$
          'phase_burst','age_burst','tau_burst','burst_total',$
          'bulge_total']
parnorm=[1.,1.,0.,$
         0.,1.,0.,$
         0.,0.,0.,$
         0.,0.,0.,0.,$
         0]

;;------------------
;;REGULARIZE DATA
;;------------------

init_lvl_data,data_mag,data_magerr,dm=dm,$
              mag_tir=mag_tir,magerr_tir=magerr_tir,$
              mag_ha=mag_ha,magerr_ha=magerr_ha,name=name

nlvl=n_elements(dm)

;;-----------------
;;DO THE FITTING
;;----------------
;gm=findgen(round(sp.nmod/1))
;pars=pars[gm,*]

fit:
t=systime(1) ;get set, go!

for ig=0,nlvl-1 do begin
   if (ig mod 20) EQ 0 then print,ig
   dmag=data_mag[ig,*]
   dmagerr=data_magerr[ig,*]
   mask=(finite(dmag) and (dmagerr GT 0) and finite(dmagerr))
   mask[12:*]=0
   if (total(mask) LT 4) then continue

   mass_out=phot_fit_compress(dmag,dmagerr,model_mags,$
                bandmask=reform(mask),$
                pars=pars,$
                parnames=parnames,$
                parnorm=parnorm,$
                galid=ig,mini_mass=mini_mass)

   out='results/cdfs/'+repstr(sp.fitname,'lvl_','lvlI'+ts(ig)+'_'+name(ig)+'_')+'.fits'
   mwrfits,mass_out,out,/create
   spawn,'gzip -f '+out
   if ig EQ 0 then all_mini=mini_mass else $
      all_mini=[all_mini,mini_mass]
endfor
print,'done in ',systime(1)-t

all_mini = struct_addtags(all_mini,['LDUST_OBSERVED'],['0.0d'])
all_mini.ldust_observed = 10^(0-mag_tir[all_mini.gal_id]/2.5)

;should add a header with all the setup info....
mwrfits,all_mini,'results/'+sp.fitname+'.fits',header_arr,/create
spawn,'gzip -f results/'+sp.fitname+'.fits'

end
