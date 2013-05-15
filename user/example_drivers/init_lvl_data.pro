;;---------------------
;;--DATA SETUP
;;----------------------

PRO init_lvl_data,mag_lvl,magerr_lvl,fc=fc,do_ha=do_ha,do_irmag=do_irmag,good=good

;;--read data input
lvlphot='/Users/bjohnson/SFR_FIELDS/Nearby/LVL/sdss_photom/sdssphot.lvl.v3_3_mms.fits'
lvl=mrdfits(lvlphot,1)
linfo=mrdfits('~/SFR_FIELDS/Nearby/LVL/wikidata_0808/info.lvl_w0808.fits',1)
lvl_ir=mrdfits('~/SFR_FIELDS/Nearby/LVL/wikidata_0808/irphot.lvl_w0808.fits',1)
lvl_ha=mrdfits('~/SFR_FIELDS/Nearby/LVL/wikidata_0808/halpha.lvl_w0808.fits',1)
lvl_uv=mrdfits('~/SFR_FIELDS/Nearby/LVL/wikidata_0808/uvphot.lvl_daleaps.fits',1)

dl07=mrdfits('~/SFR_FIELDS/Nearby/LVL/fits/irsed/lvl_dl07_mini.fits.gz',1)
ha=lvl_ha

;;--set up data array
filters='sdss_'+['u','g','r','i','z']+'0.par'
filters=['galex_FUV.par','galex_NUV.par',filters]
tt='twomass_'+['J','H','Ks']+'.par'
ii='spitzer_irac_'+['ch1','ch2','ch3']+'.par'
filters=[filters,tt,ii]
filters=filters[fc]
nfilt=n_elements(filters)

mag_lvl=[[lvl_uv.fuv_mag],[lvl_uv.nuv_mag],[lvl.u_mag],[lvl.g_mag],[lvl.r_mag],[lvl.i_mag],[lvl.z_mag],[lvl_ir.J_mag],[lvl_ir.H_mag],[lvl_ir.K_mag],[lvl_ir.I1_mag],[lvl_ir.I2_mag]]
magerr_lvl=[[sqrt(lvl_uv.fuv_magerr^2-0.05^2)],[sqrt(lvl_uv.nuv_magerr^2-0.03^2.)],[lvl.u_magerr*8],[lvl.g_magerr*2],[lvl.r_magerr*2],[lvl.i_magerr*2],[lvl.z_magerr*2],[lvl_ir.J_magerr],[lvl_ir.H_magerr],[lvl_ir.K_magerr],[lvl_ir.I1_magerr],[lvl_ir.I2_magerr]]

;mag_lvl=mag_lvl[*,fc]
;magerr_lvl=magerr_lvl[*,fc]

;;--add Halpha
;if keyword_set(do_ha) then begin
   haf=10.^(lvl_ha.halpha_flux)  ;cgs
   haf=haf*(1-lvl_ha.nii_halpha)
   ;haf=reform(haf,1,n_elements(haf)
   mag_lvl_ha=0.-2.5*alog10(haf)
   magerr_ha=2.5*ha.halpha_flux_err

   ul=where(ha.halpha_code EQ -1)
   mag_lvl_ha[ul]=mag_lvl_ha[ul]+2.5*alog10(5.)
   magerr_ha[ul]=1.086

   mag_lvl=[[mag_lvl],[mag_lvl_ha]]
   magerr_lvl=[[magerr_lvl],[magerr_ha]]
   nfilt=nfilt+1
   ;fc=[fc,13]
;endif

;;--add IR
;if keyword_set(do_irmag) then begin
   
;24 micron based tirf
   c=2.998E18  ;AA/s
   loglsol=alog10(3.839)+33
   log_mpc_to_cm=alog10(3.086)+18+6
   tirf=alog10((c/(k_lambda_eff(filterlist=['spitzer_mips_24.par']))[0])*3631.*10.^(0.-0.4*lvl_ir.m1_mag)*(1/0.074))  ;Jy*Hz
   
;add the tir from draine+li fitting if it exists
   tirf[dl07.gal_id]=alog10(dl07.f_tir_median)+loglsol+23-alog10(4.*!PI)-2.0*alog10(linfo[dl07.gal_id].dist)-2.*log_mpc_to_cm

   mag_tir=0.-2.5*tirf+2.5*23.  ;cgs
   magerr_tir=sqrt(lvl_ir.m1_magerr^2.+0.4^2.)
   magerr_tir[dl07.gal_id]=sqrt(lvl_ir[dl07.gal_id].m1_magerr^2.+0.2^2.)

   ul=where(lvl_ir.m1_magerr LE 0)  ;upper limits
   magerr_tir[ul]=1.086   ;set flux_err=flux
   mag_tir[ul]=mag_tir[ul]+2.5*alog10(3.)  ;convert from 3sigma to 1 sigma flux

   mag_lvl=[[mag_lvl],[mag_tir]]
   magerr_lvl=[[magerr_lvl],[magerr_tir]]
   nfilt=nfilt+1
   ;fc=[fc,14]
;endif


;;--do distance modulus correction
dm=5.0*alog10(linfo.dist)+25. ;where dist is in Mpc
dmarr=dm#(fltarr(nfilt)+1)
mag_lvl=mag_lvl-dmarr

;;--do MW extinction correction

;f,n are 8.24,8.2, but are already corrected?
r=[0.0,0.0,5.16,3.79,2.75,2.126,1.505,0.9,0.56,0.36,0.0,0.0,0.0,2.5,0.] ; 2nd to last one is for h-alpha, need to get value ;; got it from odonnell - other values are from ccm
a=lvl_ir.ebv#r[fc]
mag_lvl=mag_lvl-a

;;--add systematics to the errors
;calibration errors of GALEX+SDSS+Ha
syserr=([0.04,0.02,0.03,0.015,0.015,0.015,0.03,0.0,0.0,0.0,0.,0.,0.,0.1,0.00])
;estimate of error dus to aperture choices, fg stars, etc....
aperr=(0.*[0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05])
magerr_lvl=sqrt(magerr_lvl^2.+((fltarr(nlvl)+1)#syserr)^2.+((fltarr(nlvl)+1)#aperr)^2.)
magerr_lvl=magerr_lvl < 1.086

;;--cut on good magnitudes
f1=(mag_lvl GT 0)  ;require sdss+galex flux measurement
f2=(magerr_lvl GE 0) 
f3=(finite(mag_lvl))

flag=total(f1,2)+total(f2,2)+total(f3,2)
good=where(flag EQ (nfilt)*3)
nlvl=n_elements(good)

mag_lvl=mag_lvl[good,*]
magerr_lvl=magerr_lvl[good,*]
lvl=lvl[good]
linfo=linfo[good]
lvl_ir=lvl_ir[good]
lvl_ha=lvl_ha[good]
end

