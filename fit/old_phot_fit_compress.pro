;+
; NAME:
;   NAME
;
; VERSION:
;   0.0 (Oct, 2008)
;
; PURPOSE:
;   to do stuff  amazing stuff.
;
; REFERENCE:
;   Smiley, G. 2008 ApJ
;
; CATEGORY:
;   Spectral fitting
;
; CALLING SEQUENCE:
;   result=routine(arg1,arg2,arg3,[ARG4=,OUT1=,/SWITCH1, /SWITCH2])
;
; INPUTS:
;    arg1 - the first argument
;    arg2 - the second argument.  It might have a long description in
;       which case the indentation would be like this.  just like this
;       with so many spaces.
;
; OPTIONAL INPUTS:
;    arg3 - the optional input
;
; KEYWORD PARAMETERS:
;    ARG4 - the keyword that must have a value
;    SWITCH1 -  a switch that can be set
;
;
; OUTPUT:
;    Describe the returned result.
;
;    OUT1 - the other returned output
;
; COMMENTS:
;    Describe useful info
;
; REVISION HISTORY:
;    Mar 2008 - written, B. Johnson
;
;--------------------------------------------------------------------
FUNCTION phot_fit_compress,obs_mag,obs_magerr,model_mag,grid_mass,grid_sfr,$
                  galid=galid,inspect=inspect,mini_mass=mini_mass,$
                  pars=pars,parnames=parnames,parnorm=parnorm,$
                  nmag=nmag,magptr=magptr,bandmask=bandmask,$ 
                  plot=plot;,percentiles=percentiles

;window,1

sz=size(obs_mag)

nlvl=sz[1]
nfilt=sz[2]
if keyword_set(magptr) then begin
   mptr=magptr
   nfilt=nmag
endif else mptr=findgen(nfilt)

ngrid=n_elements(grid_mass)

if keyword_set(pars) then begin
   psize=size(pars)
   if psize[0] GT 1 then npars=psize[2] else npars=1
   if psize[1] NE ngrid then begin
      print,'Parameters are not in the right format - N_grid x N_par array required'
      stop
   endif
   pars=transpose(pars)
;   print,keyword_set(parnorm)
   if keyword_set(parnorm) EQ 0 then parnorm=fltarr(npars)  ;default to unnormalized parameters
endif else npars=0
npar=npars

mags=model_mag
if keyword_set(galid) then goo=galid else goo=findgen(nlvl)


for ilvl=0,nlvl-1 do begin


   chi2=find_xi_photo(mags,reform(obs_mag[ilvl,*]),reform(obs_magerr[ilvl,*]),norm=norm,mask=bandmask)
   bchi=min(chi2,imin)
   good=where(chi2-bchi LT 1416d,ngood)  ;restrict to only those models that will have finite likelihood

;stop
print,bchi,ngood

;output structure
   mass_out={$
            gal_id:0.,corr_mags:fltarr(nfilt),magerrs:fltarr(nfilt),$
            mass_median:0.,mass_bestfit:0.,delta_mag_best:fltarr(nfilt),$
                                ;photo_grid_name:' ',$
            chi2:fltarr(ngrid),$
            mass_grid:fltarr(ngrid),mass_cdf:dblarr(2,ngood),$
            sfr_grid:fltarr(ngrid),sfr_cdf:dblarr(2,ngood),$
            ssfr_grid:fltarr(ngrid),ssfr_cdf:dblarr(2,ngood)$
            
            }

   start=n_TAGS(mass_out)


   ns=strcompress(string(long(ngrid)),/remove_all)
   ngs=strcompress(string(long(ngood)),/remove_all)

   if npars GT 0 then $
      mass_out=struct_addtags(mass_out,[parnames+'_grid',parnames+'_cdf'],[strarr(npar)+'fltarr('+ns+')',strarr(npar)+'dblarr(2,'+ngs+')']) ;,mass_out



   mini_mass={$
             gal_id:0.,corr_mags:fltarr(nfilt),magerrs:fltarr(nfilt),$
             delta_mag_best:fltarr(nfilt),chibest:0.,bestfit_index:0l,$
             mass_median:0.,mass_bestfit:0.,mass_97p5:0.,mass_2p5:0.,$
             sfr_median:0.,sfr_bestfit:0.,sfr_97p5:0.,sfr_2p5:0.,$
             ssfr_median:0.,ssfr_bestfit:0.,ssfr_97p5:0.,ssfr_2p5:0. $
                                ;photo_grid_name:' '$
             }

   ministart=n_TAGS(mini_mass)

   if npars GT 0 then $
      mini_mass=struct_addtags(mini_mass,[parnames+'_median',parnames+'_bestfit',parnames+'_2p5',parnames+'_97p5'],strarr(npar*4)+'0.') ;strarr(npar)+'0.',strarr(npar)+'0.',strarr(npar)+'0.']


   ;mass_out=replicate(mass_out,nlvl)
   ;mini_mass=replicate(mini_mass,nlvl)


;stop


  mass_out[ilvl].chi2=chi2;find_xi_photo(mags,reform(obs_mag[ilvl,*]),reform(obs_magerr[ilvl,*]),norm=norm)
  mini_mass[ilvl].chibest=bchi
  mass_out[ilvl].delta_mag_best[mptr]=reform(obs_mag[ilvl,*])-(mags[*,imin]-2.5*alog10(norm[imin]))


  mass_out[ilvl].sfr_grid=grid_sfr*norm
  mass_out[ilvl].mass_grid=grid_mass*norm
  mass_out[ilvl].ssfr_grid=grid_sfr*norm/(grid_mass*norm)
  
;note theat final cdfs will have to multiplied by exp(-bchi/2) to get
;the real value out

  mass_out[ilvl].mass_cdf=transpose(build_cdf($
                      mass_out[ilvl].chi2[good]-bchi,mass_out[ilvl].mass_grid[good]))  
  mass_out[ilvl].sfr_cdf=transpose(build_cdf($
                      mass_out[ilvl].chi2[good]-bchi,mass_out[ilvl].sfr_grid[good]))
  mass_out[ilvl].ssfr_cdf=transpose(build_cdf($
                      mass_out[ilvl].chi2[good]-bchi,mass_out[ilvl].ssfr_grid[good]))


  mini_mass[ilvl].bestfit_index=imin
  mini_mass[ilvl].mass_bestfit=mass_out[ilvl].mass_grid[imin]
  mini_mass[ilvl].sfr_bestfit=mass_out[ilvl].sfr_grid[imin]
  mini_mass[ilvl].ssfr_bestfit=mass_out[ilvl].ssfr_grid[imin]


;GENERIC PARAMETERS CDF
  for ipar=0, npar-1 do begin
     nn=(parnorm[ipar] EQ 0)+(parnorm[ipar] GT 0)*norm 
     mass_out[ilvl].(start+ipar)=pars[ipar,*]*nn  
     mass_out[ilvl].(start+npar+ipar)=transpose(build_cdf($
                      mass_out[ilvl].chi2[good]-bchi,mass_out[ilvl].(start+ipar)[good],ord=ord))
     mini_mass[ilvl].(ministart+npar+ipar)=mass_out[ilvl].(start+ipar)[imin]
  endfor

  j=nearest(mass_out[ilvl].mass_cdf[1,*]/max(mass_out[ilvl].mass_cdf[1,*]),[0.0250,0.50,0.9750])
  if min(j) GE 0 then begin
     mini_mass[ilvl].mass_median=mass_out[ilvl].mass_cdf[0,j[1]]
     mini_mass[ilvl].mass_97p5=mass_out[ilvl].mass_cdf[0,j[2]]
     mini_mass[ilvl].mass_2p5=mass_out[ilvl].mass_cdf[0,j[0]]

     j=nearest(mass_out[ilvl].sfr_cdf[1,*]/max(mass_out[ilvl].sfr_cdf[1,*]),[0.0250,0.50,0.9750])
     mini_mass[ilvl].sfr_median=mass_out[ilvl].sfr_cdf[0,j[1]]
     mini_mass[ilvl].sfr_97p5=mass_out[ilvl].sfr_cdf[0,j[2]]
     mini_mass[ilvl].sfr_2p5=mass_out[ilvl].sfr_cdf[0,j[0]]

     j=nearest(mass_out[ilvl].ssfr_cdf[1,*]/max(mass_out[ilvl].ssfr_cdf[1,*]),[0.0250,0.50,0.9750])
     mini_mass[ilvl].ssfr_median=mass_out[ilvl].ssfr_cdf[0,j[1]]
     mini_mass[ilvl].ssfr_97p5=mass_out[ilvl].ssfr_cdf[0,j[2]]
     mini_mass[ilvl].ssfr_2p5=mass_out[ilvl].ssfr_cdf[0,j[0]]


;GENERIC PARAMETERS PERCENTILES
     for ipar=0,npar-1 do begin
        j=nearest(mass_out[ilvl].(start+npar+ipar)[1,*]/max(mass_out[ilvl].(start+npar+ipar)[1,*]),[0.0250,0.50,0.9750])
        mini_mass[ilvl].(ministart+ipar)=mass_out[ilvl].(start+npar+ipar)[0,j[1]]
        mini_mass[ilvl].(ministart+npar*3+ipar)=mass_out[ilvl].(start+npar+ipar)[0,j[2]]
        mini_mass[ilvl].(ministart+npar*2+ipar)=mass_out[ilvl].(start+npar+ipar)[0,j[0]]
     endfor

     if keyword_set(plot) then begin
        wset,0
        plot,mass_out[ilvl].mass_cdf[0,*],mass_out[ilvl].mass_cdf[1,*],psym=-8,/xlog,xrange=[1E7,1E12]
                                ;oplot,[1,1]*mass_median[ilvl],[0,1E10],linestyle=2
                                ;oplot,[1,1]*mass_best[ilvl],[0,1E10],linestyle=3

     ;stop
     endif
;print,mass_median
  endif else print,'No Median for '+string(goo[ilvl])

  if keyword_set(plot) then begin
     wset,1
     plot,mass_out[ilvl].mass_grid,mass_out[ilvl].chi2,psym=3,/ylog,/xlog,xrange=[1E7,1E12],$
          yrange=[0.8,1000.],title=strcompress('!6'+string(goo[ilvl]),/remove_all),xtitle='SFR',ytitle='!7v!6!u2!n'
     tmp=''
     if keyword_set(inspect) then read,tmp
  endif


endfor

;stop

;;--------
;fill summary structure
;;---------
;mini_mass.photo_grid_name=gridname
mini_mass.gal_id=goo
mini_mass.corr_mags[mptr]=transpose(obs_mag)
mini_mass.magerrs[mptr]=transpose(obs_magerr)

mini_mass.delta_mag_best=mass_out.delta_mag_best


;;------
;;fill cdf structure
;;-------
;if keyword_set(gridname)then mass_out.photo_grid_name=gridname
mass_out.gal_id=goo
mass_out.corr_mags[mptr]=transpose(obs_mag)
mass_out.magerrs[mptr]=transpose(obs_magerr)

mass_out.mass_median=mini_mass.mass_median
mass_out.mass_bestfit=mini_mass.mass_bestfit


return,mass_out

end
