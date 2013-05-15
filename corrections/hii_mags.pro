FUNCTION hii_mags,filterlist=filterlist,cloudyfile=cloudyfile,filt_lam=filt_lam,band_shift=band_shift,wave=wave,flux=flux


if keyword_set(cloudyfile) then fn=cloudyfile else $
fn=expand_path('$SPECFIT_DIR/data/')+'/HII_T40k_nH100_U-2.5_Z1.0_hires'
K=alog10(4.*!PI)+17.4771*2.   ;conversion from cloudy units                  

gg=['galex_FUV','galex_NUV']+'.par'
ss='sdss_'+['u','g','r','i','z']+'0.par'
tt='twomass_'+['J','H','Ks']+'.par'
ii='spitzer_irac_'+['ch1','ch2','ch3']+'.par'

print,keyword_set(filterlist)
if n_elements(filterlist) EQ 0 then filterlist=[gg,ss,tt,ii]


read_cloudy,fn,energy=energy,name=name,intrinsic=intrinsic,emergent=emergent,wave=wave,diffout=diffout,netout=netout,radius=radius,te=te,eden=eden,hih=hih,hiih=hiih,h2h=h2h,incident=incident


aname=name
dwave=wave

f_lambda=diffout/dwave
log_flux=alog10(f_lambda)+K  ;cloudy output for diffuse (but not the lines) is divided by sphere radius, so need to multiply by this

dd=k_lambda_to_edges(dwave)
nw=n_elements(dd)
dlambda=dd[0:nw-2]-dd[1:nw-1]


;;--get line fluxes --
linelist=expand_path('$SPECFIT_DIR/data/')+'/line_list'
fmt='(A13,A11,A11)'
close,1
openr,1,linelist
name=''
wave=''
cldyname=''
nn=''
ww=''
cc=''

while (NOT EOF(1)) do begin
  readf,1,name,wave,cldyname,format=fmt
  nn=[nn,name]
  ww=[ww,wave]
  cc=[cc,cldyname]
endwhile
close,1
cloudynames=cc
mname=name
name=aname


intrinsic_flux=fltarr(n_elements(cloudynames))
name=strtrim(name,2)
for il=0,n_elements(cloudynames) -1 do begin
   gl=(where(name EQ cloudynames[il]))[0]
   if GL GE 0 then intrinsic_flux[il]=intrinsic[gl]
endfor

good=where(intrinsic_flux GT 0 and (strpos(cloudynames,'Ca B') EQ -1))
linewave=float(ww[good])
gg=where(linewave GT 250)  ;remove the IR lines with wavelengths given in micron
linewave=linewave[gg]
lineflux=intrinsic_flux[good[gg]]
nline=n_elements(gg)


;stop

;;
;;--- put lines in diffuse spectrum with correct scaling ----
;;
; Using output with PunchLwidth=c
;see eq (76) of hazy1
;  the output diffuse has nuF_nu=nuF_nu(cont)+c/PLW*I(line)
;  what k_correct expects is lambda, F_lambda

lw=nearest(dwave,linewave)
log_flux_line=log_flux
;for iw=0,nline-1 do log_flux_line[lw[iw]]=30.+alog10(10.^(log_flux[lw[iw]]-30.)+10.^(lineflux[iw]-30.-alog10(dlambda[lw[iw]+1])))
;for iw=0,nline-1 do log_flux_line[lw[iw]]=alog10(10.^(log_flux[lw[iw]])+10.^(lineflux[iw]+alog10(1/dlambda[lw[iw]] - 1/linewave[iw])))
for iw=0,nline-1 do log_flux_line[lw[iw]]=alog10(10.^(log_flux_line[lw[iw]])+10.^(lineflux[iw]+alog10(1/dlambda[lw[iw]] - 1/linewave[iw])))

;plot,dwave,log_flux_line,/xlog,xrange=[0.9E3,12E3]                     
;pc2cm=alog10(3.085)+18
;log_flux_line=log_flux_line-alog10(4.*!PI)-2.*(alog10(10.)+pc2cm) ;erg/s/cm^2/AA at 10pc


;;- normalize by Lha  --
refptr=where(cloudynames[good[gg]] EQ 'H  1  6563A')
log_flux_line=log_flux_line-lineflux[refptr[0]]
flux=10.^log_flux_line
oo=sort(dwave)
wave=dwave[oo]
flux=flux[oo]

filt_lam=k_lambda_eff(filterlist=filterlist)

if keyword_set(band_shift) then begin
   mags =fltarr(n_elements(band_shift),n_elements(filterlist))
   for i=0,n_elements(band_shift)-1 do begin
      maggies=reform(k_project_filters(k_lambda_to_edges(wave),flux,filterlist=filterlist,band_shift=band_shift[i],/silent))
      mags[i,*]=0.-2.5*alog10(maggies)
   endfor
endif else $
   mags=0.-2.5*alog10(k_project_filters(k_lambda_to_edges(wave),flux,filterlist=filterlist,/silent))


mags=reform(mags)
;;- normalize by fha  --  superceded by the above
;ref_ha=lineflux[refptr]-alog10(4.*!PI)-2.*(alog10(10.)+pc2cm)
;ref_ha=ref_ha[0]
;mag_ha=0.-2.5*ref_ha
;mags=mags-mag_ha

;stop

return, mags

end
