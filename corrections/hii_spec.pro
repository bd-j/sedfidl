
;+
; NAME:
;  HII_SPEC
;
; VERSION:
;   1.0 (Apr, 2010)
;
; PURPOSE:
;   returns the spectrum of an HII region in L_lambda units, given the
;   name of a CLOUDY output fileset.  Treatement of emission lines is
;   in order to give energy conservation when F_lambda is integrated
;   over lambda
;
; REFERENCE:
;   Bruzual & Charlot, 2003, MNRAS
;
; CALLING SEQUENCE:
;   flux=hii_spec(wave=wave,logNly=logNly,cloudyfile=cloudyfile)
;
;
; CATEGORY:
;   Population Synthesis
;
; OPTIONAL INPUTS:
;   wave - the desired wavelength vector for output
;   cloudyfile - string giving the filename for CLOUDY output.  .lint
;                .hyd and .cont will be appended
;
; OUTPUTS:
;   flux - HII region luminosity in L_lambda units for the continuum
;          and lines.  Further more, the units are erg/s/AA for the
;          given input Q(H)
;
;
;
;

FUNCTION hii_spec,wave=wave,logNly=logNly,cloudyfile=cloudyfile

if n_elements(wave) GT 0 then iwave=wave else iwave=0.

;ideally this would result from a CLOUDY run with the host spectrum as
;the input, but this is time consuming to say the least, especailly as
;the CLOUDY results are really only used for emission line corrections
;to the magnitudes.  That said, different HII region parameters are
;available through the use of different cloudy files

if keyword_set(cloudyfile) then fn=cloudyfile else begin
   fn=expand_path('$SPECFIT_DIR/data/')+'/HII_T40k_nH100_U-2.5_Z1.0_hires'
   lognly=47.73  ;Q(H)
   K=alog10(4.*!PI)+17.4771*2.  ;conversion from cloudy units.  i.e., surface area in cm^2 of sphere             
endelse



read_cloudy,fn,energy=energy,name=name,intrinsic=intrinsic,emergent=emergent,wave=wave,diffout=diffout,netout=netout,radius=radius,te=te,eden=eden,hih=hih,hiih=hiih,h2h=h2h,incident=incident

aname=name
dwave=wave

;restrict to requested wavelength range and interpolate onto the
;requested wavelength grid
if iwave[0] GT 0 then begin
   wr=minmax(iwave)
   gd=where(dwave GT wr[0]*0.9 and dwave LT wr[1]*1.1)
   linterp,dwave[gd],diffout[gd],iwave,new_diffout,missing=0.0
   dwave=iwave
   diffout=new_diffout
endif

f_lambda=diffout/dwave
log_flux=alog10(f_lambda)+K  ;cloudy output for diffuse (but not the line list) is divided by sphere area, so need to multiply by this,  this has units of erg/s/AA

dd=k_lambda_to_edges(dwave)
nw=n_elements(dd)
dlambda=abs(dd[0:nw-2]-dd[1:nw-1])

;;--get line fluxes --
;The default line list includes all optical/UV/NIR lines with
;intensitites > I(Ha)/40, and a few extra metal and Balmer series
;lines

linelist=expand_path('$SPECFIT_DIR/data/')+'/line_list'
fmt='(A13,A11,A11)'
close,1
openr,1,linelist
name='' & wave='' & cldyname='' & nn='' & ww='' & cc=''

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

;;
;;--- put lines in diffuse spectrum with correct scaling ----
;;
;Using cloudy output with PunchLwidth=c  (PLW)
;see eq (76) of hazy1
;  the output diffuse has nuF_nu=nuF_nu(cont)+c/PLW*I(line) (erg/s)
;  what k_correct expects is lambda, F_lambda (erg/s/AA/cm^2)
;  log flux is in units of log (erg/s/AA)
;stop

lw=nearest(dwave,linewave)
log_flux_line=log_flux
for iw=0,nline-1 do begin
   if (lw[iw] GT 0 and lw[iw] LT (n_elements(dwave)-1)) then $
      log_flux_line[lw[iw]]=alog10(10.^(log_flux_line[lw[iw]])+$
                                   10.^(lineflux[iw]+alog10(1/dlambda[lw[iw]] - 1/linewave[iw])))
endfor

flux=10.^log_flux_line
oo=sort(dwave)
wave=dwave[oo]
flux=flux[oo]

return, flux

end
