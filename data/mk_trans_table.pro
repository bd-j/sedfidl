readfast,'IGM_transmission_table.dat',tt,skipline=1
zt=findgen(12)*0.5+1.5
ow=reform(tt[0,*])

bc03=k_im_read_bc03_v3(isedpath='/Volumes/dali/MIRO_DATA/POP_SYNTH/bc03/runs/z10/',isedfile='exp90_z10_mu1_tv0.ised',endian='little')
wave=bc.wave
nw=n_elements(wave)
tr=fltarr(13,nw)
tr[0,*]=wave

red=where(wave GT 1216.)

for i=0,11 do begin
   ww=ow/(1+zt[i])
   linterp,ww,reform(tt[i+1,*]),wave,trest,missing=0.0
   trest[red]=1.0
   tr[i+1,*]=trest
endfor


mwrfits,tr,'meiksin_igmtrans_rest_bc03.fits',/create

bc03=k_im_read_vb07_v3(isedpath='$SPECFIT_LIB/cb07/runs/z10/',isedfile='cb07_exp90_z10_mu1_tv0.ised',endian='little')
wave=bc.wave
nw=n_elements(wave)
tr=fltarr(13,nw)
tr[0,*]=wave

red=where(wave GT 1216.)

for i=0,11 do begin
   ww=ow/(1+zt[i])
   linterp,ww,reform(tt[i+1,*]),wave,trest,missing=0.0
   trest[red]=1.0
   tr[i+1,*]=trest
endfor


mwrfits,tr,'meiksin_igmtrans_rest_cb07.fits',/create

end
