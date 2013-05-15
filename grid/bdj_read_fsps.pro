FUNCTION bdj_read_fsps,file,extras=extras

h=1.054D-27 ;planck constant
J=!lsun/(h*!lightspeed)


str=read_spec(file)
nt=n_elements(str)
nw=n_elements(str[0].lambda)

;stop

fsps={age:fltarr(nt+1),wave:fltarr(nw),flux:fltarr(nw,nt+1)}
extras={logage:fltarr(nt),m_:fltarr(nt),sfr_yr:fltarr(nt),mgalaxy:fltarr(nt),nly:fltarr(nt)}

fsps.wave=str[0].lambda
fsps.flux[*,1:*]=str.spec*(!lightspeed/str.lambda^2) ;convert to L_sun/AA
fsps.age[1:*]=str.agegyr*1E9

extras.logage=alog10(str.agegyr*1E9)
extras.m_=10.^str.logmass
extras.sfr_yr=10.^str.logsfr
extras.mgalaxy=extras.m_ ;HACK.  This is bad. should integrate SFR?

;get the number of ionizing photons from simple integration
wvptr_ion=where(fsps.wave LT 912,nwi)
if nwi GT 0 then begin
   wion=k_lambda_to_edges(fsps.wave[wvptr_ion])
   dwion=wion[1:*]-wion[0:nwi]
   log_Nion=alog10(total((dwion*fsps.wave[wvptr_ion])#(fltarr(nt)+1)*reform(fsps.flux[wvptr_ion,*]),1))+alog10(J)
   extras.nly=log_nion
endif

return,fsps

end
