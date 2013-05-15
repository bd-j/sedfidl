tau_sf=[0.01,$
        (findgen(9)+1)*0.02,$
        (findgen(5)+2)*0.1,$
        (findgen(7)+4)*0.2,$
        (findgen(10)+5)*0.5]    ;,findgen(8)+8]
;tau_sf=[tau_sf,-13,-10,-7,-4]

imf=['salp','chab']
z=[0.02,0.2,0.4,1.0,2.5]
zname=['m32','m42','m52','m62','m72']
mu_dust=[0.3]
tau_v=[0.0,1.0]

rundir='$SPECFIT_LIB/bc03/'

count=0
bcount=0


for iz=0,n_elements(z)-1 do begin
   zdir='z'+strcompress(string(round(z[iz]*100),format='(I3.3)'),/remove_all)+'/'
   scriptname='bc03_scripts/script_bc03allimf_grid_onedust_delay_z'+$
              strcompress(string(round(z[iz]*100),format='(I3.3)'),/remove_all)
   close,1
   openw,1,scriptname
   printf,1,'#'

   for im=0,1 do begin
      for it=0,n_elements(tau_sf)-1 do begin
         ;;delayed has no negative tau_sf
         ;;loop over dust
         for imu=0,n_elements(mu_dust)-1 do begin
            for id=0,n_elements(tau_v)-1 do begin

               outname='bc03'+imf[im]+'_delay'+$
                    strcompress(string(round(tau_sf[it]*1000),format='(I6.6)'),/remove_all)+$
                    '_z'+strcompress(string(round(z[iz]*100),format='(I3.3)'),/remove_all)+$
                    '_mu'+strcompress(string(round(mu_dust[imu]*10),format='(I2.2)'),/remove_all)+$
                    '_tv'+strcompress(string(round(tau_v[id]*10),format='(I2.2)'),/remove_all)

               if file_test(rundir+zdir+outname+'.ised') EQ 0 then begin 
                  count=count+1
                  printf,1,'$bc03/csp_galaxev << !'
                  printf,1,'bc2003_hr_'+zname[iz]+'_'+imf[im]+'_ssp.ised'
                  printf,1,'Y'  ;add dust
                  printf,1,string(tau_v[id],format='(F4.1)')
                  printf,1,string(mu_dust[imu],format='(F4.1)')
                  printf,1,'4'  ;delayed burst
                  printf,1,string(tau_sf[it],format='(F5.2)')
                  printf,1,'N'  ; no gas recycling
                  printf,1,'20.0' ;end after 20 Gyr
                  printf, 1, outname
                  printf,1,'!'

               endif else bcount=bcount+1
               
            endfor
         endfor
      endfor
   endfor

   printf,1,'mv bc03*_delay*_z'+strcompress(string(round(z[iz]*100),format='(I3.3)'),/remove_all)+'* '+$
          rundir+zdir
   close,1

endfor


end
