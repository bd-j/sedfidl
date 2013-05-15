asciiroot='$SPECFIT_DIR/generate_lib/cb07_scripts/exp'
tau_sf=[0.01,(findgen(10)+1)*0.02,(findgen(10)+1)*0.2,(findgen(10)+5)*0.5,findgen(8)+8]
;tau_sf=[tau_sf,-13,-10,-7,-4]

z=[0.02,0.2,0.4,1.0,2.5]
zname=['m32','m42','m52','m62','m72']
mu_dust=[0.3]
tau_v=[0.0,1.0]

rundir='$SPECFIT_LIB/cb07/'
cb=k_im_read_cb07_v3(isedpath=expand_path('$SPECFIT_DIR')+'/sedlib/cb07/runs/z10/',isedfile='cb07_exp0_z10_mu1_tv0.ised',/silent,endian='little')
ages=cb.age


count=0
bcount=0

for iz=0,n_elements(z)-1 do begin
   zdir='z'+strcompress(string(z[iz]*100,format='(I3.3)'),/remove_all)+'/'
   close,1
   openw,1,'cb07_scripts/cb07salp_csp_grid_onedust_z'+$
         strcompress(string(z[iz]*10,format='(I2)'),/remove_all)
   printf,1,'#'
   for it=0,n_elements(tau_sf)-1 do begin
      
      ;if negative, make the tabulated SFH for BC03
      if tau_sf[it] LT 0 then begin
         na=asciiroot+strcompress(string(tau_sf[it]*1000,format='(I6)'),/remove_all)+'_sfh.ascii'
         if n_elements(sfhfile) GT 0 then sfhfile=[sfhfile,na] else sfhfile=na
         nn=n_elements(sfhfile)
         sfrs=1/abs(tau_sf[it]*1E9)*exp(0.-ages/(tau_sf[it]*1E9))
         close,2
         openw,2,sfhfile[nn-1]
         for iii=0,220 do printf,2,ages[iii],sfrs[iii]
         close,2
      endif

      ;loop over dust
      for imu=0,n_elements(mu_dust)-1 do begin
         for id=0,n_elements(tau_v)-1 do begin

            outname='cb07salp_exp'+$
                    strcompress(string(tau_sf[it]*1000,format='(I6.6)'),/remove_all)+$
                    '_z'+strcompress(string(z[iz]*100,format='(I3.3)'),/remove_all)+$
                    '_mu'+strcompress(string(mu_dust[imu]*10,format='(I2.2)'),/remove_all)+$
                    '_tv'+strcompress(string(tau_v[id]*10,format='(I2.2)'),/remove_all)

            if file_test(rundir+zdir+outname+'.ised') EQ 0 then begin 
               count=count+1
               printf,1,'$bc03/csp_galaxev << !'
               printf,1,'cb2007_hr_stelib_'+zname[iz]+'_salp_ssp.ised'
               printf,1,'Y'     ;add dust
               printf,1,string(tau_v[id],format='(F4.1)')
               printf,1,string(mu_dust[imu],format='(F4.1)')
               if tau_sf[it] LT 0 then begin
                  printf,1,'6'  ;user
                  printf,1,sfhfile[nn-1]
               endif else begin
                  printf,1,'1'  ;exponential decline
                  printf,1,string(tau_sf[it],format='(F4.1)')
                  printf,1,'N'  ; no gas recycling
                  printf,1,'20.0' ;end after 20 Gyr
               endelse
               printf, 1, outname
               printf,1,'!'

            endif else bcount=bcount+1
            
         endfor
      endfor
   endfor

   printf,1,'mv cb07salp_exp*_z'+strcompress(string(z[iz]*100,format='(I2.2)'),/remove_all)+'* '+$
          rundir+zdir
   close,1
endfor

end
