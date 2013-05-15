;tau_sf=[1,3,5,7,9,12,10]
z=[0.02,0.2,0.4,1.0,2.5]
zname=['m32','m42','m52','m62','m72']
;,'m52','m62','m72']
mu_dust=[0.1]
tau_v=[0.0,1.0]


rundir='~/POP_SYNTH/bc03/runs/cb07/bursts/'
tburst=0.04 ;tau for the delayed burst model in Gyr
oroot='cb07_burst100Myr'
count=0
bcount=0

close,1
openw,1,'csp_grid_burst100'
printf,1,'#'

for iz=0,n_elements(z)-1 do begin
   for imu=0,n_elements(mu_dust)-1 do begin
      for id=0,n_elements(tau_v)-1 do begin

         outname=oroot+$
                 '_z'+strcompress(string(z[iz]*10,format='(I2)'),/remove_all)+$
                 '_mu'+strcompress(string(mu_dust[imu]*10,format='(I2)'),/remove_all)+$
                 '_tv'+strcompress(string(tau_v[id]*10,format='(I2)'),/remove_all)

        if file_test(rundir+outname+'.ised') EQ 0 then begin
           count=count+1
           printf,1,'$bc03/csp_galaxev << !'

           printf,1,'cb2007_hr_stelib_'+zname[iz]+'_chab_ssp.ised'
           printf,1,'Y'         ;add dust
           printf,1,string(tau_v[id],format='(F4.1)')
           printf,1,string(mu_dust[imu],format='(F4.1)')
           printf,1,'2'         ;delayed burst
           printf,1,string(tburst)
;           printf,1,'N'         ; no gas recycling
;           printf,1,'20.0'      ;end after 20 Gyr
;           printf,1,tcut         ;cutoff time for burst, in Gyr
           printf, 1, outname
           printf,1,'!'
;printf,1,'rm e*color*'
;printf,1,'rm e*sindx*'
;printf,1,'rm e*ABmag'
        endif else bcount=bcount+1

     endfor
   endfor
endfor

close,1



end
