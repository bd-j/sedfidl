FUNCTION get_wave_array,vers_sps=vers_sps,imftype=imftype

;if n_elements(vers_sps) EQ 0 then vers_sps='cb07'
;models=file_search('$SPECFIT_LIB/'+vers_sps+'/z100/*delay*mu03_tv00.ised')

;if vers_sps EQ 'cb07' then $
;   cb=

if keyword_set(vers_sps) EQ 0 then vers_sps='cb07'
if keyword_set(imftype) EQ 0 then imftype='chab'
sfhtype='delay'
sfhroot=vers_sps+imftype+'_'+sfhtype+'001000'
cb=bdj_read_sps_zt(1.0,1E9,0.3,0.0,sfhroot,wave=wave,vers_sps=vers_sps,/zdir)

return,wave

end
