FUNCTION get_d4n_array,lambda,flux,flux_ivar=flux_ivar,ivar=ivar

sz=size(flux)
szl=size(lambda)

bluelim=[3850., 3950.]
redlim=[4000., 4100.]

wave=k_lambda_to_edges(lambda)
dlambda=(wave[1:*]-wave[0:szl[1]-2])#(fltarr(sz[2])+1)

iblue=where(wave gt bluelim[0] and wave le bluelim[1], nblue)
ired=where(wave gt redlim[0] and wave le redlim[1], nred)

bluec=total(flux[iblue,*]*dlambda[iblue,*],1)
redc=total(flux[ired,*]*dlambda[ired,*],1)

d4000n=redc/bluec

return, d4000n

end
