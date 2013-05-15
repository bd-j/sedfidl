;NAME:
;    BC_VALS_NEW
;
;PURPOSE: 
;    To return a set of magnitudes for a Bruzual and Charlot (2003)
;    model, for each time step.  Basically a wrapper on
;    k_project_filters. 
;
;INPUTS:
;    wave - [N_wave] vector of wavelngths in angstroms
;    flux - [N_model x n_wave] array of model spectral fluxes,
;           usually given in L_sun/angstrom/M_sun

FUNCTION bc_vals_new,wave,flux,filters,band_shift=band_shift,d4n=d4n

lambda=k_lambda_to_edges(wave)
;normalize to return absolute magnitudes per solar mass
;conv=alog10(3.826)+33-alog10(4*!PI)-2.0*alog10(3.086)-38
maggies=k_project_filters(lambda,flux,filterlist=filters,band_shift=band_shift,/silent)
mags=0.-2.5*alog10(maggies)
;for i=0,n_elements(d4n)-1 do $
if arg_present(d4n) then d4n=get_d4n_array(lambda,flux)

return,mags

end
