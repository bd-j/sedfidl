;+
; NAME:
;   POWER_ZTD
;
; PURPOSE: 
;   to return, in the form a_lambda/a_5500, a power law attenuation
;   curve 
;    
; CALLING SEQUNECE:
;   a_lambda=power_ztd(wave[,pars])
;  
; INPUTS: 
;   WAVE - wavelength vector (angstroms)
;
; OPTIONAL INPUTS:
;   PARS - a vector whose first element gives alpha, the power law
;          exponent (positive for atttenuation that decreases with wavelength)
;
; OUTPUT: 
;   ALAMBDA - A/A(5500) 
;
;------------------



FUNCTION power_ztd,wave,pars

if n_elements(pars) GT 0 then alpha=pars[0] else alpha=0.7

taulambda=(5500./wave)^(alpha)

return,taulambda

end
