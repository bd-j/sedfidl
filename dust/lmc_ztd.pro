
;+
; NAME:
;   LMC_ZTD
;
; PURPOSE: 
;   to return, in the form a_lambda/a_5500, the LMC attenuation/extinction
;   curve as parameterized by  Fitzpatrick 1986? 
;
; CALLING SEQUNECE:
;   a_lambda=lmc_ztd(wave[,pars])
;  
; INPUTS: 
;   WAVE - wavelength vector (angstroms)
;
; OPTIONAL INPUTS:
;   PARS - a vector whose first element gives R_v, the ratio of total
;          to selective extinction (defaults to 3.1) and whose second
;          element is a switch that if true causes the extinction
;          curve for 30 Dor to be used 
;
; OUTPUT: 
;   ALAMBDA - A/A(5500) 
;
;------------------
; from Fitzpatrick 1986, with by default R_V = 3.1 (standard ISM)
; choice of LMC average or 30 Doradus

FUNCTION EXTINCTION_LMC, wave, pars

mic=wave/1E4
if n_elements(pars) GT 0 then pars[0]=R_v else R_v=3.1
if n_elements(pars) GT 1 then if pars[1] then doradus=1
 x = 1. / mic

if keyword_set(doradus) then begin
   c0 = -2.19
   c1 = 1.39
   c2 = 1.49
   c3 = 4.606
   c4 = 0.894
   c5 = 0.43
endif else begin
   c0 = -0.69
   c1 = 0.89
   c2 = 2.55
   c3 = 4.608
   c4 = 0.994
   c5 = 0.50
endelse


ratio_e = c0 + c1 * x + c2 / ((x - c3^2. / x)^2. + c4^2.)
w = where(x ge 5.9, cc)
if cc ge 1 then begin
   y = x[w] - 5.9
   ratio_e[w] = ratio_e[w] + c5 * (0.539 * y^2. + 0.0564 * y^3.)
endif

a_l_v = ratio_e / r_v + 1.

w = where(x lt 3., cc)
if cc ge 1 then for i = 0, cc - 1 do begin
   a_l_v[w] = cardelli_ztd(1/x[w]*1E4,[3.1])
endfor


return, a_l_v

end
