;+
; NAME:
;   CARDELLI_ZTD
;
; PURPOSE: 
;   to return, in the form a_lambda/a_5500, the Milky Way attenuation
;   curve as parameterized by Cardelli, Clayton and Mathis 1988? 
;
; CALLING SEQUNECE:
;   a_lambda=cardelli_ztd(wave[,pars])
;  
; INPUTS: 
;   WAVE - wavelength vector (angstroms)
;
; OPTIONAL INPUTS:
;   PARS - a vector whose first element gives R_v, the ratio of total
;          to selective extinction (defaults to 3.1)
;
; OUTPUT: 
;   ALAMBDA - A/A(5500) 
;
;------------------

FUNCTION CARDELLI_ZTD, wave,pars


mic=wave/1E4
if keyword_set(pars) then R_v=pars[0] else R_v=3.1

x_sup = 10.
x_inf = 0.3
x = 1. / mic
a = x * 0.
b = x * 0.

w1 = where(x ge 1.1 and x le 3.3, cc1)
w2 = where(x ge x_inf and x lt 1.1, cc2)
w3 = where(x gt 3.3 and x le 8., cc3)
w4 = where(x gt 8. and x le x_sup, cc4)
wsh = where(x gt x_sup, ccsh)
wlg = where(x lt x_inf,cclg)
if (cc1 ge 1) then begin
   y = x(w1) - 1.82
   a(w1) = 1. + 0.17699 * y - 0.50447 * y^2. - 0.02427 * y^3. + 0.72085 * y^4. $
          + 0.01979 * y^5. - 0.77530 * y^6. + 0.32999 * y^7.
   b(w1) = 1.41338 * y + 2.28305 * y^2. + 1.07233 * y^3. - 5.38434 * y^4. $
          - 0.62251 * y^5. + 5.30260 * y^6. - 2.09002 * y^7.
endif
if (cc2 ge 1) then begin
   y = (x(w2))^1.61
   a(w2) = 0.574 * y
   b(w2) = -0.527 * y
endif
if (cc3 ge 1) then begin
   fa = x(w3) * 0.
   fb = x(w3) * 0.
   ou = where(x(w3) gt 5.9, cou)
   if (cou ge 1) then begin
      y = x(w3(ou)) - 5.9
      fa(ou) = -0.04473 * y^2. - 0.009779 * y^3.
      fb(ou) = 0.2130 * y^2. + 0.1207 * y^3.
   endif
   a(w3) = 1.752 - 0.316 * x(w3) - 0.104 / ((x(w3) - 4.67)^2. + 0.341) + fa
   b(w3) = -3.090 + 1.825 * x(w3) + 1.206 / ((x(w3) - 4.62)^2. + 0.263) + fb
endif
if (cc4 ge 1) then begin
   y = x(w4) - 8.
   a(w4) = -1.073 - 0.628 * y + 0.137 * y^2. - 0.070 * y^3.
   b(w4) = 13.670 + 4.257 * y - 0.420 * y^2. + 0.374 * y^3.
endif
;if (ccsh ge 1) then begin
;   ;extrapolate
;   y=[0.3,0.4]-1.82
;   aa=0.574*y*r_v
;   bb=-0.527*y
   
   ;print, 'wavelength domain not implemented !'
   ;goto, fin
;endif

abs_av = a + b / r_v

return, abs_av

fin:
end
