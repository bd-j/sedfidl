;+
; NAME:
;   PEI_ZTD
;
; PURPOSE: 
;   to return, in the form a_lambda/a_5500, the attenuation
;   curves of Pei 1992 
;
; CALLING SEQUNECE:
;   a_lambda=pei_ztd(wave[,pars])
;  
; INPUTS: 
;   WAVE - wavelength vector (angstroms)
;
; OPTIONAL INPUTS:
;   PARS - a vector whose first element gives the choice of MW, LMC,
;          or SMC.  
;
; OUTPUT: 
;   ALAMBDA - A/A(5500) 
;
; HISTORY:
;   bdj - adapted from code of H. Roussel, 04/2011
;
; NOTES:
;   not fully implemented. one would like to
;   make this the default SMC curve.
;
; from Pei 1992
; mic : wavelength in microns
; gal : 'mw', 'lmc' or 'smc'
; output : A(lambda) / A(B)
;---------

FUNCTION PEI_ZTD, wave,pars

mic=wave/1E4
if n_elements(pars[0]) GT 0 then choix=pars[0] else begin
;   print, 'Choose aprs= 0|1|2 (mw|lmc|smc)'
;   goto, fin
choix=2 ;default to SMC
endelse

;case gal of
;'mw': choix = 0
;'lmc': choix = 1
;'smc': choix = 2
;else: begin
;       print, 'Choose gal= ''mw'', ''lmc'' or ''smc''.'
;       goto, fin
;      end
;endcase

a_i = fltarr(3,6)
l_i = fltarr(3,6)
b_i = fltarr(3,6)
n_i = fltarr(3,6)
a_i[0,*] = [165.,  14.,  0.045, 0.002, 0.002, 0.012]
l_i[0,*] = [0.047, 0.08, 0.22,  9.7,   18.,   25.]
b_i[0,*] = [90.,   4.00, -1.95, -1.95, -1.80, 0.00]
n_i[0,*] = [2.0,   6.5,  2.0,   2.0,   2.0,   2.0]

a_i[1,*] = [175.,  19.,  0.023, 0.005, 0.006, 0.020]
l_i[1,*] = [0.046, 0.08, 0.22,  9.7,   18.,   25.]
b_i[1,*] = [90.,   5.50, -1.95, -1.95, -1.80, 0.00]
n_i[1,*] = [2.0,   4.5,  2.0,   2.0,   2.0,   2.0]

a_i[2,*] = [185.,  27.,  0.005, 0.010, 0.012, 0.030]
l_i[2,*] = [0.042, 0.08, 0.22,  9.7,   18.,   25.]
b_i[2,*] = [90.,   5.50, -1.95, -1.95, -1.80, 0.00]
n_i[2,*] = [2.0,   4.0,  2.0,   2.0,   2.0,   2.0]

abs_ab = mic * 0.
mic_5500=5500./1E4
norm=0.
for i = 0, 5 do begin
   norm=norm + a_i[choix,i] / ((mic_5500 / l_i[choix,i])^n_i[choix,i] + $
            (l_i[choix,i] / mic_5500)^n_i[choix,i] + b_i[choix,i])
   abs_ab = abs_ab + a_i[choix,i] / ((mic / l_i[choix,i])^n_i[choix,i] + $
            (l_i[choix,i] / mic)^n_i[choix,i] + b_i[choix,i])
endfor

return, abs_ab/norm

fin:
end
