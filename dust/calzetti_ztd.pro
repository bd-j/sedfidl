;+
;
; NAME:
;   CALZETTI_ZTD
;
; PURPOSE:
;   given a wavelength vector lambda, return the Calzetti 2000
;   attenuation curve in terms of a_lambda/a_5500 
;
; CALLING SEQUENCE: 
;   alambda=calzetti_ztd(wave[,pars])
;
; INPUTS: 
;   WAVE - wavelength vector (angstroms)
;
; OPTIONAL INPUTS:
;   PARS - a vector whose first element gives R_v, the ratio of total
;          to selective extinction.  Defaults to 4.05 
;
; OUTPUT: 
;   ALAMBDA - A/A(5500) 
;
; NOTES:
;   extrapolates past the region defined by Calzetti 2000, following
;   the method employed in Hyper-Z, Bolzonella et al.
;-----------

FUNCTION calzetti_ztd,wave,pars
if n_elements(pars) GT 0 then R_v=pars[0] else r_v=4.05
av=1.

;for extrapolations...
p11=1/0.11
ff11=2.659*(-2.156+1.509*p11-0.198*p11^2.+0.011*p11^3.0)+r_v
p12=1/0.12
ff12=2.659*(-2.156+1.509*p12-0.198*p12^2.+0.011*p12^3)+r_v
slope1=(ff12-ff11)/100.
ff99=2.659*(-1.857+1.040/2.19)+r_v
ff100=2.659*(-1.857+1.040/2.2)+r_v
slope2=(ff100-ff99)/100.

;do it
x=1E4/wave
ff= (wave GE 6300. and wave LE 22000)*(2.659*(-1.857+1.040*x)+r_v)
ff+=(wave GE 1200. and wave LT 6300)*(2.659*(-2.156+1.509*x-0.198*x^2.+0.011*x^3.)+r_v)
ff+=(wave LT 1200.)*(ff11+(wave-1100.)*slope1)
ff+=(wave GT 22000.)*(ff99+(wave-21900.)*slope2)

bad=where(ff LT 0, nbad)
if nbad GT 0 then ff[nbad]=0. 

alambda=av*ff/r_v/0.999479
            
return,alambda
        
END
