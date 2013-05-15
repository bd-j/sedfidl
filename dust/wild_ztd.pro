;+
;
; NAME:
;   WILD_ZTD
;
; PURPOSE:
;   given a wavelength vector lambda, return the Wild 2010
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

FUNCTION wild_ztd,wave,pars

nw=n_elements(wave)
wave=[wave,5500]

if n_elements(pars) GT 0 then R_v=pars[0] else r_v=4.05
isb=0  ;mu_* < 3x10^8
x=0 ;x=b/a-0.6
y=0 ;y=(log sSFR)+9.5

n=20 ;softening

lb_eff=[2175.0d,3000.0d,8000.0d]
c=[[0.15,0.00,0.20,0.70,0.00,0.40,1.10,0.40,-0.10],$
   [0.20,0.90,0.10,1.10,0.30,-0.20,1.30,0.60,0.10]]
c=reform(c[isb,*])
s=fltarr(4)
s[0]=c[0]+c[1]*x+c[2]*y
s[1]=c[3]+c[4]*x+c[5]*y
s[2]=c[6]+c[7]*x+c[8]*y
s[3]=1.6
lb1=lb_eff[0]
lb2=(lb1^s[1]/lb_eff[1]^(s[1]-s[2]))^(1/s[2])
lb3=(lb2^s[2]/lb_eff[2]^(s[2]-s[3]))^(1/s[3])

q=( (wave/lb1)^(n*s[0])+(wave/lb1)^(n*s[1])+(wave/lb2)^(n*s[2])+$
    (wave/lb3)^(n*s[3]) )^(0-1.0/n)

alambda=(q/q[nw])[0:nw-1]     
wave=wave[0:nw-1]   

return,alambda
        
END
