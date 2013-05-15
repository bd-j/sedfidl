;+
; NAME:
;   CONROY_ZTD
;
; PURPOSE: 
;   to return, in the form a_lambda/a_5500, the attenuation curve as
;   descibed by Conroy et al. 2010, including a reduced UV bump 
;
; CALLING SEQUNECE:
;   a_lambda=conroy_ztd(wave[,pars])
;  
; INPUTS: 
;   WAVE - wavelength vector (angstroms)
;
; OPTIONAL INPUTS:
;   PARS - a vector whose first element gives R_v, the ratio of total
;          to selective extinction (defaults to 3.1) and whose second
;          elemnet gives the 'bump strength' (defaults 0.6)
;
; OUTPUT: 
;   ALAMBDA - A/A(5500) 
;
;------------------

FUNCTION conroy_ztd,lambda,pars

if n_elements(pars) GT 0 then R_v=pars[0] else R_v=3.1
if n_elements(pars) GT 1 then bump=pars[1] else bump=0.6

R=R_v
A_v=1.

x=1./lambda*1E4
nx=n_elements(x)
a=fltarr(nx)
b=a

;IR
ir=where(x GE 0.3 and x LT 1.1,nir)
if nir GT 0 then begin 
   a[ir]=0.574*(x[ir]^1.61)
   b[ir]=(-0.527)*(x[ir]^1.61)
endif

;optical
opt=where(x GE 1.1 and x LT 3.3,nopt)
if nopt GT 0 then begin
   y=x[opt]-1.82
   a[opt] = 1+0.177*y-0.504*y^2-0.0243*y^3+0.721*y^4+0.0198*y^5-0.775*y^6+0.330*y^7
   b[opt]= 1.413*y+2.283*y^2+1.072*y^3-5.384*y^4-0.622*y^5+5.303*y^6-2.090*y^7
endif

;NUV
nuv=where(x  GE 3.3 and x LT 5.9,nn)
if nn GT 0 then begin
   fa=(3.3/x[nuv])^6.*(-0.0370+0.0469*bump-0.601*bump/R+0.542/R)
   a[nuv]=1.752-0.316*x[nuv]-(0.104*Bump/((x[nuv]-4.67)^2+0.341))+fa
   b[nuv]=(-3.09)+1.825*x[nuv]+(1.206*Bump/((x[nuv]-4.62)^2+0.263))
endif

;FUV
fuv=where(x GE 5.9 and x LT 8.0,nf)
if nf GT 0 then begin
   fa= (-0.0447)*(x[fuv]- 5.9)^2.-0.00978*(x[fuv]- 5.9)^3.
   fb= 0.213*(x[fuv]- 5.9)^2.+0.121*(x[fuv]- 5.9)^3.
   a[fuv]=1.752-0.316*x[fuv]-(0.104*Bump/((x[fuv]-4.67)^2+0.341))+fa
   b[fuv]=(-3.09)+1.825*x[fuv]+(1.206*Bump/((x[fuv]-4.62)^2+0.263))+fb
endif

alam=(a+b/R)*A_v

;XUV
xuv=where(x GE 8.0,nx)
if nx GT 0 then begin
   x8=8.0
   fa= (-0.0447)*(x8- 5.9)^2.-0.00978*(x8- 5.9)^3.
   fb= 0.213*(x8- 5.9)^2.+0.121*(x8- 5.9)^3.
   af=1.752-0.316*x8-(0.104*Bump/((x8-4.67)^2+0.341))+fa
   bf=(-3.09)+1.825*x8+(1.206*Bump/((x8-4.62)^2+0.263))+fb

   a8=(af+bf/R)*A_v
   alam[xuv]=(x8/x[xuv])^(-1.3)*a8

endif

;stop

return,alam


end
