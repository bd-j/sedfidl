FUNCTION INTIND,lam,func,lo,hi

  ;perform integral over spectrum for index computation, vectorised

  ;---------------------------------------------------------------!
  ;---------------------------------------------------------------!
  nwave=n_elements(lam)
  if (size(func))[0] GT 1 then nspec=(size(func))[2] else nspec=1
  ;take care of the ends
  L1 = (nearest(lam,lo,/below)<(nwave-4))>0
  L2 = (nearest(lam,hi,/below)<(nwave-2))>2
  f1 = reform((func[L1+1,*]-func[L1,*])/(lam[L1+1]-lam[L1])*(lo-lam[L1])+$
              func[L1,*])
  f2 = reform((func[L2+1,*]-func[L2,*])/(lam[L2+1]-lam[L2])*(hi-lam[L2])+$
              func[L2,*])

    IF (L1 EQ L2) THEN $
     intind = (f2+f1)/2.*(hi-lo) $
  ELSE BEGIN
     intind = TOTAL( ((lam[L1+2:l2]-lam[L1+1:L2-1])#(fltarr(nspec)+1))* $
                     (func[L1+2:L2,*]+func[L1+1:L2-1,*])/2.,1)
     intind = intind + (lam[L1+1]-lo)*(f1+func[L1+1,*])/2.
     intind = intind + (hi-lam[L2])*(f2+func[L2,*])/2.
  ENDELSE

return,intind

END 

;!------------------------------------------------------------!
;!------------------------------------------------------------!
;!------------------------------------------------------------!

FUNCTION GETINDX,lambda,spec,index_names=index_names,indexdef=indexdef

  ;routine to calculate indices from an input spectrum
  ;indices are defined in data/allindices.dat
  ;converted from the FSPS Fortran routine.  needs to be vectorised
;!!!!!!! FORTRAN ROUTINE EXPECTS INPUT TO BE AA and F_NU !!!!!!!!
;!!!!!!! BUT, THAT IS PROBABLY WRONG, IT WILL WORK IFF !!!!!!!
; !!!!!! SPECTRA ARE IN AA and F_LAMBDA
;  !---------------------------------------------------------------!
;  !---------------------------------------------------------------!
  nw=n_elements(lambda)
  if (size(spec))[0] GT 1 then nspec=(size(spec))[2] else nspec=1

  ;;LOAD THE INDEX DEFINITIONS
  readfast,'$SPECFIT_DIR/data/allindices.dat',indexdef,header,skipline=8
  nindsps=(size(indexdef))[2]
  indices = fltarr(nindsps,nspec)

  FOR j=0,nindsps-1 DO BEGIN

     ;;!blue continuum
     cb = intind(lambda,spec,indexdef[2,j],indexdef[3,j])
     cb = cb / (indexdef[3,j]-indexdef[2,j])     
     lb = (indexdef[2,j]+indexdef[3,j])/2.

     IF (indexdef[6,j] NE 4) THEN BEGIN
        ;;!red continuum 
        cr = intind(lambda,spec,indexdef[4,j],indexdef[5,j])
        cr = cr / (indexdef[5,j]-indexdef[4,j])
        lr = (indexdef[4,j]+indexdef[5,j])/2.
        ;;!compute integral(fi/fc)
        ;;!NB: fc here is a linear interpolation between the red and blue.
        ;intfifc = intind(lambda,spec/((cr-cb)/(lr-lb)*(lambda-lb) + cb),$
        ;                 indexdef[0,j],indexdef[1,j])

        larr=lambda#(fltarr(nspec)+1)
        cbarr=(fltarr(nw)+1)#cb
        fact=(fltarr(nw)+1)#((cr-cb)/(lr-lb))
        fc=(larr-lb)*fact+cbarr
        ;;vectorised
        intfifc=intind(lambda,spec/fc,indexdef[0,j],indexdef[1,j])
     
     ENDIF 
     IF (indexdef[6,j] EQ 4) THEN BEGIN
        intfifc = intind(lambda,spec,indexdef[1,j],indexdef[2,j])
        intfifc = intfifc/(indexdef[2,j]-indexdef[1,j])
     ENDIF

     IF (indexdef[6,j] EQ 1.) THEN BEGIN
        ;;!compute magnitude
        indices[j,*] = $
           -2.5*ALOG10(intfifc/(indexdef[1,j]-indexdef[0,j]))
     ENDIF 
     IF (indexdef[6,j] EQ 2.) THEN $ ;        !compute EW (in Ang)
        indices[j,*] = (indexdef[1,j]-indexdef[0,j]) - intfifc
     IF (indexdef[6,j] EQ 3.) THEN BEGIN
;        !compute Dn4000
;        !NB: this only works with cr and cb computed above b/c 
;        !the wavelength intervals for blue and red are the 
;        !same for these indices, otherwise a slightly 
;        !different cr and cb would have to be computed
;; N.B. the original formulation for Dn4000 was inconsistent with the
;; indices defined above.  Dn4000 is defined as the ratio of the
;; integrals \int_l1^l2 f_nu d\lambda

        cr = intind(lambda,spec*(lambda^2/!lightspeed#(fltarr(nspec)+1)),indexdef[4,j],indexdef[5,j])
        cb = intind(lambda,spec*(lambda^2/!lightspeed#(fltarr(nspec)+1)),indexdef[2,j],indexdef[3,j])
        ;cr = intind(lambda,spec,indexdef[4,j],indexdef[5,j])
        ;cb = intind(lambda,spec,indexdef[2,j],indexdef[3,j])
        indices[j,*] = cr/cb
     ENDIF     
     IF (indexdef[6,j] EQ 4.) THEN $;!compute magnitude from a flux ratio
        indices[j,*] = -2.5*ALOG10(intfifc/cb)
 
     ;set dummy values for indices defined off of the wavelength grid
     IF (indexdef[4,j] GT lambda[nw-1]) then indices[j,*] = 999.0
     IF (indexdef[2,j] LT lambda[0]) then indices[j,*] = 999.0

  ENDFOR

return,indices

END 
