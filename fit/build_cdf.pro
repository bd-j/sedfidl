;+
; NAME:
;   BUILD_CDF
;
; VERSION:
;   0.0 (Oct, 2008)
;
; PURPOSE:
;   to do stuff  amazing stuff.
;
; REFERENCE:
;   Smiley, G. 2008 ApJ
;
; CATEGORY:
;   Spectral fitting
;
; CALLING SEQUENCE:
;   result=routine(arg1,arg2,arg3,[ARG4=,OUT1=,/SWITCH1, /SWITCH2])
;
; INPUTS:
;    arg1 - the first argument
;    arg2 - the second argument.  It might have a long description in
;       which case the indentation would be like this.  just like this
;       with so many spaces.
;
; OPTIONAL INPUTS:
;    arg3 - the optional input
;
; KEYWORD PARAMETERS:
;    ARG4 - the keyword that must have a value
;    SWITCH1 -  a switch that can be set
;
;
; OUTPUT:
;    Describe the returned result.
;
;    OUT1 - the other returned output
;
; COMMENTS:
;    Describe useful info
;
; REVISION HISTORY:
;    Mar 2008 - written, B. Johnson
;
;--------------------------------------------------------------------

FUNCTION build_cdf,chi2,parval,ord=ord

w=exp(0.-chi2/2d)

ord=sort(parval)

cdf=total(w[ord],/cumulative)
if keyword_set(renormalize) then cdf=cdf/max(cdf)
par=parval[ord]

return,[[par],[cdf]]

end


