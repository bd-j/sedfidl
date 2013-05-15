;+
; NAME:
;   struct2lines
; PURPOSE:
;   Convert a struct to string array containing keyword-value pairs 
; CALLING SEQUENCE:
;   lines= struct2lines(str)
; INPUTS:
;   str - structure where keywords are tag names, values are values
; OUTPUTS:
;   lines - [N] string array 
; COMMENTS:
;   All structure elements treated as type string
; REVISION HISTORY:
;   May 7, 2008, MRB, NYU
;   Modified to handle fields with multiple elements
;-
;------------------------------------------------------------------------------
function struct2lines_bdj, struct

ntags=n_tags(struct)
lines=''

;lines=strarr(ntags)
names=tag_names(struct)
for i=0L, ntags-1L do begin
   line=strlowcase(names[i])+' '+string(struct.(i)[0])
   lines=[lines,line]
   nel=n_elements(struct.(i))
   if nel GT 1 then begin
      for ie=1l,nel-1L do begin
         line=strlowcase(names[i])+'_'+ts(ie+1)+' '+string(struct.(i)[ie])
         lines=[lines,line]
      endfor
   endif 
endfor

lines=strupcase(lines[1:*])

return, lines

end
