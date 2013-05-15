;there is a faster array based way to code this.  Or the NR version
;using bisection....

FUNCTION nearest,array,value,above=above,below=below

j=where(finite(array),count)
if count EQ 0 then message,'no finite values of array'
array=array[j]

nv=n_elements(value)
mm=dblarr(nv)

for iv=0l,nv-1 do begin
  diff=double(array-value[iv])
  nr=n_elements(array)
  restrict=lindgen(nr)

  if keyword_set(above) then restrict=where(diff GE 0,nr) ;nearest value above the target
  if keyword_set(below) then restrict=where(diff LE 0,nr) ;nearest value below the target
  if nr EQ 0 then restrict=lindgen(n_elements(array))
  mm[iv]=(restrict[where(abs(diff[restrict]) EQ min(abs(diff[restrict])))])[0]
;stop
endfor

if nv EQ 1 then mm=mm[0]

return,j[mm]

end
