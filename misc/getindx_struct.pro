FUNCTION getindx_struct,wave,flux,fluxerr=fluxerr,errors=errors
;wrapper on getindx to put results in a structure

nw=n_elements(wave)
if (size(flux))[0] GT 1 then nobj=(size(flux))[2] else nobj=1

readcol,'$SPECFIT_DIR/data/indexlist',i1,i2,b1,b2,r1,r2,type,indname,f='F,F,F,F,F,F,F,A',/silent
indices=getindx(wave,flux,indexdef=indexdef)
;nn=nearest(indexdef[2,*],b1)

nind=n_elements(indname)
tagnames=indname
tagvals=strarr(nind)+'0.0'
if keyword_set(errors) then begin
   tagnames=[tagnames,indname+'_err']
   tagvals=[tagvals,tagvals]
endif

ii={id:''}
ii=struct_addtags(ii,tagnames,tagvals)
ii=replicate(ii,nobj)
for i=0,nind-1 do begin
   ii.(i+1)=reform(indices[i,*])
   if keyword_set(errors) then begin
      ii.(2*i+1)=reform(indices_errors[i,*])
   endif
endfor

return,ii

end
