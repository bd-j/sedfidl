PRO cumul_array,array,inds,vals

mx=max(inds)

h1=histogram(inds,reverse_indices=ri1,omin=om)
h2=histogram(h1,reverse_indices=ri2,min=1)

if ri2[1] GT ri2[0] then begin
   vec_inds=ri2[ri2[0]:ri2[1]-1]
   array[om+vec_inds]=vals[ri1[ri1[vec_inds]]]
endif
for j=1,n_elements(h2) do begin
   if ri2[j+1] EQ ri2[j] then continue
   vec_inds=ri2[ri2[j]:ri2[j+1]-1]
   vinds=om+vec_inds
   vec_inds=rebin(ri1[vec_inds],h2[j],j+1,/sample) + $
            rebin(transpose(lindgen(j+1)),h2,[j],j+1,/sample)
   array[vinds]=array[vinds]+total(vals[ri1[vec_inds]],2)
endfor

end
