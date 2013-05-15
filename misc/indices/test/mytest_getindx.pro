files=file_search('./*fits')
nf=n_elements(files)

i=0
allindices=fltarr(nf,26)

for i=0,nf-1 do begin
   spec=mrdfits(files[i],0,hdr,/silent)
   nw=n_elements(spec)
   wave=(dindgen(nw)-sxpar(hdr,'CRPIX1')+1)*sxpar(hdr,'CDELT1')+sxpar(hdr,'CRVAL1')
   wave=wave*1E10
   indices=getindx_struct(wave,spec)
    if i EQ 0 then print,'name ';,string(indexdef[0,*])
   ;print,strcompress(sxpar(hdr,'OBJECT'),/remove_all);,indices[0:25]
   indices.id=strcompress(sxpar(hdr,'OBJECT'),/remove_all)
   if i GT 0 then allindices=[allindices,indices] else allindices=indices
endfor

myreadcol,'mindex.out',name,cn1,cn2,ca4227,g4300,fe4383,ca4455,fe4531,fe4668,hbeta,fe5015,mg1,mg2,mgb,fe5270,fe5335,fe5406,fe5709,fe5782,nad,tio1,tio2,hda,hga,hdf,hgf,format='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F'

iiref=[[cn1],[cn2],[ca4227],[g4300],[fe4383],[ca4455],[fe4531],[fe4668],[hbeta],[fe5015],[mg1],[mg2],[mgb],[fe5270],[fe5335],[fe5406],[fe5709],[fe5782],[nad],[tio1],[tio2],[hda],[hga],[hdf],[hgf]]

choose=[findgen(6),7]
allindices=allindices[choose]

end
