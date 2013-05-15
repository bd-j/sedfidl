FUNCTION wavg,val,err,nsig=nsig,niter=niter,dev=dev,good=index

f=where(finite(val) and finite(1/err))

np=n_elements(val)
ivar=1./err^2.
index=findgen(n_elements(val))

newval=val[f]
newivar=ivar[f]
index=index[f]

;check these formulas
aver=total(newval*newivar)/total(newivar)
sigma=1./sqrt(total(newivar))
sigma2=sqrt(total((val-aver)^2.*ivar)/total(ivar))

if keyword_set(nsig) then begin
  if keyword_set(niter) EQ 0 then niter=1
  for iter=0,niter-1 do begin
    good = where (abs(newval-aver) LT nsig*sigma)
    index=index[good]
    newval=newval[good]
    newivar=newivar[good]
    np=n_elements(good)
    aver=total(newval*newivar)/total(newivar)
    sigma=stddev(newval-aver)
    sigma=1./sqrt(total(newivar))
  endfor
endif

dev=[sigma,sigma2]

return,aver



end
