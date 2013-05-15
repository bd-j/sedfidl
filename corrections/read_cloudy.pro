PRO read_cloudy,root,energy=energy,name=name,intrinsic=intrinsic,emergent=emergent,wave=wave,diffout=diffout,netout=netout,radius=radius,te=te,eden=eden,hih=hih,hiih=hiih,h2h=h2h,incident=incident


openr,15,root+'.lint'
header=''
readf,15,format='(A)',header

e=0.
n=''
i=0.
ii=0.

energy=0.
name=''
emergent=0.
intrinsic=0.
k=0.
while (~ EOF(15)) do begin
  readf,15,format='(F12,A14,F9,F9)',e,n,i,ii
  energy=[energy,e]
  name=[name,n]
  intrinsic=[intrinsic,i]
  emergent=[emergent,ii]
  k=k+1
endwhile
close,15

readcol,root+'.cont',ryd,incident,trans,diffout,netout,ref,tot,/silent

wave=911.6/ryd

readcol,root+'.hyd',radius,te,hden,eden,hih,hiih,h2h,/silent


end
