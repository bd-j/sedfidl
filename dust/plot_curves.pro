wave=findgen(3000)*3.+1E3

ca=calzetti_ztd(wave)
po7=power_ztd(wave,[0.7])
po13=power_ztd(wave,[1.3])
co=conroy_ztd(wave)


plot,wave,po13,yrange=[0.1,10],/ylog
oplot,wave,po7,color=255
oplot,wave,ca,color=djs_icolor('blue')
oplot,wave,co,color=djs_icolor('green')

end
