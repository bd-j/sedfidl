FUNCTION draw_params_Ztd,nmod,tausf=tausf,tauv=tauv,mu=mu,met=met,grid_sfh=grid_sfh,grid_tauv=grid_tauv,grid_mu=grid_mu,redshifts=redshifts,max_age=max_age


maxtausf=12  ;Gyr
mintausf=0.
;the desired, existing bc03 grid parameters
if n_elements(grid_sfh) EQ 0 then grid_sfh=[(findgen(10)+1)*0.2,(findgen(10)+5)*0.5,findgen(5)+8,-13,-10,-7,-4,-2,-1,-0.5]
metgrid=[0.02,0.2,0.4,1.0,2.5]

;generate random numbers
x=randomu(seed,nmod) ;for gamma_SF
y=randomu(seed,nmod) ;for tau_v
z=randomn(seed,nmod) ;for mu
k=randomu(seed,nmod) ;for metals
h=randomu(seed,nmod) ;for redshifts
l=randomu(seed,nmod) ;for ages

;;----gamma.----
ns=n_elements(grid_sfh)
srange=minmax(grid_sfh)
ns=srange[1]-srange[0]
sfh_ptr=fix(x*ns)+srange[0]
tausf=sfh_ptr
;tausf=sfh_grid[nearest(sfh_grid,sfh_ptr)]

;;-----tau_v----
;get something to match LVL based on 24/Ha
line=(findgen(1350)+1)*0.003

;c=line-(alog(cosh(1.5*line-6.7)))/1.5+alog(cosh(-6.7))/1.5 ;da Cunha, p=1-tanh(1.5*tv-6.7)
;c=1.05*line-(alog(cosh(4.5*line-4.0)))/4.5+alog(cosh(-4.0))/4.5 ;match lvl (v1)

c=1.1*line-(alog(cosh(5.5*line-3.0)))/5.5+alog(cosh(-3.0))/5.5 ;match lvl (v2)
c=c/max(c)
j=nearest(c,y)
tauv=line[j]


;;----- mu -----
;want gaussian around 0.3 with 0.1 scatter

bar=0.3
sigma=0.1  ;(v2) (v1 value was 0.2)
z=(z*sigma)+bar
bad=where(z LT 0.1 or z GT 1,nbad)
while nbad GT 0 do begin
   zz=randomn(seed,nbad)
   z[bad]=(zz*sigma)+bar
   bad=where(z LT 0.1 or z GT 1,nbad)
endwhile
mu=z

;;--- metallicity ----
;uniform in Z?
met=k*max(metgrid)
;no, weight towards lower z to match LVL gas phase better, max is 2.5 solar
line=(findgen(830)+1)*0.003
c=1.05*line-(alog(cosh(5.0*line-3.0)))/5.0+alog(cosh(-3.0))/5.0 ;match lvl (v1,v2)
c=c/max(c)
j=nearest(c,k)
met=line[j]

;;--- redshift ---
;;weighted by volume and distance modulus and mass
;;function or luminosity function, or taken from spectroscopic redshifts
;;p(z)=(z/z0)*exp(0.-(1/(z0/z+1/20.)))  need to integrate this
;;p(z)=(z/z0)*exp(0.-(z/z0))  integral is the incomplete gamma function
;z0=1.  
;line=(findgen(fix(700/z0))+1)*0.01
;c=z0*igamma(2,line);+0.05*line
;c=c/max(c)
;j=nearest(c,h)
;z=line(j)*z0

max_z=0.
z=h*max_z ;flat

;;--- age ---
if n_elements(max_age) EQ 0 then max_age=13.7E9 ;~age of the universe
if keyword_set(redshifts) then begin 
   H_0=72.*3.241E-20*3.156E16   ;in Gyr^{-1}
   grid_age_universe=13.7-lookback(grid_redshift,0.27,0.73)/H_0
   max_age=grid_age_universe
endif
age=l*max_age ;flat prior on ages (up to the age of the universe at that redshift)

;;----snap to grid------
;take random parameter values and re-bin them onto the available model grids
f=nearest(grid_sfh,tausf)
tausf_reg=grid_sfh[f]
tauv_reg=tauv
mu_reg=mu
met_reg=met>min(metgrid)<max(metgrid)
age_reg=age
z_reg=z

return,[[tausf_reg],[tauv_reg],[mu_reg],[met_reg],[age_reg],[z_reg]]


end
