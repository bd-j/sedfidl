;+
; NAME:
;  RECONSTRUCT_SFH
;
; VERSION:
;   1.0 (Apr, 2010)
;
; PURPOSE:
;   to add together the fluxes and other quantities from structures
;   containg information about hosts models and burst models
;
;
; REFERENCE:
;   Bruzual & Charlot, 2003, MNRAS
;
; CALLING SEQUENCE:
;   
;
;
; CATEGORY:
;   Population Synthesis
;
; INPUTS:
;   hosts - structure containing the output from photo_grid_Ztd
;   bursts - structure containing the output from burst_grid_Ztd
;   
; OPTIONAL INPUTS:

FUNCTION DelayedBurst,tau,age,times=times,dt=dt
;sfr ~ tau^(-2)te^(-t/tau)
  dd=dt
  aa=double(age/tau)
  ndt=n_elements(dt)
  bar_sfr=fltarr(n_elements(aa),ndt)
  for i=0,ndt-1 do begin ;$ ;vectorize?
     dt=double(dd[i]/tau)
     bar_sfr[*,i]=(igamma(2,aa,/double)-igamma(2,(aa-dt)>0,/double))/(dt*tau)
  endfor
  dt=dd
  return,bar_sfr
end

FUNCTION BoxcarBurst,tau,age,times=times,dt=dt
;sfr(t)=1/tau for t<tau, sfr(t)=0 for t>tau
;integral is the sum of a bunch of of ramp functions, which are the
;integrals of offset Heaviside functions with opposite sign.  Joy.
  ndt=n_elements(dt)
  bar_sfr=fltarr(n_elements(age),ndt)
  for i=0,ndt-1 do begin ;vectorize?
     int_sfr=age*(age GT 0) - (age-tau)*(age-tau GT 0) - $
             (age-dt[i])*(age-dt[i] GT 0) + (age-Dt[i]-tau)*(age-Dt[i]-tau GT 0)
     bar_sfr[*,i]=int_sfr/dt[i]/tau
  endfor
  return,bar_sfr
end

FUNCTION ExponentialDecay,tau,age,dt=dt,times=times
  ;sfr~1/tau*exp(-t/tau)
  dd=dt
  aa=double(age/tau)
  ndt=n_elements(dd)
  bar_sfr=fltarr(n_elements(aa),ndt)
  for i=0,ndt-1 do  begin ;$
     dt=double(dd[i]/tau)
     bar_sfr[*,i]=(exp(0.-((aa-dt)*((age-dd[i]) GT 0)))-exp(0.-aa))/(dt*tau)*(2.*(tau GT 0)-1.)
  endfor

  dt=dd

  return,bar_sfr
end


FUNCTION Reconstruct_SFH,grid,times=times,dt=dt,scale=scale
  if keyword_set(scale) EQ 0 then scale=1
  ;timescales over which to get the average sfr:
  if n_elements(dt) EQ 0 then dt=[1E7,1E8,3E8] 
  ndt=n_elements(dt)
  nmod=n_elements(grid)
  nburst=n_elements(grid[0].burst_weight)
  
  ;parse the burst shape
  delayed=(strpos(strlowcase(grid.burst_func),'delay') GE 0 $
           and grid.burst_age GT 0)
  boxcar=(strpos(strlowcase(grid.burst_func),'delay') LT 0 $
          and grid.burst_age GT 0)
  ndelay=total(delayed)
  nbox=total(boxcar)
;find a number within the string, call that tau
  tau=float(stregex(grid.burst_func,'[0123456789]+Myr',/extract))*1E6

  bsfr=fltarr(nburst,nmod,ndt)
  junk=bsfr
  bsfr=reform(bsfr,[nburst*nmod,ndt],/overwrite)
  junk[*,*,0]=grid.burst_age
  junk=reform(junk,[nburst*nmod,ndt],/overwrite)
  if ndelay GT 0 then begin
     sub=where(delayed)
     
     bsfr[sub,*]=DelayedBurst(tau[sub],$
                       (grid.burst_age)[sub],dt=dt)
  endif

  if nbox GT 0 then begin
     sub=where(boxcar)
     bsfr[sub,*]=BoxcarBurst(tau[sub],$
                             (grid.burst_age)[sub],dt=dt)
  endif

  bsfr=reform(bsfr,[nburst,nmod,ndt],/overwrite)*$
       rebin(reform(grid.burst_weight,[nburst,nmod,1]),[nburst,nmod,ndt])
  bsfr=total(bsfr,1)


  hsfr=ExponentialDecay(grid.tau_sf*1E9,grid.age_host,dt=dt);*host_norm

  sfr=hsfr+scale*bsfr


return,sfr


end
