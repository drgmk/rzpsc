;; model of something occulting a star

;; fraction of total stellar brightness given u
function fstar,u
  
  nu = N_ELEMENTS(u)
  fstar = FLTARR(nu)
  
  tmp = WHERE(ABS(u) gt 1,COMPLEMENT=in,NCOMPLEMENT=nin)
  bef = WHERE(u lt -1,nbef)
  aft = WHERE(u gt 1,naft)
  if nbef gt 0 then fstar[bef] = FLTARR(nbef)+1.0 ; if before transit then no dimming
  if naft gt 0 then fstar[aft] = FLTARR(naft)     ; if after transit then dimmed 0

  if nin gt 0 then begin
     fstar[in] = 1/!pi*(acos(u[in])-u[in]*sqrt(1-u[in]^2)) ; eq 3 of Winn+06

     ;; differentiated w.r.t. u (dudr=1/rstar)
     dfdu = FLTARR(nu)
     dfdu[in] = -2*!pi*sqrt(1-u[in]^2)
  endif
  
  return,fstar
end

function occult,t,tc,du,tau,dudt

  ;;convert time to units of stellar radii
  u = dudt * ( t - tc ) + du/2.
  nu = N_ELEMENTS(u)
  f = FLTARR(nu)
  fin = fstar(u)
  fout = fstar(-u+du)
  for i=0,nu-1 do begin
     fvis = fin[i]+fout[i]
     f[i] = fvis + (1-fvis)*(1-tau)
  endfor
  return,f
end

;; uniform optically thick clump that is larger than star
nu = 1e3
u = linarrw(min=-2,max=4,n=nu) ; distance from star center along distance of travel in u
tau = 0.5                      ; uniform optical depth of cloud
du = 2                        ; width of cloud in u (2=diameter of star)

;; plot,u,fstar(u)
;; oplot,u,fstar(-u+du),linestyle=1
;; oplot,u,occult(u,0,du,tau,3),thick=4

!P.MULTI=[0,1,2]

readcol,'all-lc.txt',t,f,o,format='d,f,a'
plot,t-2453900,f,psym=6,xrange=[90,96],symsize=0.1,thick=2
u = linarrw(min=90,max=180,n=nu)
oplot,u,occult(u,90.64,2,0.58,3.5)
oplot,u,occult(u,92,2,0.45,1.5)
oplot,u,occult(u,94,2,0.2,1.5)
oplot,u,occult(u,94.6,2,0.88,5.5)

plot,t-2453900,f,psym=6,xrange=[160,180],symsize=0.2,thick=2
oplot,u,occult(u,169.2,2,0.58,3.5)
oplot,u,occult(u,164.85,2,0.45,1.5)
oplot,u,occult(u,174,2,0.2,1.5)
oplot,u,occult(u,165.4,2,0.88,5.5)

!P.MULTI=0

end
