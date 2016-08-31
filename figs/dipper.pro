;; model of something occulting a star

;; function for cloud optical depth as two 1D gaussians, distance units are stellar radii
function taufunc_g,x,hwhm1,hwhm2,taumax
  nx = N_ELEMENTS(x)
;  if taumax gt 1 then print,'tau (WARNING): asked for depth greater than 1'
  lo = WHERE(x le 0,nlo,COMPLEMENT=hi,NCOMPLEMENT=nhi)
  g = FLTARR(N_ELEMENTS(x))
  if nlo gt 0 then g[lo] = GAUSSIAN(x[lo],[taumax,0,2*hwhm1/2.35])
  if nhi gt 0 then begin        ; for large hwhm2 just make flat
     if hwhm2 ge 100 then g[hi] = FLTARR(nhi)+1.0 else g[hi] = GAUSSIAN(x[hi],[taumax,0,2*hwhm2/2.35])
  endif
  return,g
end

;; function for cloud optical depth with a tail
function taufunc_b,x,q,a,t,stot

  if t eq 0 then return,FLTARR(N_ELEMENTS(t)) ; nothing has happened yet
  
  dx = x[1]-x[0]                 ; bin sizes
  c = 0.5                        ; constant to convert beta to diameter
  dminbeta = c / 0.5             ; blowout size
  mstar = 1.                     ; stellar mass
  n = 2*!pi*sqrt(mstar)/a^(3/2.) ; mean motion, assuming some stellar mass
  maxDM = n * t                  ; max DM for bound particles (i.e. at dminbeta)
  maxx = maxDM * a               ; max distance behind clump particles can get to now
  Df = x / FLOAT(a)              ; distance behind clump
  deltaf = dx / FLOAT(a)         ; width of small region behind clump
  dmax = c*n * t / Df
  dmin = c*n * t / (Df+deltaf)
  fmtq = 5. - 3 * q

  ;; normalisation, K is the whole lot in front
  if q gt 5/3. then K = -1. * stot / dminbeta^fmtq
 
  sig = K * ( dmax^fmtq - dmin^fmtq )

  tmp = WHERE(x gt maxx)              ; no particles beyond where beta=0.5
  if tmp[0] ne -1 then sig[tmp] = 0.0 ;
  tmp = WHERE(x lt 0)                 ; or in front of the clump
  if tmp[0] ne -1 then sig[tmp] = 0.0 ;
    
  return,sig/dx                 ; 1-D optical depth
end

;; fraction of star visible
function fstarvis,x,taux,test=TEST

  ok = WHERE(x gt -1 and x lt 1)     ; only integrate over stellar surface
  xint = [-1,x[ok],1]                ; x is in units of stellar radii, add -1 and 1 points
  taum1 = INTERPOL(taux,x,-1)        ;
  taup1 = INTERPOL(taux,x,1)         ;
  tau = [taum1,taux[ok],taup1] < 1.0 ; optical depth no greater than 1 in integration
  
  ;; uniform case, better to compute amount blocked
  int = tau * 2 * sqrt(1-xint^2)
  fvis = 1.0 - INT_TABULATED( xint, int, /DOUBLE ) / !pi

  if KEYWORD_SET(test) then begin
     !P.MULTI=[0,2,1]
     plot,x,taux,xrange=MINMAX(x),yrange=[0,1]
     oplot,[0,0],MINMAX(taux),linestyle=2
     oplot,[1,1],MINMAX(taux),linestyle=1
     oplot,[-1,-1],MINMAX(taux),linestyle=1
     oplot,xint,tau,thick=4
     plot,xint,int,xrange=MINMAX(x),yrange=[0,2]
     !P.MULTI=0
;     wait,0.05
  endif

  return,fvis
end

;; do the dips
function dipper,t,tc,hwhm1,hwhm2,taumax,dxdt,nx=nx,test=test

  nt = N_ELEMENTS(t)                   ; number of times desired
  ncl = N_ELEMENTS(tc)                 ; number of clumps
  if N_ELEMENTS(hwhm1) ne ncl or N_ELEMENTS(hwhm2) ne ncl or $ ; sanity check
     N_ELEMENTS(taumax) ne ncl or N_ELEMENTS(dxdt) ne ncl then $
        STOP,'dipper: need same number of params for each input'
  ;; if N_ELEMENTS(nx) eq 0 then nx = 100 ; number of points across star
  ;; x = linarrw(min=-1,max=1,n=nx)       ; x array for cloud, x is in stellar radii
  ;;  xcl = FLTARR(ncl,nt)
  ;;  for i=0,ncl-1 do xcl[i,*] = dxdt[i]*(t-tc[i]) ; position at each time for each clump
  ;; fs = FLTARR(nt)
  ;; for i=0,nt-1 do begin         ; add up tau's geometrically, redo at each time step
  ;;    tau = FLTARR(nx)
  ;;    for j=0,ncl-1 do tau += taufunc_g(x+xcl[j,i],hwhm1[j],hwhm2[j],taumax[j])
  ;;    fs[i] = fstarvis(x,tau,test=test)
  ;; endfor
  ;; return,fs; > 0.15               ; minimum to account for scattered light

  ;; use C version
   str = STRARR(ncl*5)
   for i=0,ncl-1 do str[5*i:4+i*5] = [tc[i],hwhm1[i],hwhm2[i],taumax[i],dxdt[i]]
   spawn,'/Users/grant/astro/code/c/dippers/dip_gauss1d '+STRN(ncl)+' '+STRJOIN(str,' ')+' '+STRJOIN(t,' '),tmp
   tmp = STRSPLIT(tmp,' ',/EXTRACT)
   fs = FLOAT(tmp)
   return,fs            ; > 0.15

end
