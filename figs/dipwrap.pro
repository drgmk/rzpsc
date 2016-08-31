;; wrapper to get light curves

function dipwrap,x,p
  if N_ELEMENTS(p) eq 0 then return,FLTARR(N_ELEMENTS(x))+1
  n = N_ELEMENTS(p)/5
  in = INDGEN(n)*5
  tc = p[in+0]
  hwhm1 = p[in+1]
  hwhm2 = p[in+2]
  taumax = p[in+3]
  dxdt = p[in+4]
  f = dipper(x,tc,hwhm1,hwhm2,taumax,dxdt,nx=50)
  return,f
end
