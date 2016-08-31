;; look for reflectance spectrum signature

bs = ['UJ','BJ','VJ','RC','IC']
v = 13.67

;; dip colours
b = 0.81+v
u = 0.14 + b
r = v - 0.72
i = v - 1.33
dipms = [u,b,v,r,i]

;; normal colours
b = 0.7 + v
u = 0.2 + b
r = v - 0.7
i = v - 1.1
starms = [u,b,v,r,i]

fdip = FLTARR(5)
fstar = FLTARR(5)
for j=0,4 do begin
   fdip[j] = mag2jy(dipms[j],bs[j])
   fstar[j] = mag2jy(starms[j],bs[j])
endfor

;; restore,'~/astro/debrisdisks/rz-psc/RZ-Psc.xdr'
;; fstar = FLTARR(5)
;; for j=0,4 do begin
;;    tmp = WHERE(starsynth.band eq bs[j])
;;    fstar[j] = starsynth.fnu[tmp]
;; endfor

wavs = getweffs(bs)

plot,wavs,fstar
oplot,wavs,fdip

plot,wavs,fdip/fstar

end
