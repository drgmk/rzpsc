;; which clumps to fit
yr = 2008
jdcnv,yr,5,1,0,t0         ; 1st May of each year
day = 150
t0 += day
print,t0

;; get clumps
readcol,'fitdata.txt',j,a,b,c,d,e,FORMAT='i,d,f,f,f,f',count=nl
delvarx,p,nos
extra = 50
tspan = [t0-extra,t0+extra] ; 20 days either side

;; read data
readcol,'all-lc.txt',t,f,o,format='d,f,a'
tmp = WHERE(t gt tspan[0] and t lt tspan[1])
t = t[tmp]
f = f[tmp]
o = o[tmp]

;; get params for clumps in this time range
for i=0,nl-1 do begin
   if a[i] gt MIN(t) and a[i] lt MAX(t) then begin
      push,p,[a[i],b[i],c[i],d[i],e[i]]
      push,nos,j[i]
   endif
endfor
ncl = N_ELEMENTS(p)/5
in = INDGEN(ncl)*5              ; starting index of each clump in p

plot,t,f,psym=1,symsize=2.0,thick=1;,xrange=[2454060,2454080]
t1 = linarrw(min=MIN(t),max=MAX(t),n=4e3)

oplot,t1,dipwrap(t1,p),linestyle=2
for i=0,ncl-1 do print,nos[i],p[INDGEN(5)+in[i]]
;end

fixcl = INTARR(ncl)+0           ; fix all by default
;; unfix = clfit                   ; decide which clumps to fit
;; for i=0,ncl-1 do begin
;;    tmp = WHERE(nos[i] eq unfix,ntmp)
;;    if ntmp gt 0 then fixcl[i] = 0
;; endfor
pi = replicate({fixed:1, relstep:0.1, limited:[1,0], limits:[0.D,0]},N_ELEMENTS(p))
pi[in+1].limits[0:1] = 0.1                                         ; min edge sharpness
pi[in+2].limits[0:1] = 0.1                                         ; min edge sharpness
pi[in].relstep = 1e-7                                              ; smaller step for dates
for i=0,ncl-1 do pi[INDGEN(5)+5*i].fixed = fixcl[i]                ; free up some clumps
pi[in+4].fixed = 1                                                 ; re-ensure we fix velocities
res = mpfitfun('dipwrap',t,f,0.02,p,parinfo=pi,yfit=fit,status=st);,xtol=1e-5)
f1 = dipwrap(t1,res)
oplot,t1,f1

;; sort and output params
srt = SORT(res[in])
for i=0,ncl-1 do print,STRING([nos[i],res[INDGEN(5)+5*srt[i]]],FORMAT='(i3," ",d12.3," ",f5.2," ",f5.2," ",f5.2," ",f5.2)')

end

;; re write out params
readcol,'fitdata.txt',j,a,b,c,d,e,FORMAT='i,d,f,f,f,f',count=nl
delvarx,p
for i=0,N_ELEMENTS(a)-1 do push,p,[a[i],b[i],c[i],d[i],e[i]]
srt = SORT(a)
for i=0,nl-1 do print,STRING([i,a[srt[i]],b[srt[i]],c[srt[i]],d[srt[i]],e[srt[i]]],FORMAT='(i3," ",d12.3," ",f5.2," ",f5.2," ",f5.2," ",f5.2)')

end
