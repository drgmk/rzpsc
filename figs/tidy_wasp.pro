;; detrend WASP lc

;; knobs
cotemp = 0.25                   ; min fraction of obs co-temporal with RZ Psc to keep lc
stdlim = 0.04                   ; max fractional dispersion allowed to keep lc (after first correction)
sigdecorlim = 0.1               ; max dispersion in decorrelation vector to keep data point
isplit = [0,3179]
isplit = [3180,4690]

dir = '../data/wasp-rz-psc-50closest/'
dir = '../data/wasp-rz-psc-5-13mag-50closest/'

f1 = '1SWASP J010942.05+275701.9.fits' ; RZ Psc itself

;; get lc for RZ Psc
fits_read,f1,tmp,pri,exten_no=0
jdref = SXPAR(pri,'JD_REF','oh noes!')
fits_read,f1,data,head,exten_no=1
t1 = tbget(head,data,'TMID')/86400.0+jdref
f1 = tbget(head,data,'FLUX2')
tam1 = tbget(head,data,'TAMFLUX2')
e1 = tbget(head,data,'FLUX2_ERR')
t1 = t1[isplit[0]:isplit[1]]
f1 = f1[isplit[0]:isplit[1]]
nt = N_ELEMENTS(t1)

;; get others
fn = FILE_SEARCH(dir+'*fits',COUNT=nf)
ts = DBLARR(nf,nt)
fs = DBLARR(nf,nt)
es = DBLARR(nf,nt)
med = FLTARR(nf)
sigref = FLTARR(nf)
ok = INTARR(nf)

for i=0,nf-1 do begin

   fits_read,fn[i],tmp,pri,exten_no=0
   jdref = SXPAR(pri,'JD_REF','oh noes!')
   fits_read,fn[i],data,head,exten_no=1
   t = tbget(head,data,'TMID')/86400.0+jdref
   f = tbget(head,data,'FLUX2')
   e = tbget(head,data,'FLUX2_ERR')
   for j=0,nt-1 do begin
      tmp = WHERE(t eq t1[j],ntmp)
      if ntmp eq 1 then begin
         ts[i,j] = t[tmp]
         fs[i,j] = f[tmp]
         es[i,j] = e[tmp]
      endif
   endfor

   tmp = WHERE(fs[i,*] ne 0,nok) ; ignore curves that are all zero
   if nok eq 0 then continue

;   print,nok/FLOAT(nt)
   med[i] = MEDIAN(fs[i,tmp])           ; light curve median
   sigref[i] = STDDEV(fs[i,tmp])/med[i] ; light curve fractional dispersion

   if nok/FLOAT(nt) lt cotemp then continue
   ok[i] = 1

   plot,fs[i,tmp]/med[i],psym=3,yrange=[0.5,2]
;   wait,0.01

   counter,i,nf,'done'
endfor   

!P.MULTI=[0,2,5]
plot,sigref
histogauss,sigref,a

;; get initial normalisation vector
fdecor = FLTARR(nt)
sigdecor = FLTARR(nt)
for i=0,nt-1 do begin
   qt = WHERE(fs[*,i] ne 0 and ok eq 1,nq)
   fdecor[i] = TOTAL( fs[qt,i]/med[qt] ) / FLOAT(nq) ; iniital decorrelation vector
   sigdecor[i] = STDDEV( fs[qt,i]/med[qt] )          ; dispersion in vector
endfor

plot,fdecor,psym=3,/yst
plot,sigdecor

;; create initial corrected light curves
fscor = FLTARR(nf,nt)
fssig = FLTARR(nf)
for i=0,nf-1 do begin
   fscor[i,*] = fs[i,*]/fdecor                                 ; corrected fluxes
   tmp = WHERE(fscor[i,*] ne 0)                                ;
   if tmp[0] eq -1 then continue                               ;
   fssig[i] = ROBUST_SIGMA(fscor[i,tmp]/MEDIAN(fscor[i,tmp])) ; fractional dispersion of corrected
endfor
fs = fscor                      ; update light curves
f1 /= fdecor                    ;

plot,fssig
histogauss,fssig,a

;; now ignore curves that have a large dispersion remaining
tmp = WHERE(fssig gt stdlim)
ok[tmp] = 0

;; now get updated decorrelation curve
fdecor = FLTARR(nt)
sigdecor = FLTARR(nt)
nused = INTARR(nt)
for i=0,nt-1 do begin
   qt = WHERE(fs[*,i] ne 0 and ok eq 1,nq)
   nused[i] = nq
   fdecor[i] = TOTAL( fs[qt,i]/med[qt] ) / FLOAT(nq) ; decorrelation vector
   sigdecor[i] = STDDEV( fs[qt,i]/med[qt] )          ; dispersion in vector
endfor

plot,fdecor,psym=3,/yst
plot,sigdecor

;; set badly constrained points to 0 -> NaN
tmp = WHERE(sigdecor gt 0.1)
fdecor[tmp] = 0.

;; create final corrected light curves
fscor = FLTARR(nf,nt)
medcor = FLTARR(nf)
sigcor = FLTARR(nf)
for i=0,nf-1 do begin
   fscor[i,*] = fs[i,*]/fdecor                                 ; corrected fluxes
   tmp = WHERE(fscor[i,*] ne 0)                                ;
   if tmp[0] eq -1 then continue                               ;
   medcor[i] = MEDIAN(fscor[i,tmp])
   sigcor[i] = ROBUST_SIGMA(fscor[i,tmp]/medcor[i]) ; fractional dispersion of corrected
endfor
f1 /= fdecor

plot,sigcor
histogauss,sigcor,a
!P.MULTI=0

;end

;; plot,t1,f1cor/MEDIAN(f1cor),yrange=[-1,27],/xst,/yst,psym=3;,xrange=[2453960,2454105];[2453160,2453270]
;; for i=0,nf-1 do begin
;;    if ok[i] eq 0 then continue
;;    oplot,t1,fscor[i,*]/med[i]+1+i/2.,psym=3
;; endfor
;; oplot,t1,sigdecor*2-1,psym=6

xr = [0.9,1.1]
loadct,0
nbins = 100
yall = FLTARR(nbins)
delvarx,histtot
plot,[0],[0],xrange=xr,yrange=[0,65],/xst,/yst
sunset_colors
for i=0,nf-1 do begin
   if ok[i] eq 0 then continue
   tmp = WHERE( fscor[i,*] ne 0 and FINITE(fscor[i,*]) )
   histxy,hist=fscor[i,tmp]/medcor[i],x=x,y=y,min=xr[0],max=xr[1],nbins=nbins
   push,histtot,fscor[i,tmp]/medcor[i]
   oplot,x,SMOOTH(y,2),psym=10,color=i*5
   for j=0,nbins-1 do yall[j] += y[j]
endfor
histogauss,histtot,a,/NOPLOT
xs = linarrw(min=xr[0],max=xr[1],n=nbins)
oplot,xs,a[0]*exp( -( (xs-a[1])/(sqrt(2)*a[2]) )^2 )/nbins/(xr[1]-xr[0])

;; cleanup and save
tmp = WHERE(FINITE(f1))
t1 = t1[tmp]
f1 = f1[tmp]
save,file='wasp1-lc.xdr',f1,t1

end
