;; montage of interesting light curves

!P.MULTI=[0,1,2]
!P.FONT=0
thk = 6
ssz = 0.35
csz = 1.5
distinct_colors

yl = [0.22,0.55,0.99]
yl[1] = (yl[2]-yl[0])/2. + yl[0]

set_plot,'ps'
device,file='lcs.eps',xsize=8,ysize=2,/inch,/time,/col,/enca

xr = 99                         ; days
yr = [0.71,1.04]
yn = [' ','0.8',' ','0.9',' ','1.0']

;; KIC 8462852
x0 = 1490
readcol,'~/astro/collab/tabethaboyajian/kic8462852/data/KIC8462852_normalized_flux.txt',t,f,format='d,d'
plot,[0],[0],yrange=yr,xrange=[0,xr],/xst,/yst,position=[0.07,yl[1],0.99,yl[2]],ytitle='brightness  ',xtickname=REPLICATE(' ',60),xthick=thk,ythick=thk,yticklen=0.01,xticklen=0.08,ytickname=yn
oplotcircw,t-x0,f,/fill,symsize=ssz,col=1
legend,['KIC 8462852','(space)'],/bot,/lef,box=0,margin=0

;; RZ Psc
x0 = 2453170.
readcol,'~/astro/doc/rz-psc/figs/all-lc.txt',t,f,format='d,d'
plot,[0],[0],yrange=yr,xrange=[0,xr],/xst,/yst,position=[0.07,yl[0],0.99,yl[1]],ytitle='  Relative',xtitle='Days',xthick=thk,ythick=thk,yticklen=0.01,xticklen=0.08,ytickname=yn
oplotcircw,t-x0,f,/fill,symsize=ssz,col=1
legend,['RZ Psc','(ground)'],/bot,/lef,box=0,margin=0

;; EPIC 203937317
;; readcol,'~/astro/collab/meganansdell/k2-dippers/data/k2/203937317.txt',t,f,format='d,d'
;; plot,[0],[0],yrange=yr,xrange=2092+[0,xr],/xst,/yst
;; oplot,t,f,psym=6,symsize=ssz
;; legend,['EPIC 203937317'],/bot,/rig,box=0

;; EPIC 205519771
;; readcol,'~/astro/collab/meganansdell/k2-dippers/data/k2/205519771.txt',t,f,format='d,d'
;; plot,[0],[0],yrange=yr,xrange=2092+[0,xr],/xst,/yst
;; oplot,t,f,psym=6,symsize=ssz
;; legend,['EPIC 205519771'],/bot,/rig,box=0

device,/close


device,file='lc-cadence.eps',xsize=8,ysize=3,/inch,/time,/col,/enca

ssz = 0.5
xr = 12                         ; days
yr = [0.0,1.1]
yn = [' ','0.2',' ','0.6',' ','1.0']

x0 = 2453987.d
readcol,'~/astro/doc/rz-psc/figs/all-lc.txt',t,f,format='d,d'

plot,[0],[0],yrange=yr,xrange=[0,xr],/xst,/yst,position=[0.07,yl[1],0.99,yl[2]],ytitle='brightness  ',xtickname=REPLICATE(' ',60),xthick=thk,ythick=thk,yticklen=0.01,xticklen=0.08,ytickname=yn
oplotcircw,t-x0,f,/fill,symsize=ssz,col=1
legend,['x0='+STRN(x0,format='(I)')],/bot,/lef,box=0

;; RZ Psc
plot,[0],[0],yrange=yr,xrange=[0,xr],/xst,/yst,position=[0.07,yl[0],0.99,yl[1]],ytitle='  Relative',xtitle='Days - x0 (two deep WASP events in 2006 season)',xthick=thk,ythick=thk,yticklen=0.01,xticklen=0.08,ytickname=yn

;; actual
x0 = 2454062.d
oplotcircw,t-x0,f,/fill,symsize=ssz,col=1
legend,['x0='+STRN(x0,format='(I)')],/bot,/lef,box=0

;; simulated
;; readcol,'~/astro/doc/rz-psc/figs/fitdata.txt',j,a,b,c,d,e,FORMAT='i,d,f,f,f,f',count=nl
;; delvarx,p
;; for j=0,nl-1 do if a[j] gt x0 and a[j] lt x0+xr then push,p,[a[j],b[j],c[j],d[j],e[j]]
;; x = linarrw(min=0,max=xr,n=xr*4)

;; tmp = WHERE(x gt 10.3 and x lt 18,n)
;; x1 = linarrw(min=x[tmp[0]],max=x[tmp[-1]],n=n*3)

;; oplotcircw,x,dipwrap(x0+x,p),col=5,thick=4,/fill,symsize=ssz ; 6h cadence
;; oplotcircw,x1,dipwrap(x0+x1,p),col=8,thick=2,/fill,symsize=0.3 ; 2h cadence

legend,['RZ Psc'],/bot,/rig,box=0,margin=0

;; EPIC 203937317
;; readcol,'~/astro/collab/meganansdell/k2-dippers/data/k2/203937317.txt',t,f,format='d,d'
;; plot,[0],[0],yrange=yr,xrange=2092+[0,xr],/xst,/yst
;; oplot,t,f,psym=6,symsize=ssz
;; legend,['EPIC 203937317'],/bot,/rig,box=0

;; EPIC 205519771
;; readcol,'~/astro/collab/meganansdell/k2-dippers/data/k2/205519771.txt',t,f,format='d,d'
;; plot,[0],[0],yrange=yr,xrange=2092+[0,xr],/xst,/yst
;; oplot,t,f,psym=6,symsize=ssz
;; legend,['EPIC 205519771'],/bot,/rig,box=0

device,/close
set_plot,'x'
!P.MULTI=0


end
