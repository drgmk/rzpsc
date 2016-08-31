;; test beta particle model

q = 10.6/6.
a = 0.5
stot = 1.

nx = 1e5
x=linarrw(min=-1,max=100,n=nx)
dx = x[1]-x[0]

nt = 10
t=logarrw(min=1,max=10,n=nt)

!P.MULTI=[0,1,2]

tau = FLTARR(nx,nt)

;; tau
plot,[0],[0],yrange=[1e-2,40],xrange=MINMAX(x),/xst,/yl
for i=0l,nt-1 do begin
   tau[*,i] = taufunc_b(x,q,a,t[i],stot)
   oplot,x,tau[*,i]
endfor

;; total surface area
stot = FLTARR(nt)
for i=0,nt-1 do stot[i] = TOTAL(tau[*,i]*dx)
plot,t,stot,psym=1,/yst
oplot,t,stot

;end

!P.MULTI=0

nxc = 1000
xc = linarrw(min=-2,max=100,n=nxc)
fs = FLTARR(nxc,nt)
for i=0,nxc-1 do begin
   for j=0,nt-1 do begin
      fs[i,j] = fstarvis(5.5*(x-xc[i]),tau[*,j])
;   wait,0.1
   endfor
endfor

distinct_colors,n_colors=12
plot,[0],[0],yrange=[0.2,1],/xst,xrange=[-1,5];MINMAX(xc)
for i=0,nt-1 do oplot,xc,fs[*,i],color=i,thick=5

end

readcol,'all-lc.txt',t,f,o,format='d,f,a'
plot,t,f,psym=1,xrange=[2454060,2454080]
oplot,xc+2454065.3d,fs[*,3]

end
