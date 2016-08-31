;; Clump radius vs. a plot

a = logarrw(min=0.05,max=20,n=1e3)
tcl = [2,5,12,30]
ncl = N_ELEMENTS(tcl)

distinct_colors
set_plot,'ps'
!P.FONT=0
csz = 1.8
thk = 8
device,file='rcl-a.eps',/enc,/col,xsize=8,ysize=6,/inche,/times

xr = minmax(a)
yr = [0.5,300]
plot,[0],[0],xrange=xr,yrange=yr,/xl,/yl,xtitle='Semi-major axis (au)',ytitle='Clump radius (R!dSun!n)',xthick=thk,ythick=thk,charsize=csz,position=[0.13,0.135,0.99,0.99],/yst,/xst;,ytickformat='exponentsonaxis';,xtickformat='exponentsonaxis'

;; RZ Psc
dr = [0.3,0.5]
polyfill,[dr[0],dr[1],dr[1],dr[0]],[yr[0],yr[0],yr[1],yr[1]],/line_fill,orientation=45,color=6,thick=thk
xyouts,0.43,0.6,'RZ Psc',charsize=csz,orientation=90

;; HAe/Be
dr = [0.6,1.2]
polyfill,[dr[0],dr[1],dr[1],dr[0]],[yr[0],yr[0],yr[1],yr[1]],/line_fill,orientation=45,color=2,thick=thk
xyouts,1.07,0.6,'UXOrs',charsize=csz,orientation=90

;; something to do with orbit circumference, 213 is a/RSun and 1rad around orbit
;oplot,a,a*213.7,thick=thk,color=2 ; clump goes 1 radian
;xyouts,0.65,150,'1rad',orientation=34,color=2,charsize=1.4

oplot,a,a*21.37,thick=thk,color=8 ; clump goes a/10
xyouts,4,93,'a/10',orientation=34.5,color=8,charsize=1.4

res = float(10)                              ; 11:10 inside and 10:9 outside
rin = 213*a*((res/(res-1))^(2./3.)-1)       ; how far out a 10:9 resonant orbit is from a
rout = 213*a*(1-(res/(res+1))^(2./3.))      ; ditto for 10:11
;oplot,a,0.5*(rin+rout)/2.,thick=thk,color=5 ; halve since will be sheared in both directions
;xyouts,1.35,65,'sheared in 2 orbits',charsize=1.4,color=5,orientation=34.5
;oplot,a,,thick=thk,color=5

col = [1,4,7,9]
for i=0,ncl-1 do begin

   rcl = 1.85*tcl[i]*sqrt(1./a) - 1.
   oplot,a,rcl,thick=thk,color=col[i]
   xyouts,0.1,1.1*interpol(rcl,a,0.1),STRN(tcl[i])+' days',color=col[i],charsize=1.4,orientation=-21
   
endfor

device,/close
set_plot,'x'

end
