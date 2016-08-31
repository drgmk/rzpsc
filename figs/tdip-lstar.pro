;; picture showing difference between dippers, UXOrs, and RZ Psc

distinct_colors
set_plot,'ps'
!P.FONT=0
csz = 1.8
thk = 8
device,file='tdip-lstar.eps',/enc,/col,xsize=8,ysize=6,/inche,/times

xr = [5e-3,50]
yr = [0.1,40]

plot,[0],[0],/xst,/yst,xrange=xr,/xl,yrange=yr,/yl,xtitle=textoidl('L_{star} (L_{Sun})'),ytitle=textoidl('t_{dim} (days)'),xthick=thk,ythick=thk,charsize=csz,position=[0.13,0.135,0.99,0.99],xtickformat='exponentsonaxis'

;; dippers
lr = [0.01,1.5]
tr = [0.5,2]
polyfill,[lr[0],lr[1],lr[1],lr[0]],[tr[0],tr[0],tr[1],tr[1]],/line_fill,orientation=45,color=6,thick=thk
xyouts,0.05,2,'dippers',charsize=csz

;; UXOrs
lr = [5,20]
tr = [14,30]
polyfill,[lr[0],lr[1],lr[1],lr[0]],[tr[0],tr[0],tr[1],tr[1]],/line_fill,orientation=45,color=2,thick=thk
xyouts,0.05,2,'UXOrs',charsize=csz



device,/close
set_plot,'x'

end
