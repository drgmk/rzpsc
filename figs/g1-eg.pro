;; show some properties of the dust model

!P.FONT=0
csz = 1.5
thk1 = 8
thk = 12
distinct_colors,n_colors=12
set_plot,'ps'
device,file='g2-eg.eps',/enca,xsize=8,ysize=6,/inches,/times,/col
xr = [-6,3]

plot,[0],[0],xr=xr,yrange=[0,1],xtitle='time',ytitle='Relative flux',xthick=thk1,ythick=thk1,charsize=csz,/xst,/yst,position=[0.1,0.1,0.99,0.99]

;; how things change with parameters
x = linarrw(min=xr[0],max=xr[1],n=1e2)
cols = [1,4,9]

;; varying g1
g1 = [0.3,2,4];logarrw(min=0.1,max=2,n=4)
legend,['du/dt=1',textoidl('\sigma_{in}='+STRING(g1,format='(f3.1)'))],/bot,/lef,box=0,linestyle=[-1,2,2,2],charsize=csz,thick=thk,margin=0,colors=[0,cols]
for i=0,N_ELEMENTS(g1)-1 do begin
   f = dipper(x,0+g1[i],g1[i],100,1,1)
   oplot,x,f,linestyle=2,thick=thk,color=cols[i]
endfor

;; varying dxdt
dxdt = 2*[0.5,0.25,0.125]
legend,[textoidl('\sigma_{in}=1'),'du/dt='+STRING(dxdt,format='(F4.2)')],/lef,/center_legend,box=0,linestyle=[-1,0,0,0],charsize=csz,thick=thk,margin=0,colors=[0,cols]
for i=0,N_ELEMENTS(dxdt)-1 do begin
   f = dipper(x,0+0.01,0.01,100,1,dxdt[i])
   oplot,x,f,thick=thk,col=cols[i]
endfor

;; show similar curves for different params
;; f = dipper(x,0,0.01,2000,1,1)
;; oplot,x,f,linestyle=0,thick=thk,color=1
;; f = dipper(x,0.5,0.5,2000,1,1)
;; oplot,x,f,linestyle=0,thick=thk,color=1
;; f = dipper(x,1/1.3,1,5000,1,1.3)
;; oplot,x,f,linestyle=0,thick=thk,color=2
;; f = dipper(x,1,10,5000,1,10)
;; oplot,x,f,linestyle=0,thick=thk,color=3
;; f = dipper(x,1,20,5000,1,20)
;; oplot,x,f,linestyle=2,thick=thk,color=4
;; f = dipper(x,1,50,5000,1,50)
;; oplot,x,f,linestyle=2,thick=thk,color=5

axis,xrange=xr,/xst,xthick=thk1,xtickname=REPLICATE(' ',60)

;; show that it gives the same result as the Winn analyitic result
;oplot,x,1/!pi*(acos(x)-x*sqrt(1-x^2)),linestyle=2,thick=thk/2

device,/close
set_plot,'x'

end
