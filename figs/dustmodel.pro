;; show some properties of the dust model

!P.FONT=0
csz = 1.5
thk1 = 8
thk2 = 4
thk = 12
distinct_colors,n_colors=12
set_plot,'ps'
device,file='model.eps',/enca,xsize=8,ysize=6,/inches,/times,/col
xr = [-4,4]
yr = [-1,1]

plot,[0],[0],xrange=xr,yrange=yr,position=[0.02,0.12,0.98,0.99],xst=5,yst=5
axis,xrange=xr,xtitle=textoidl('x (units of R_{star})'),xthick=thk1,charsize=csz,xst=1
;axis,yaxis=0,yrange=yr;,ytickname=REPLICATE(' ',60),ythick=thk1,charsize=csz

;; the star
;tvcircle,0.3,0.25,0.6,/data,thick=thk,color=7
smx = 0.3
polyfill,[-1,-1,1,1],[yr[0],smx,smx,yr[0]],/line_fill,orientation=45,color=7
oplot,[-1,-1],[yr[0],smx],thick=thk,color=7
oplot,[1,1],[yr[0],smx],thick=thk,color=7
xyouts,0,0.21,'star',alignment=0.5,charsize=csz,color=7,charthick=thk
arrow,0.2,0.23,1,0.23,color=7,thick=thk1,/data,hsize=500
arrow,-0.2,0.23,-1,0.23,color=7,thick=thk1,/data,hsize=500

;; the Gaussian
x = linarrw(min=xr[0]-2,max=xr[1]+2,n=200)
g = taufunc_g(x,1,0.15,1)
oplot,x+0.6,g-0.9,thick=thk,color=5
oplot,x-1.5,g-0.9,thick=thk,color=5
oplot,x+3,g-0.9,thick=thk,color=5
arrow,0.7,1.0-0.9,1.7,1.0-0.9,thick=thk1,/data,hsize=500,color=5
xyouts,-2.95,-0.6,'optical depth',orientation=62,color=5,charsize=csz

;; light curve
lc = dipper(x,0,0.15,1,0.8,1)
oplot,x,lc-0.01,thick=thk,color=1
xyouts,-3,0.9,'light curve',alignment=0.5,charsize=csz,color=1

device,/close
set_plot,'x'

end
