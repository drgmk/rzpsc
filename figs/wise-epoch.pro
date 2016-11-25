;; show WISE photometry of RZ Psc

readcol,'all-wise.csv',mjd,w1,w1e,w2,w2e,w3,w3e,w4,w4e,format='f,f,f,f,f,f,f,f,f'

mjd += 2400000.5
daycnv,mjd,y,m,d,h
mjd = y + m/12. + d/12./30 + h/12./30./24.

set_plot,'ps'
!P.FONT=0
thk = 4
csz = 1.7
ssz = 0.2
distinct_colors
device,file='wise-var.eps',/col,/enc,/inch,xsize=12,ysize=3,/times
!P.MULTI = [0,4,1]

xr = [2009.7,2016.2]
yr = [9.05,9.55]

;; total/star flux ratios
rw1 = 1.44
rw2 = 2.42
rw3 = 22.1
rw4 = 76.3

xp = linarrw(min=0.05,max=0.99,n=5)

xarr = linarrw(min=4,max=10,n=100)

plot,[0],[0],psym=1,xrange=xr,yrange=yr,xtitle='years since 2010 Jan 13 ',ytitle='W1',xthick=thk,ythick=thk,charsize=csz,position=[xp[0],0.13,xp[1],0.99],/xst,/yst
oplot,mjd,w1,psym=6,color=1,symsize=ssz,thick=thk

xr = minmax(w2[WHERE(w2 gt 0)]) * [0.99,1.01]
plot,[0],[0],psym=1,xrange=xr,yrange=yr,xtitle='W2',xthick=thk,ythick=thk,charsize=csz,position=[xp[1],0.13,xp[2],0.99],ytickname=REPLICATE(' ',60),/xst,/yst
oplot,xarr,xarr*(1-1/rw1)/(1-1/rw2)+4.66,linestyle=2,color=2,thick=thk
oplot,w2,w1,psym=6,color=1,symsize=ssz,thick=thk

xr = minmax(w3[WHERE(w3 gt 0)]) * [0.99,1.01]
plot,[0],[0],psym=1,xrange=xr,yrange=yr,xtitle='W3',xthick=thk,ythick=thk,charsize=csz,position=[xp[2],0.13,xp[3],0.99],ytickname=REPLICATE(' ',60),/xst,/yst
oplot,xarr,xarr*(1-1/rw1)/(1-1/rw3)+7.33,linestyle=2,color=2,thick=thk
oplot,w3,w1,psym=6,color=1,symsize=ssz,thick=thk

xr = minmax(w4[WHERE(w4 gt 0)]) * [0.99,1.01]
plot,[0],[0],psym=1,xrange=xr,yrange=yr,xtitle='W4',xthick=thk,ythick=thk,charsize=csz,position=[xp[3],0.13,xp[4],0.99],ytickname=REPLICATE(' ',60),/xst,/yst
oplot,xarr,xarr*(1-1/rw1)/(1-1/rw2)+6.7,linestyle=2,color=2,thick=thk
oplot,w4,w1,psym=6,color=1,symsize=ssz,thick=thk

device,/close
set_plot,'x'

!P.MULTI=0

end
