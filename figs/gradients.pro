;; look at lc gradients

readcol,'all-lc.txt',t1,f,o,format='d,f,a'

tmp = WHERE((o eq 'WASP' or o eq 'KELT') and f lt 1.15)
t1 = t1[tmp]
f = f[tmp]
o = o[tmp]
wasp = WHERE(o eq 'WASP',complement=kelt)
e = f
e[wasp] = 0.02
e[kelt] = 0.02

dt = t1[1:-1] - t1[0:-2]
tmp = WHERE(dt gt 0.5) + 1
tbr = t1[tmp]
tbr = [min(t1),tbr,max(t1)]
;tbr = tbr[SORT(tbr)]
nbr = N_ELEMENTS(tbr)

!P.MULTI = [0,30,25]
plot,t1,f,psym=3
times = DBLARR(nbr)-99
slopes = FLTARR(nbr)-99
uncs = FLTARR(nbr)-1
chisq = FLTARR(nbr)-1
meanf = FLTARR(nbr)+100
sunset_colors
for i=0,nbr-2 do begin
   in = WHERE(t1 ge tbr[i] and t1 lt tbr[i+1],nin)
   if nin le 5 then begin
      plot,[t1[in]],[f[in]],/xst,/yst
      continue
   endif
   fit = LINFIT(t1[in]-mean(t1[in]),f[in],MEASURE_ERR=e[in],SIGMA=coef,chisq=chi2)
;   fit = ROBUST_LINEFIT(t1[in]-mean(t1[in]),f[in],yfit,sig,coef)
;   chi2 = TOTAL( ((f[in] - ((t1[in]-mean(t1[in]))*fit[1] + fit[0]))/(FLTARR(nin)+0.02))^2 )
   col = 50*chi2/float(nin-2)
   plot,t1[in],f[in],color=col,psym=1,/xst,/yst,yrange=[0,1.1],xrange=mean(t1[in])+[-.25,.25],xtickname=replicate(' ',60),ytickname=replicate(' ',60)
  oplot,t1[in],(t1[in]-mean(t1[in]))*fit[1] + fit[0],color=col
   chisq[i] = chi2/float(nin-2)
   if chisq[i] lt 100 then begin
      times[i] = mean(t1[in])
      uncs[i] = coef[1]
      slopes[i] = fit[1]
      meanf[i] = min(f[in])
   endif

endfor
!P.MULTI=0

distinct_colors,n_colors=12
set_plot,'ps'
!P.FONT=0
csz = 1.8
thk = 8
device,file='gradients.eps',/enc,/col,xsize=8,ysize=6,/inche,/times
xr = [-5.4,2.9]
yr = [-0.09,1.15]
plot,[0],[0],xrange=xr,yrange=yr,xthick=thk,ythick=thk,charsize=csz,xtitle='Gradient (1/day)',ytitle='Minimum flux',position=[0.13,0.135,0.99,0.99],/yst,/xst
noise = 0.3
y = linarrw(min=0,max=1,n=500)
gr = [0.3,1.,10,100]
col = [11,9,3,7]
for i=0,N_ELEMENTS(gr)-1 do begin
   oplot,sqrt(8.7/gr[i])*(1-y),y,thick=thk,color=col[i]
   oplot,-sqrt(8.7/gr[i])*(1-y),y,thick=thk,color=col[i]
endfor
for i=0,nbr-1 do begin
   if uncs[i]/abs(slopes[i]) lt 0.25 then begin
      oplotcircw,slopes[i],meanf[i],symsize=0.05/uncs[i],col=1,thick=6
   endif else oplotcircw,[slopes[i]],[meanf[i]],symsize=0.5,thick=2,col=2,/fill
;   xyouts,slopes[i],meanf[i],STRN(times[i],format='(F10.2)'),orientation=90
endfor
legend,textoidl(['0.3au','1au','10au','100au']),colors=col,/top,/lef,box=0,margin=0,charsize=1.4,linestyle=0,thick=thk
device,/close
set_plot,'x'

;; estimate semimajor axes
sma = 8.7 * ( (1-meanf)/slopes )^2

set_plot,'ps'
!P.FONT=0
csz = 1.8
thk = 8
device,file='sma.eps',/enc,/col,xsize=8,ysize=6,/inche,/times
xr = [2e-3,200]
plot,[0],[0],xrange=xr,yrange=yr,xthick=thk,ythick=thk,charsize=csz,xtitle=textoidl('a_{circ} (au)'),ytitle='Minimum flux',position=[0.13,0.135,0.99,0.99],/xl,/yst,/xst,xtickformat='exponentsonaxis'
;polyfill,[xr[0],xr[0],5e-3,5e-3],[yr[0],yr[1],yr[1],yr[0]],/line_fill,orientation=45,color=7,thick=4
;oplot,5e-3*[1,1],yr,color=7,thick=thk
oplot,8.7 * ( (1-(y))/noise )^2,y,thick=thk,color=3
oplot,8.7 * ( (1-(y))/1. )^2,y,thick=thk,color=6
;tmp = REVERSE(WHERE(8.7 * ( (1-y)/42. )^2 gt xr[0]))
;polyfill,[xr[0],8.7 * ( (1-y[tmp])/42. )^2,5e-3,xr[0]],[0,y[tmp],0],/line_fill,orientation=45,color=7,thick=4
;oplot,8.7 * ( (1-y)/42. )^2,y,thick=thk,color=6
for i=0,nbr-1 do begin
   if uncs[i]/abs(slopes[i]) lt 0.25 then begin
      oplotcircw,sma[i],meanf[i],symsize=0.05/uncs[i],col=1,thick=6
;      xyouts,sma[i],meanf[i],times[i]
   endif else oplotcircw,[sma[i]],[meanf[i]],symsize=0.5,thick=2,col=2,/fill
endfor
legend,['Detection limit','Gradient=1/day'],/cen,/lef,box=0,charsize=1.4,linestyle=0,thick=thk,colors=[3,6],margin=0
oplot,5e-3*[1,1],yr,thick=thk,linestyle=2,color=FSC_COLOR('gray')
device,/close
set_plot,'x'

;; some stats
fs = [0.,0.5,0.8,1.15]
n = N_ELEMENTS(fs)
for i=0,n-2 do begin
   tmp = WHERE(meanf gt fs[i] and meanf le fs[i+1] and slopes lt -0.03,nlo)
   tmp = WHERE(meanf gt fs[i] and meanf le fs[i+1] and slopes ge -0.03,nhi)
   print,fs[i],fs[i+1],nlo,nhi
endfor

end
