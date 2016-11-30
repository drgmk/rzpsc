;; look for periodicity, discrete autocorrelation seems the way to go

;; to use sigma-clipped mean in DACF compile dfc1.pro first!

function linear_cdf,x
  return,((x-10.)/(140.-10) > 0. ) < 1.
end

minpt = 2                      ; min number of pairs to use

loadct,0
set_plot,'ps'
!P.FONT=0
csz = 1.5
thk = 8
thk1 = 7
sqsz = 0.05
sqsig = 3
fn = 'auto.eps'
device,file=fn,/enc,/col,xsize=7,ysize=13,/inches,/times
xr = [10,155]
yrng = [0,560]
plot,[0],[0],xrange=xr,yrange=yrng,xthick=thk,ythick=thk,charsize=csz,/xst,/yst,xtitle='Lag (days)',ytitle='DACF power + offset',position=[0.1,0.06,0.94,0.95],xticklen=0.01,ytickname=REPLICATE(' ',60)

;; polyfill,[65,65,80,80],[yrng[0],yrng[1],yrng[1],yrng[0]],color=230,thick=2;,orientation=90,/LINE_FILL
axis,xaxis=1,xrange=xr,/xst,charsize=csz,xthick=thk,xticklen=0.01,xtitle='Lag (days)'
;; axis,xaxis=0,xrange=xr,/xst,charsize=csz,xthick=thk,xticklen=0.01,xtickname=REPLICATE(' ',60)

;; get data
readcol,'all-lc.txt',t1,f1,o,format='d,f,a'
tmp = WHERE(f1 lt 1.1 and t1 gt 2448012.5 and $
            ~( o eq 'AAVSO' and t1 < 2456413.5) and o ne 'Plates')
tmp = WHERE((o eq 'WASP' or o eq 'KELT') and f1 lt 1.15)
t1 = t1[tmp]
f1 = f1[tmp]

;; "modern" data, e.g. from 1989 when Herbst monitoring began or 2004 when WASP starts
;new = 2004
;getyeardata,t1,f1,new,/AFTER,tnew,fnew,told,fold
;dcf,tnew,fnew,tnew,fnew,lagnew,cornew,minlag=xr[0],maxlag=xr[1],minpt=minpt,numf=ROUND(MAX(tnew)-MIN(tnew))
;oplot,lagnew,cornew+195

;; old data
;dcf,told,fold,told,fold,lagold,corold,minlag=xr[0],maxlag=xr[1],minpt=minpt,numf=ROUND(MAX(told)-MIN(told))
;oplot,lagold,corold+190

;; all data
;dcf,t1,f1,t1,f1,lag3,cor3,minlag=xr[0],maxlag=xr[1],minpt=minpt,numf=2*ROUND(MAX(t1)-MIN(t1))
ok = WHERE(FINITE(cor3))
oplot,lag3[ok],cor3[ok]+410,thick=thk1
meanclip,cor3[ok],mn,sd
tmp = WHERE(cor3[ok] gt sqsig*sd or cor3[ok] lt -sqsig*sd)
oplot,lag3[ok[tmp]],cor3[ok[tmp]]+410,psym=6,symsize=sqsz,thick=thk1
xyouts,15,425,'All data',charsize=1.2,alignment=0

;; excluding 2004/6 data
;getyeardata,t1,f1,[2006],t456,f456,tout,fout
;dcf,tout,fout,tout,fout,lag4,cor4,minlag=xr[0],maxlag=xr[1],minpt=minpt,numf=2*ROUND(MAX(t1)-MIN(t1))
ok = WHERE(FINITE(cor4))
oplot,lag4[ok],cor4[ok]+380,thick=thk1
meanclip,cor4[ok],mn,sd
tmp = WHERE(cor4[ok] gt sqsig*sd or cor4[ok] lt -sqsig*sd)
oplot,lag4[ok[tmp]],cor4[ok[tmp]]+380,psym=6,symsize=sqsz,thick=thk1
xyouts,15,385,'excluding 2006',charsize=1.2,alignment=0

;; data on a yearly basis
distinct_colors,n_colors=12
yr0 = 2003
j=0
yoff = 5.
mx = 0
cors = PTRARR(100,/ALLOCATE_HEAP)
nb = 14
bin = FLTARR(nb)
binl = linarrw(min=5,max=145,n=nb+1)
for i=0,2015-yr0 do begin

   col = j mod 11 + 1
   yr = yr0 + i
   print,yr
   getyeardata,t1,f1,yr,t,f,nin=nin
   if nin le 1 then CONTINUE
;   tmp = WHERE(f gt 0.8)
;   f[tmp] = 1.0
   
   dcf,t,f,t,f,lag,cor,minlag=xr[0],maxlag=xr[1],minpt=minpt,numf=2.*(MAX(t)-MIN(t))
   *cors[i] = cor
   ok = WHERE(FINITE(cor),nok)
   sd = STDDEV(cor[ok])
   meanclip,cor[ok],mn,sd
   ntmp = 0
   if MAX(cor,/NAN) gt 4*sd or nok gt 100 then begin
      yoff += (mx - MIN(cor,/NAN) + 5) > 5.
      oplot,lag[ok],cor[ok]+yoff,color=col,thick=thk1
      tmp = WHERE(cor[ok] gt sqsig*sd or cor[ok] lt -sqsig*sd,ntmp)
      if tmp[0] ne -1 then oplot,lag[ok[tmp]],cor[ok[tmp]]+yoff,psym=6,color=col,symsize=sqsz,thick=thk1
      xyouts,15,yoff+3,STRN(yr),color=col,charsize=1.3;,alignment=1
      mx = MAX(cor,/NAN)
      j++
   endif

   ;; go through bins and add 1 if there's a positive peak in it
   if ntmp gt 0 then begin
      for k=0,nb-1 do begin
         for l=0,ntmp-1 do begin
            if lag[ok[tmp[l]]] gt binl[k] and lag[ok[tmp[l]]] le binl[k+1] then begin
               if cor[ok[tmp[l]]] gt 0 then begin
                  bin[k] += 1
                  break      ; only add 1 per year to each bin
               endif
            endif
         endfor
      endfor
   endif
      
   ;; for k=0,nb-1 do begin
   ;;    tmp = WHERE(lag gt binl[k] and lag le binl[k+1])
   ;;    if tmp[0] ne -1 and FINITE(MEAN(cor[tmp],/NAN)) then bin[k] += MEAN(cor[tmp],/NAN)
   ;; endfor

endfor

binc = (binl[0:-2]+binl[1:-1])/2.
yoff += (mx - MIN(bin) + 0.1) > 3.
oplot,binc,4*bin+520,psym=10,thick=thk,color=1
oplot,xr,FLTARR(nb+1)+520,linestyle=2,thick=thk,color=1
xyouts,100,543,'Yearly peak counting',charsize=1.2,color=1
xyouts,157,517,'0',charsize=csz,color=1
xyouts,157,537,'5',charsize=csz,color=1
xyouts,156,557,'10',charsize=csz,color=1

device,/close
set_plot,'x'
spawn,'pstopdf auto.eps'

delvarx,cumbin
for i=0,nb-1 do push,cumbin,replicate(binc[i],bin[i])
ksone,cumbin,'linear_cdf',d,p,/plot

end

!P.MULTI=[0,7,7]
yrng = [-5,20]
;plot,lag3,cor3,yrange=yrng
;legend,['all'],/top,/rig
;plot,lag4,cor4,yrange=yrng
;legend,['excl 2006'],/top,/rig
;plot,lagnew,cornew,yrange=yrng
;legend,['>'+STRN(new)],/top,/rig
;plot,lagold,corold,yrange=yrng
;legend,['<'+STRN(new)],/top,/rig
for i=0,99 do begin
   if N_ELEMENTS(*cors[i]) eq 0 then CONTINUE
   tmp = WHERE(FINITE(*cors[i]),nok)
   if nok gt 0 and MIN(*cors[i],/NAN) ne MAX(*cors[i],/NAN) then begin
      plot,lag,*cors[i],thick=1,yrange=yrng
;      oplot,lag,*cors[i],psym=6,symsize=0.1,thick=thk1
      legend,[STRN(yr0+i),STRN(nok),STRN(STDDEV(*cors[i],/NAN))],/top,/rig,box=0
      oplot,[68,68],yrng,linestyle=1
      oplot,[80,80],yrng,linestyle=1
   endif
endfor
!P.MULTI=0

end

;; using the model instead
readcol,'fitdata.txt',j,a,b,c,d,e,FORMAT='i,d,f,f,f,f',count=nl
delvarx,p
for i=0,N_ELEMENTS(a)-1 do push,p,[a[i],b[i],c[i],d[i],e[i]]
treg = linarrw(min=MIN(t1),max=max(t1),n=2e5)
fmod = dipwrap(treg,p)
b = A_CORRELATE(fmod,lag)
oplot,lag*dt,b,linestyle=2,thick=4

end
