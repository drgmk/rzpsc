;; show yearly light curves

pro plotem,t,f,o,os

  distinct_colors,n_colors=12
  for i=0,N_ELEMENTS(os)-1 do begin
     col = i
     if os[i] eq 'KELT' then col = 4
     if os[i] eq 'WASP' then col = 10
     a = WHERE(o eq os[i],na)
;     if os[i] eq 'Home94' then fill = 0 else fill = 1
;     if na gt 0 then oplotcircw,t[a],f[a],symsize=0.3,col=col+1,fill=1
     if na gt 0 then oplot,t[a],f[a],symsize=0.1 ,col=col+1,psym=6,thick=4
  endfor
end

readcol,'all-lc.txt',t1,f,o,format='d,f,a'

;; uncomment this and comment chunk below for just WASP/KELT lcs
tmp = WHERE(o eq 'WASP' or o eq 'KELT')
t1 = t1[tmp]
f = f[tmp]
o = o[tmp]
years = [2004,2006,2007,2008,2009,2010,2011,2012,2013,2014]
cols = [5,11]

os = (o[UNIQ(o,SORT(o))])
no = N_ELEMENTS(os)

event = [2453126.5+52,$         ; 2004
         2453210.25,$
         2453240,$
         2453260.25,$
         
         2453996.5,$            ; 2006
;         2454022,$
         2454068.25,$
;         2454101,$

         2454391.75,$           ; 2007
         2454396,$
         2454415.5,$
         2454435.,$
         2454457.75,$

         2454739.25,$           ; 2008
;         2454789.5,$
         2454829.75,$
         2454846,$
         2454859.25,$
         2454883.75,$

         2455100.25,$           ; 2009
         2455122.75,$
         2455152.5,$
         2455177,$
         2455213,$

         2455481.25,$           ; 2010
         2455529.25,$
;         2455507,$
         2455571,$

         2455842.5,$            ; 2011
         2455864.25,$
         2455884,$
         2455894.5,$
         2455967,$

         2456195,$              ; 2012
         2456209.25,$
         2456223.25,$
         2456236.25,$
         2456263.25,$

         2456568,$              ; 2013
         2456581,$
         2456600.25,$
         2456633.25,$
         2456653.25,$
         2456662,$

         2456933.25,$           ; 2014
         2456947.25,$
         2456975.25,$
        2456995]

win   = [5.,$                   ; 2004
         5,$
         2,$
         7.25 ,$
         
         7,$                    ; 2006
;         2,$
         12,$
;         2.5,$

         2.5,$                    ; 2007
         2,$
         3,$
         3,$
         8.5,$

         6,$                    ; 2008
;         4,$
         2.5,$
         3,$
         2,$
         2,$

         4,$                    ; 2009
         2.5,$
         5.5,$
         5.2,$
         2,$

         2,$                    ; 2010
         2,$
;         12.3,$
         3,$

         2.5,$                  ; 2011
         2,$
         2,$
         3,$
         5,$

         3,$                    ; 2012
         2,$
         2,$
         2,$
         2.5,$

         2,$                    ; 2013
         2,$
         4,$
         3,$
         2,$
         4,$

         3,$                    ; 2014
         3,$
         3,$
         10]

nevents = N_ELEMENTS(event)

;; model
readcol,'fitdata.txt',j,a,b,c,d,e,FORMAT='i,d,f,f,f,f',count=nl
 
set_plot,'ps'
!P.MULTI=[0,5,10]
!P.FONT=0
csz = 1.4
thk = 8
fn = 'yearly-zoom.eps'
print,fn
device,file=fn,/enc,/col,xsize=6,ysize=8,/inches,/times
   
yp = REVERSE(linarrw(min=0.04,max=0.995,n=11))
off = 0.005

for l=0,9 do begin

   jdcnv,years[l],5,1,0,jd0     ; 1st May of each year
   inyr = WHERE(event gt jd0 and event lt jd0+365.25,nevyr)
   print,years[l],jd0,nevyr
   yr = [0.,1.19]
;   if l eq 0 then yr = [0.6,1.1]

   xoff1 = 0.005
   for m=0,nevyr-1 do begin

      event1 = event[inyr[m]]
      win1 = win[inyr[m]]
      
      tlo = event1-win1
      thi = event1+win1
      xr = [tlo,thi]

      in = WHERE(t1 gt tlo and t1 lt thi and t1 gt jd0 and t1 lt jd0+365.25,nin)
      if nin eq 0 then continue ; stop,'no data in window'

      xoff2 = xoff1 + 0.2 *win[inyr[m]]/4.
;      if m eq nevyr-1 then xoff2 = 0.99
      plot,[0],[0],xrange=xr,yrange=yr,/yst,/xst,position=[xoff1,yp[l+1]+off,xoff2,yp[l]],xtickname=REPLICATE(' ',60),ytickname=REPLICATE(' ',60),xticklen=0.05,thick=thk,xtickinterval=1,xminor=1,ytickinterval=0.25

      xoff1 = xoff2 + off

      loadct,0,/silent
      col = 150
      thk1 = 2
      oplot,xr,[1,1],linestyle=1,color=col,thick=thk1
      oplot,xr,[1,1]*0.75,linestyle=1,color=col,thick=thk1
      oplot,xr,[1,1]*0.5,linestyle=1,color=col,thick=thk1
      oplot,xr,[1,1]*0.25,linestyle=1,color=col,thick=thk1

      plotem,t1[in],f[in],o[in],os
      if m eq 0 then legend,[STRING(years[l],FOR='(I4)')],/bot,/rig,margin=0,box=0,charsize=0.5
;      if m eq 0 then legend,[STRING(jd0,FOR='(F14.1)')],/bot,/lef,margin=0,box=0,charsize=0.5

      ;; add model
      delvarx,p
      for j=0,nl-1 do if a[j] gt tlo and a[j] lt thi then push,p,[a[j],b[j],c[j],d[j],e[j]]
      x = linarrw(min=0,max=365.25,n=2e3)
;      oplot,x,dipwrap(x+t0here+t0,p),color=5,thick=4
      
;      oplot,a[41]-t0 + (FINDGEN(100)-50)*(130) - t0here,FLTARR(100)+0.1,psym=1
      
      if l eq 9 and m eq nevyr-1 then begin
         axis,xaxis=0,xtitle='Days + offset',xrange=xr-xr[0],/xst,charsize=csz
      endif 

      if l eq 6 and m eq nevyr-1 then axis,yaxis=1,ytitle='Relative Flux',yrange=yr,/yst,charsize=csz,ytickinterval=0.25

   endfor
endfor

!P.MULTI=0
device,/close
set_plot,'x'
spawn,'pstopdf '+fn

end
