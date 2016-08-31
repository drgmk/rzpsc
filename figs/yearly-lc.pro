;; show yearly light curves

pro plotem,t,f,o,os

  distinct_colors,n_colors=12
  for i=0,N_ELEMENTS(os)-1 do begin
     col = i
     if os[i] eq 'KELT' then col = 4
     if os[i] eq 'WASP' then col = 10
     a = WHERE(o eq os[i] and t gt 0 and t lt 365.25,na)
;     if os[i] eq 'Home94' then fill = 0 else fill = 1
;     if na gt 0 then oplotcircw,t[a],f[a],symsize=0.4,col=col+1,fill=1
     if na gt 0 then oplot,t[a],f[a],symsize=0.1 ,col=col+1,psym=6,thick=4
  endfor
end

readcol,'all-lc.txt',t1,f,o,format='d,f,a'

;; uncomment this and comment chunk below for just WASP/KELT lcs
tmp = WHERE(o eq 'WASP' or o eq 'KELT')
t1 = t1[tmp]
f = f[tmp]
o = o[tmp]
firstyr = 2004
ndec = 1.
nyrs = 11.
cols = [5,11]

;; daycnv,MIN(t1),yr
;; firstyr = 10*FLOOR(yr/10.)
;; daycnv,MAX(t1),yr
;; lastyr = 10*CEIL(yr/10.)
;; nyrs = 10.
;; ndec = (lastyr-firstyr)/FLOAT(nyrs)

os = (o[UNIQ(o,SORT(o))])
no = N_ELEMENTS(os)
;cols = FINDGEN(no)+1

;; model
readcol,'fitdata.txt',j,a,b,c,d,e,FORMAT='i,d,f,f,f,f',count=nl

openw,lun,'../yearlyfigs.tex',/GET_LUN
for k=0,ndec-1 do begin

   yearzero = firstyr + k*nyrs
   jdcnv,yearzero,5,1,0,t0      ; 1st May of each year
   t = t1 - t0
   
   set_plot,'ps'
   !P.MULTI=[0,1,nyrs]
   !P.FONT=0
   csz = 1.5
   thk = 8
   fn = 'yearly-'+STRN(yearzero,FORMAT='(I4)')+'.eps'
   print,fn
   device,file=fn,/enc,/col,xsize=6,ysize=8,/inches,/times
   
   yp = REVERSE(linarrw(min=0.04,max=0.995,n=11))

   i = 0
   for l=0,nyrs-1 do begin
      
      xr = [0,365.25]
      yr = [0.,1.4]
      xr = [20,325.25] ; uncomment these
      yr = [0.,1.3]
      t0here = 365.25*l
      tmp = WHERE(t-t0here gt 0 and t-t0here lt 365.25,nin)
      if nin eq 0 then begin
         continue               ; uncomment this too
      endif
      plot,[0],[0],xrange=xr,yrange=yr,/yst,/xst,position=[0.07,yp[i+1],0.995,yp[i]],xtickname=REPLICATE(' ',60),ytickname=REPLICATE(' ',60),xticklen=0.05,thick=thk

      oplot,xr,[1,1],linestyle=1
      oplot,xr,[1,1]*0.5,linestyle=1

      plotem,t-t0here,f,o,os
      jdcnv,yearzero+l,5,1,0,jd0 ; 1st May of each year
      legend,[STRING(yearzero+l,FOR='(I4)')],/bot,/lef,margin=0,box=0,charsize=1.
      legend,[STRING(jd0,FOR='(F14.1)')],/bot,/lef,margin=0,box=0,charsize=0.5

      ;; add model
      delvarx,p
      for j=0,nl-1 do if a[j] gt t0+t0here and a[j] lt t0+t0here+366 then push,p,[a[j],b[j],c[j],d[j],e[j]]
      x = linarrw(min=0,max=365.25,n=2e3)
;      oplot,x,dipwrap(x+t0here+t0,p),color=5,thick=4
      
;      oplot,a[41]-t0 + (FINDGEN(100)-50)*(130) - t0here,FLTARR(100)+0.1,psym=1
      
      if i eq 9 then legend,os,colors=cols,box=0,/bot,/rig,psym=6,symsize=FLTARR(no)+.1,textcolors=cols,thick=thk,margin=0.,charsize=0.53
      
      if l eq nyrs-1 then axis,xaxis=0,xtitle='Days since 1 May each year',xrange=xr,/xst,charsize=csz
      if l eq nyrs-1 then axis,yaxis=0,ytitle='Relative Flux',yrange=yr,/yst,charsize=csz

      i += 1
   endfor
   
   !P.MULTI=0
   device,/close
   set_plot,'x'
   spawn,'pstopdf '+fn
;spawn,'rm '+fn

   printf,lun,'\begin{figure*}'
   printf,lun,'\begin{center}'
   printf,lun,'\hspace{-0.25cm} \includegraphics[width=1\textwidth]{figs/'+fn+'}'
   printf,lun,'\caption{V-band photometry for RZ Psc. Rows begin on 1 May of the year shown in each panel, with the initial Julian date also shown.}\label{fig:yearly}'
   printf,lun,'\end{center}'
   printf,lun,'\end{figure*}'
   printf,lun,''
   
endfor
free_lun,lun

end
