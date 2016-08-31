;; collect all data and flatten a bit

dir = '../data/'

;; normalise roughly to 1, will do properly afterwards

;; WASP
restore,dir+'wasp1-s1-lc.xdr'
t = t1
f = f1/22.
o = STRARR(N_ELEMENTS(t))+'WASP'
restore,dir+'wasp1-s2-lc.xdr'
t = [t,t1]
f = [f,f1/20.5]
o = [o,STRARR(N_ELEMENTS(t1))+'WASP']

;; Catalina
readcol,dir+'result_web_file31zGE3.csv',i,m,e,r,de,d,format='l,f,f,d,d,d'
f1 = 10^(-0.4*m)/2.82e-5
d += 2400000.5d
tmp = WHERE(f1 lt 1.3)
t = [t,d[tmp]]
f = [f,f1[tmp]]
o = [o,STRARR(N_ELEMENTS(d[tmp]))+'Catalina']

;; KELT public
;; readcol,dir+'KELT_N02_lc_030387_V01_west_raw_lc.tbl',d,m,e,format='d,f,f'
;; readcol,dir+'KELT_N02_lc_034612_V01_east_raw_lc.tbl',d1,m1,e1,format='d,f,f'
;; f1 = 10^(-0.4*m)/4.7e-7
;; f2 = 10^(-0.4*m1)/5e-7
;; t = [t,d,d1]
;; f = [f,f1,f2]
;; o = [o,STRARR(N_ELEMENTS([d,d1]))+'KELT']

;; KELT from Josh/Joey
;; readcol,dir+'RZ_Psc_KELT.combined.raw',d1,m1,e1,format='d,f,f'
;; f1 = 10^(-0.4*m1)/4.7e-7
;; t = [t,d1]
;; f = [f,f1]
;; o = [o,STRARR(N_ELEMENTS([d1]))+'KELT']
readcol,dir+'lc_30387_mag.data.raw',d1,m1,e1,format='d,f,f'
readcol,dir+'lc_34612_mag.data.raw',d2,m2,e2,format='d,f,f'
f1 = 10^(-0.4*m1)/4.7e-7
f2 = 10^(-0.4*m2)/4.9e-7
t = [t,d1,d2]
f = [f,f1,f2]
o = [o,STRARR(N_ELEMENTS([d1,d2]))+'KELT']

;; Home MVS, from http://www.4pisysteme.de/observatory/observatory_3.html
;; readcol,dir+'home/home-data.txt',t1,m1,t2,m2,t3,m3,t4,m4,t5,m5,t6,m6,t7,m7,t8,m8,format='d,a,d,a,d,a,d,a,d,a,d,a,d,a,d,a',delimiter='	'
;; ttmp = [t1,t2,t3,t4,t5,t6,t7,t8]+2400000d
;; mtmp = [m1,m2,m3,m4,m5,m6,m7,m8]
;; nt = N_ELEMENTS(ttmp)
;; unc = INTARR(nt)
;; lim = INTARR(nt)
;; unc = WHERE( STREGEX(mtmp,':',/BOOL) )
;; lim = WHERE( STREGEX(mtmp,'L',/BOOL) )
;; m = FLTARR(nt)
;; for i=0,nt-1 do m[i] = FLOAT( STR_REPLACE(STR_REPLACE(mtmp[i],':',''),'L','') )
;; f1 = 10^(-0.4*m)/1.27e-5
;; t = [t,ttmp]
;; f = [f,f1]
;; o = [o,STRARR(N_ELEMENTS(ttmp))+'Home94']

;; Home + Gurtler, Sonneburg and Harvard plates
readcol,'../data/RZP-HASO.txt',d,m,format='d,f'
d += 2400000d
f1 = 10^(-0.4*m)/1.27e-5
t = [t,d]
f = [f,f1]
o = [o,STRARR(N_ELEMENTS(d))+'Plates']

;; Zajtseva
readcol,dir+'zajtseva/data.txt',t1,m1,format='d,f'
t1 += 2400000
f1 = 10^(-0.4*m1)/2.4e-5
t = [t,t1]
f = [f,f1]
o = [o,STRARR(N_ELEMENTS(t1))+'Zajtseva85']

;; Karetnikov, can get from http://www.konkoly.hu/IBVS/0701.html
readcol,dir+'Karetnikov73.txt',jd,hr,mn,v,format='d,f,f,f'
jd += 2441549.5d
f1 = 10^(-0.4*v)/2.4e-5
t = [t,jd]
f = [f,f1]
o = [o,STRARR(N_ELEMENTS(jd))+'Karetnikov73']

;; Kardopolov
readcol,dir+'Kardopolov.txt',jd,v,format='d,f'
jd += 2443000.d
f1 = 10^(-0.4*v)/2.2e-5
t = [t,jd]
f = [f,f1]
o = [o,STRARR(N_ELEMENTS(jd))+'Kardopolov80']

;; Kiselev 1991
readcol,dir+'Kiselev91.txt',jd,v,format='d,f'
jd += 2447000.d
f1 = 10^(-0.4*v)/2.15e-5
t = [t,jd]
f = [f,f1]
o = [o,STRARR(N_ELEMENTS(jd))+'Kiselev91']

;; Potravnov 2014
readcol,dir+'Potravnov14.txt',jd,v,format='d,f'
jd += 2400000.d
f1 = 10^(-0.4*v)/2.3e-5
t = [t,jd]
f = [f,f1]
o = [o,STRARR(N_ELEMENTS(jd))+'Potravnov14']

;; Herbst 1999
readcol,dir+'herbst99.txt',jd,v,format='d,f',delimiter='|'
jd += 2440000.d
f1 = 10^(-0.4*v)/2.15e-5
t = [t,jd]
f = [f,f1]
o = [o,STRARR(N_ELEMENTS(jd))+'Herbst99']

;; AAVSO, "tidied" with tidy_aavso.pro
readcol,dir+'aavso-tidy.txt',jd,mag,format='d,f'
f1 = 10^(-0.4*mag)/2.2e-5
t = [t,jd]
f = [f,f1]
o = [o,STRARR(N_ELEMENTS(jd))+'AAVSO']

;; ASAS
readcol,dir+'asas.txt',hjd,m1,format='d,f'
hjd += 2450000.d
f1 = 10^(-0.4*m1)/2.312e-5
t = [t,hjd]
f = [f,f1]
o = [o,STRARR(N_ELEMENTS(hjd))+'ASAS']

srt = SORT(t)
t = t[srt]
f = f[srt]
o = o[srt]

!P.MULTI=[0,1,2]
;plot,t,f,psym=3,/xst            ; original curve

;; now tidy up
readcol,'fitdata.txt',j,a,b,c,d,e,FORMAT='i,d,f,f,f,f',count=nl
delvarx,p
for i=0,N_ELEMENTS(a)-1 do push,p,[a[i],b[i],c[i],d[i],e[i]]
;oplot,t,dipwrap(t,p),linestyle=2

;plot,t,f-dipwrap(t,p),psym=3

;; filter by survey and individual years
jdcnv,1840,5,1,0,t0             ; 1st May 1840, leave earlier data alone
jdcnv,1974,5,1,0,t1974          ; 1st May 1974, treat earlier data roughly
os = o[UNIQ(o,SORT(o))]
clipsig = 3.0                   ; for sigma clipping
for i=0,N_ELEMENTS(os)-1 do begin

   print,os[i]

   thiso = WHERE(o eq os[i])
   
   ttmp = t[thiso]
   ftmp = f[thiso]-1;dipwrap(t[thiso],p) ; data are now around zero

   for j=0,200 do begin

      ;; select only data near median
      thisyr = WHERE(ttmp gt t0+365.25d*j and ttmp le t0+365.25d*(j+1))
      if os[i] eq 'Plates' or t0+365.25d*j lt t1974 then begin
         thisyruse = WHERE(ttmp gt t0+365.25d*j and ttmp le t0+365.25d*(j+1) and ftmp gt -1. and ftmp lt 1.,ntmp)
         print,'PLates or <1974'
      endif else begin
         thisyruse = WHERE(ttmp gt t0+365.25d*j and ttmp le t0+365.25d*(j+1) and ftmp gt -0.5 and ftmp lt 0.2,ntmp)
         print,'gt 1974'
      endelse
      if ntmp lt 2 then CONTINUE; do nothing if too few points
      tyr = ttmp[thisyruse]
      fyr = ftmp[thisyruse]
      tspan = MAX(tyr)-MIN(tyr)
      plot,tyr,fyr,psym=1
      legend,[os[i],STRN(1840+j)],/top,/rig

      ;; fine-tune median to 1, ntmp=1e5 means we do this for all data
      if ntmp lt 1e5 or os[i] eq 'Plates' then begin
         meanclip,fyr,mn,clipsig,subs=ok
         med = MEDIAN(fyr[ok])
         f[thiso[thisyr]] /= med + 1
         oplot,MINMAX(tyr),med*[1,1]
         oplot,tyr[ok],fyr[ok],psym=6         
         plot,t[thiso[thisyr]],f[thiso[thisyr]],psym=1

      ;; poly fitting, not used in the end
      endif else begin
         
;         ply = POLY_FIT(tyr-MEAN(tyr),fyr,3,yfit)
         ply = ROBUST_POLY_FIT(tyr-MEAN(tyr),fyr,2,yfit,clipsig)
         oplot,tyr,POLY(tyr-MEAN(tyr),ply),linestyle=2
;      plot,tyr,fyr-POLY(tyr-MEAN(tyr),ply),psym=3
;      oplot,MINMAX(tyr),[0,0],linestyle=2     

      ;; now fix this part of the lc
         f[thiso[thisyr]] /= POLY(ttmp[thisyr]-MEAN(tyr),ply) + 1

         plot,tyr,f[thiso[thisyr]],psym=1,/yst
      endelse
      
      oplot,MINMAX(tyr),[1,1],linestyle=2     
      
      wait,0.01
   endfor
endfor

openw,lun,'all-lc.txt',/GET_LUN
printf,lun,'JD,flux,obs'
for i=0,N_ELEMENTS(t)-1 do printf,lun,STRING(t[i],format='(d)')+','+STRING(f[i],format='(f)')+','+o[i]
free_lun,lun

!P.MULTI=0
end
