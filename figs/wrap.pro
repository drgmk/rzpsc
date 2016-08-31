;; try phase folding light curves

;; WASP (from nexsci)
restore,'wasp1-s1-lc.xdr'
t = t1
f = f1/22.
o = STRARR(N_ELEMENTS(t))+'WASP'
restore,'wasp1-s2-lc.xdr'
t = [t,t1]
f = [f,f1/20.5]
o = [o,STRARR(N_ELEMENTS(t1))+'WASP']

;; WASP from Amy
readcol,'transit_stplim6mfile_0_objid_91.txt',d,m,e
f1 = 10^(-0.4*m)/2.9e-5
d += 2400000.5
t = [t,d]
f = [f,f1]
o = [o,STRARR(N_ELEMENTS(d))+'WASP']

;; Catalina
readcol,'result_web_file31zGE3.csv',i,m,e,r,de,d
f1 = 10^(-0.4*m)/2.82e-5
d += 2400000.5
tmp = WHERE(f1 lt 1.3); throw out some super high points
t = [t,d[tmp]]
f = [f,f1[tmp]]
o = [o,STRARR(N_ELEMENTS(d))+'Catalina']

;; KELT
readcol,'KELT_N02_lc_030387_V01_west_raw_lc.tbl',d,m,e,format='d,f,f'
readcol,'KELT_N02_lc_034612_V01_east_raw_lc.tbl',d1,m1,e1,format='d,f,f'
f1 = 10^(-0.4*m)/4.7e-7
f2 = 10^(-0.4*m1)/5e-7
t = [t,d,d1]
f = [f,f1,f2]
o = [o,STRARR(N_ELEMENTS([d,d1]))+'KELT']

;; Home MVS
readcol,'home/home-data.txt',t1,m1,t2,m2,t3,m3,t4,m4,t5,m5,t6,m6,t7,m7,t8,m8,format='f,a,f,a,f,a,f,a,f,a,f,a,f,a,f,a',delimiter='	'
ttmp = [t1,t2,t3,t4,t5,t6,t7,t8]+2400000
mtmp = [m1,m2,m3,m4,m5,m6,m7,m8]
nt = N_ELEMENTS(ttmp)
unc = INTARR(nt)
lim = INTARR(nt)
unc = WHERE( STREGEX(mtmp,':',/BOOL) )
lim = WHERE( STREGEX(mtmp,'L',/BOOL) )
m = FLTARR(nt)
for i=0,nt-1 do m[i] = FLOAT( STR_REPLACE(STR_REPLACE(mtmp[i],':',''),'L','') )
f1 = 10^(-0.4*m)/1.27e-5
;t = [t,ttmp]
;f = [f,f1]
;o = [o,STRARR(N_ELEMENTS(ttmp))+'Home']

;; Zajtseva
readcol,'zajtseva/data.txt',t1,m1,format='f,f'
t1 += 2400000
f1 = 10^(-0.4*m1)/2.4e-5
;t = [t,t1]
;f = [f,f1]
;o = [o,STRARR(N_ELEMENTS(t1))+'Zajtzeva']

;; AAVSO
readcol,'aavso1.csv',jd,mag
;tmp = WHERE(jd gt 2455093)
tmp = WHERE(jd gt 0)
jd = jd[tmp]
mag = mag[tmp]
f1 = 10^(-0.4*mag)/2.12e-5
t = [t,jd]
f = [f,f1]
o = [o,STRARR(N_ELEMENTS(jd))+'AAVSO']

plot,t,f,psym=3,/xst
openw,lun,'all-lc.txt',/GET_LUN
printf,lun,'JD,flux,obs'
for i=0,N_ELEMENTS(t)-1 do printf,lun,STRING(t[i],format='(d)')+','+STRING(f[i],format='(f)')+','+o[i]
free_lun,lun

end

dstep = 0.1
dr = [5,100]
nd = (dr[1]-dr[0])/dstep+1
s = linarrw(min=dr[0],max=dr[1],n=nd)

mmin = TOTAL(ABS(f[1:nt-1]-f[0:nt-2]))
m = FLTARR(nd)
for i=0,nd-1 do begin

   twrap = t mod s[i]
   srt = SORT(twrap)
   twrap = twrap[srt]
   fwrap = f[srt]
   plot,twrap,fwrap
   m[i] = TOTAL(ABS(fwrap[1:nt-1]-fwrap[0:nt-2]))
;   wait,1

endfor

plot,s,m
oplot,dr,[1,1]*mmin,linestyle=2



end
