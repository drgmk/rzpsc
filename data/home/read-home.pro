;; read in MVS data

readcol,'home/home-data.txt',t1,m1,t2,m2,t3,m3,t4,m4,t5,m5,t6,m6,t7,m7,t8,m8,format='f,a,f,a,f,a,f,a,f,a,f,a,f,a,f,a',delimiter='	'

t = [t1,t2,t3,t4,t5,t6,t7,t8]
mtmp = [m1,m2,m3,m4,m5,m6,m7,m8]
nt = N_ELEMENTS(t)
unc = INTARR(nt)
lim = INTARR(nt)

unc = WHERE( STREGEX(mtmp,':',/BOOL) )
;unc[tmp] = 1
lim = WHERE( STREGEX(mtmp,'L',/BOOL) )
;lim[tmp] = 1

m = FLTARR(nt)
for i=0,nt-1 do m[i] = FLOAT( STR_REPLACE(STR_REPLACE(mtmp[i],':',''),'L','') )

end
