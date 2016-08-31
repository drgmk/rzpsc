;; tidy AAVSO data

;readcol,'aavso-topcat.csv',jd,mag,unc,format='d,f,f',delimiter=','

readfmt,'aavso-topcat.dat','2X,d13,2X,f6,6X,a6,8X,a5,11X,a4,1X,a3,2X,a3,1X,a11,1x,a11,1X,a13',jd,mag,unc,hqunc,band,obs,com,comp1,comp2,chart
nl = N_ELEMENTS(jd)
unc1 = FLTARR(nl)
for i=0,nl-1 do begin
   if STRTRIM(unc[i],2) eq '""' then unc1[i] = !VALUES.F_NAN else unc1[i] = FLOAT(unc[i])
endfor

keep = INTARR(nl) + 1

;; only V band
tmp = WHERE(STRMID(band,0,1) ne 'V')
keep[tmp] = 0

obss = obs[UNIQ(obs,SORT(obs))]
nobss = N_ELEMENTS(obss)
nobs = INTARR(nobss)
for i=0,nobss-1 do begin
   tmp = WHERE(obs eq obss[i],n)
;   plot,jd[tmp],mag[tmp],psym=1,/yst
;   legend,[obss[i]],/top,/rig
   nobs[i] = n
   print,obss[i],n
   if n le 7 then keep[tmp] = 0
endfor


;; MTH tends to be a bit too faint, bump up by 0.3 mag (typical 12.0 to 11.7)
tmp = WHERE(obs eq 'MTH')
mag[tmp] -= 0.3

;; ignore measurements from POX, SPK, SXN (large scatter with few measurements)
tmp = WHERE(obs eq 'POX' or obs eq 'SPK' or obs eq 'SXN')
keep[tmp] = 0

ki = WHERE(keep eq 1,nk)
jdkeep = jd[ki]
magkeep = mag[ki]
obskeep = obs[ki]

;; output
openw,lun,'aavso-tidy.txt',/GET_LUN
printf,lun,'JD,mag,obs'
for i=0,nk-1 do printf,lun,STRING(jdkeep[i],format='(d)')+','+STRING(magkeep[i],format='(f)')+','+obskeep[i]
free_lun,lun

end
