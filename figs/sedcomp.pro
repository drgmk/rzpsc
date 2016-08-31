;; plot some different disk SEDs for comparison

c = 2.9979e14                    ; in um/s
pk = 1e-11                       ; put all at this lam.Flam at 2MASS J

thk = 5
csz = 1.4
!P.FONT=0

distinct_colors
set_plot,'ps'
device,file='spcomp.eps',/col,bits_per_pixel=8,/times,xsize=8,ysize=6,/inches,/enca

xr = [0.7,200]
yr = [5e-14,2e-11]
plot,[0],[0],xrange=xr,yrange=yr,/xl,/yl,/xst,/yst,xtitle=textoidl('wavelength (\mum)'),ytitle=textoidl('\lambdaF_\lambda (W m^{-2})'),xthick=thk,ythick=thk,charsize=csz,charthick=thk,position=[0.12,0.11,0.99,0.99]

;; Taurus median from D'Alessio et al 1999, this is a bit too cool to look very good
readcol,'taurus-median.txt',wav,loglamflam,hi,lo,format='f,f,f,f'
lamflam = 10^loglamflam
lamflamlo = 10^lo
lamflamhi = 10^hi
norm = pk / lamflam[5]
lamflam *= norm
lamflamlo *= norm
lamflamhi *= norm
;oplot,wav,lamflam,psym=1,thick=thk
;oplot,wav,lamflamlo,psym=1,thick=thk
;oplot,wav,lamflamhi,psym=1,thick=thk
;polyfill,[wav,REVERSE(wav)],[lamflamlo,REVERSE(lamflamhi)],color=200,noclip=0
;oplot,wav,lamflam,thick=thk

axis,yaxis=1,yrange=yr,/yst,ytickname=REPLICATE(' ',60),ythick=thk

;; rough photosphere
wav = logarrw(min=3,max=500,n=10)
oplot,wav,3.7e-11*wav^(-3),linestyle=2,thick=4

;; DM Tau
;; readcol,'dm-tau.txt',band,mag,format='a,f'
;; wav = getweffs(band)
;; fnu = mag
;; for i=0,N_ELEMENTS(mag)-1 do fnu[i] = mag2jy(mag[i],band[i])
;; tmp = tidy_cassis('cassis_yaaar_spcfw_27184640.fits',wav1,fnu1)
;; srt = SORT(wav1)
;; wav1 = wav1[srt]
;; fnu1 = fnu1[srt]
;; wav2 = logarrw(min=6,max=25,n=20)
;; fnu2 = interpol(fnu1,wav1,wav2)
;; wav = [wav,wav2,70.,100,160]
;; fnu = [fnu,fnu2,0.7,0.8,0.79]
;; wav1 = [wav,wav1,70.,100,160]
;; fnu1 = [fnu,fnu1,0.7,0.8,0.79]
;; lamflam = fnu / 1e26 * c / wav     ; convert to W/m^2
;; lamflam1 = fnu1 / 1e26 * c / wav1  ; convert to W/m^2
;; norm = pk / lamflam[2] / 2.
;; lamflam *= norm                    ; normalise
;; lamflam1 *= norm
;; oplot,wav,lamflam,psym=2,thick=10,color=FSC_COLOR('orange') ; plot W/m^2
;; srt = SORT(wav)
;; wav1 = wav(srt)
;; lamflam1 = lamflam(srt)
;; oplot,wav1,lamflam1,thick=thk,color=FSC_COLOR('orange')        ; plot W/m^2

;; GM Aur
readcol,'gm-aur.txt',band,fnumjy,format='a,f',comment='#'
wav = getweffs(band)
fnu = fnumjy/1e3
tmp = tidy_cassis('cassis_yaaar_spcfw_26141696.fits',wav1,fnu1)
srt = SORT(wav1)
wav1 = wav1[srt]
fnu1 = fnu1[srt]
wav2 = logarrw(min=6,max=25,n=20)
fnu2 = interpol(fnu1,wav1,wav2)
wav = [wav,wav2]
fnu = [fnu,fnu2]
srt = SORT(wav)
wav = wav[srt]
fnu = fnu[srt]
lamflam = fnu / 1e26 * c / wav     ; convert to W/m^2
norm = pk / lamflam[4] / 1.5
lamflam *= norm                    ; normalise
oplot,wav,lamflam,psym=2,thick=10,color=3 ; plot W/m^2
oplot,wav,lamflam,thick=thk,color=3        ; plot W/m^2

;; HD 166191
restore,'HIP89046.xdr'
tmp = WHERE(~STREGEX(phot.band,'_',/BOOLEAN))
lamflam = phot.flux[tmp] / 1e26 * c / getweffs(phot.band[tmp])       ; convert to W/m^2
tmp1 = WHERE(phot.band[tmp] eq '2MR1J')
lamflam *= pk / lamflam[tmp1[0]]
wav = getweffs(phot.band[tmp])
srt = SORT(wav)
wav = wav[srt]
lamflam = lamflam[srt]
oplot,wav,lamflam,thick=thk,color=2                ; plot W/m^2
oplot,wav,lamflam,symsize=1.2,thick=6,col=2,psym=4   ; plot W/m^2

;; HD 113766A
readcol,'hd113766a.txt',wav,fnu,format='f,f'
wav1 = logarrw(min=6,max=25,n=20)
tmp = WHERE(wav gt 5 and wav lt 36,complement=ok)
fnu1 = interpol(fnu[tmp],wav[tmp],wav1)
wav = [wav[ok],wav1]
fnu = [fnu[ok],fnu1]
srt = SORT(wav)
wav = wav[srt]
fnu = fnu[srt]
lamflam = fnu / 1e26 * c / wav                                     ; convert to W/m^2
lamflam *= pk / lamflam[2]                                         ; normalise
oplot,wav,lamflam,thick=thk,color=7                                ; plot W/m^2
oplot,wav,lamflam,psym=1,thick=8,symsize=1.2,col=7                 ; plot W/m^2

;; RZ Psc
restore,'~/astro/doc/rz-psc/sed/RZ-Psc.xdr'
tmp = WHERE(~STREGEX(phot.band,'_',/BOOLEAN))
lamflam = phot.flux[tmp] / 1e26 * c / getweffs(phot.band[tmp])       ; convert to W/m^2
tmp1 = WHERE(phot.band[tmp] eq '2MR2J')
lamflam *= pk / lamflam[tmp1[0]]
wav = getweffs(phot.band[tmp])
srt = SORT(wav)
wav = wav[srt]
lamflam = lamflam[srt]
ok = WHERE(phot.lim[srt[tmp]] eq 0,complement=lim)
oplot,wav[tmp[ok]],lamflam[tmp[ok]],thick=thk,color=1                    ; plot W/m^2
usersym,[-1,0,1,-1],[1,-1,1,1],fill=1,thick=6 ; downward pointing triangle for limits
oplot,wav[tmp[lim]],lamflam[tmp[lim]],symsize=1.4,thick=6,col=1,psym=8
oplotcircw,wav[tmp[ok]],lamflam[tmp[ok]],symsize=1.2,thick=6,/fill,col=1 ; plot W/m^2

legend,['RZ Psc','GM Aur','HD 166191','HD 113766A'],/bot,/lef,box=0,psym=[8,2,4,1],charsize=csz,charthick=thk,colors=[1,3,2,7],thick=[7,10,6,8],symsize=REPLICATE(1.1,4)
;polyfill,color=200,[0.415,0.48,0.48,0.415],[6.25,6.25,8.5,8.5]*1e-15*2.75

device,/close
set_plot,'x'

end
