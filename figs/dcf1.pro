PRO dcf,ta,ra,tb,rb,lag,cor,numf=numf, $
        erra=erra_in,errb=errb_in, $
        numpt=numpt,minpt=minpt, $
        minlag=minlag,maxlag=maxlag
;+
; NAME:
;       dcf
;
;
; PURPOSE:
;       compute the Edelson & Krolik discrete correlation function
;
;
; CATEGORY:
;       time series analysis tools
;
;
; CALLING SEQUENCE:
;       dcf,ta,ra,tb,rb,lag,cor,numf=numf, $
;           erra=erra_in,errb=errb_in, $
;           numpt=numpt,minpt=minpt, $
;           minlag=minlag,maxlag=maxlag
;
;
;
; INPUTS:
;       ta,ra: time and rate of the first light curve
;       tb,rb: time and rate of the second light curve
;
; KEYWORD PARAMETERS:
;       erra, errb: Arrays containing the measurement errors of
;             ra and rb, respectively. If given, this is taken into
;             account by using eq. (3) of Edelson & Krolik.
;             NOTE: There are lots of problems with interpreting the
;             DCF computed this way instead of using Edelson & Krolik,
;             eq. (2), and it is NOT recommended to give erra and errb
;             (see, e.g., White & Peterson, 1994, PASP, 106, 879)
;       minlag: minimum lag to consider
;       maxlag: maximum lag to consider
;       numf:   number of lag bins for which to compute the DCF
;       minpt:  minimum number of data points for a DCF value to be
;               computed (default: 10)
;
; OUTPUTS:
;       lag: lag value
;       cor: value of the DCF, or NAN if it cannot be determined
;            (<minpt )
;
;
; OPTIONAL OUTPUTS:
;       numpt:  number of data points entering each lag estimate
;
;
; RESTRICTIONS:
;     The lightcurves should be stationary, i.e., prewhitening by
;     subtracting off a low order polynomial from the data might be
;     required. See Welsh (1999, PASP 111, 1347)
;
; PROCEDURE:
;     Implementation of the formulae of Edelson & Krolik, ApJ 333,
;     646, 1988. You should read this paper as well as the papers
;     cited above before using this procedure.
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
; $Log: dcf.pro,v $
; Revision 1.3  2005/11/17 14:17:59  wilms
; force computation in double for means and standard deviations to make
; computations more robust.
;
; Revision 1.2  2005/11/17 12:36:47  wilms
; bug fix: the first cor bin could include data not within the
; respective lag bin
;
; Revision 1.1  2005/11/17 12:24:58  wilms
; initial release
;
; 
;-


    ;; characteristic error of time series A
    IF (n_elements(erra_in) EQ 0) THEN BEGIN 
        erra=0.
    ENDIF ELSE BEGIN 
        IF (n_elements(erra_in) EQ 1) THEN BEGIN 
            erra=erra_in[0]
        ENDIF ELSE BEGIN 
            erra=mean(erra_in,/double)
        ENDELSE 
    ENDELSE 

    ;; characteristic error of time series B
    IF (n_elements(errb_in) EQ 0) THEN BEGIN 
        errb=0.
    ENDIF ELSE BEGIN 
        IF (n_elements(errb_in) EQ 1) THEN BEGIN 
            errb=errb_in[0]
        ENDIF ELSE BEGIN 
            errb=mean(errb_in,/double)
        ENDELSE 
    ENDELSE 

    ;; minimum number of points to consider a CCF point valid
    IF (n_elements(minpt) EQ 0) THEN minpt=10
    IF (minpt LT 0) THEN BEGIN 
        message,'Minpt must be positive'
    ENDIF 

    amean=mean(ra,/double)
    asig=stddev(ra,/double)

    bmean=mean(rb,/double)
    bsig=stddev(rb,/double)

    tmp = WHERE(ra gt 0.5 and ra lt 1.2)
    meanclip,ra[tmp],amean,3,subs=ok ; alternative calculation of std dev and mean
;    amean = median(ra[tmp[ok]])
    asig = stddev(ra[tmp[ok]])
    tmp = WHERE(rb gt 0.5 and rb lt 1.2)
    meanclip,rb[tmp],bmean,3,subs=ok
;    bmean = median(rb[tmp[ok]])
    bsig = stddev(rb[tmp[ok]])
    
    na=n_elements(ra)
    nb=n_elements(rb)

    IF (na NE n_elements(ta)) THEN BEGIN
        message,'ta and ra must have same length'
    ENDIF 
    IF (nb NE n_elements(tb)) THEN BEGIN
        message,'tb and rb must have same length'
    ENDIF 

    ;; number of lags
    IF (n_elements(numf) EQ 0) THEN numf=long(min([na,nb])/10)

    udcf=dblarr(na,nb)
    dt=dblarr(na,nb)

    FOR i=0,na-1 DO BEGIN 
        udcf[i,*]=(ra[i]-amean)*(rb[*]-bmean)
        dt[i,*]=ta[i]-tb[*]
    ENDFOR 
    udcf=udcf/sqrt((asig^2.-erra^2.)*(bsig^2.-errb^2.))

    IF (n_elements(minlag) EQ 0) THEN minlag=min(dt)
    IF (n_elements(maxlag) EQ 0) THEN maxlag=max(dt)

    lag=minlag+(maxlag-minlag)*dindgen(numf+1)/numf
    ;; ...slow version
    ;; corr=dblarr(numf)
    ;; numpt=lonarr(numf)
    ;; FOR i=0,n_elements(lag)-2 DO BEGIN 
    ;;    ndx=where(dt GE lag[i] AND dt LT lag[i+1],num)
    ;;    IF (ndx[0] EQ -1) THEN BEGIN 
    ;;        corr[i]=!values.d_nan
    ;;        numpt[i]=0
    ;;    ENDIF ELSE BEGIN 
    ;;        corr[i]=mean(udcf[ndx])
    ;;        numpt[i]=num
    ;;    ENDELSE 
    ;;ENDFOR 

    ;; now sort all coefficients as a function of the lag
    ;; (makes finding the individual elements easier and
    ;; significantly speeds up the code)
    ndx=sort(dt)
    dt=dt[ndx]
    udcf=udcf[ndx]
    npt1=n_elements(dt)-1L
    cor=replicate(!values.d_nan,numf)
    numpt=lonarr(numf)
    ;; skip to first lag bin
    sta=0L
    WHILE(sta LT npt1 AND dt[sta] LT lag[0]) DO BEGIN 
        sta=sta+1
    ENDWHILE 
    sto=sta

    FOR i=0L,n_elements(lag)-2 DO BEGIN 
        ;; search next element in dt that is greater than the 
        ;; current lag
        WHILE (sto LT npt1 AND dt[sto] LT lag[i+1]) DO BEGIN 
            sto=sto+1L
        ENDWHILE 
        ;; compute the mean correlation if we have any udcf elements
        IF (sta LT sto) THEN BEGIN 
            IF (sto EQ npt1 AND dt[sto] LT lag[i+1]) THEN BEGIN 
                cor[i]=mean(udcf[sta:sto],/double)
;               cor[i]=MAX(udcf[sta:sto]) ;try this instead
                numpt[i]=sto-sta
            endif ELSE BEGIN 
                cor[i]=mean(udcf[sta:sto-1L],/double)
;                cor[i]=MAX(udcf[sta:sto-1L])
                numpt[i]=sto-sta-1L
            ENDELSE 
            sta=sto
        ENDIF
    ENDFOR 

    lag=(lag[0:numf-1]+lag[1:numf])/2.

    ndx=where(numpt LE minpt)
    IF (ndx[0] NE -1) THEN BEGIN 
        cor[ndx]=!values.d_nan
    ENDIF 
END 
