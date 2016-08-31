;; get data from a certain year, starting 1st May

pro getyeardata,ts,fs,yr,tin,fin,tout,fout,AFTER=AFTER,nin=nin,nout=nout

  nyr = N_ELEMENTS(yr)
  nt = N_ELEMENTS(ts)
  jdcnv,yr,5,1,0,t0             ; 1st May of each year
  keep = FLTARR(nt)
  if KEYWORD_SET(after) then begin
     in = WHERE(ts gt t0[0])
     if in[0] ne -1 then keep[in] = 1
  endif else begin
     for i=0,nyr-1 do begin
        in = WHERE(ts gt t0[i] and ts le t0[i]+365.25)
        if in[0] ne -1 then keep[in] = 1
     endfor
  endelse
  in = WHERE(keep eq 1,nin,COMPLEMENT=out,NCOMPLEMENT=nout)
  tin = ts[in]
  fin = fs[in]
  tout = ts[out]
  fout = fs[out]

end
