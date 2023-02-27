;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; asd_to_instrument_v5 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro asd_to_instrument_v5, instrument_band_file, asd_spec_file, instrument_spec_file, integer_gain=integer_gain, micron=micron, header=header

  ; instrument_band_file is in the MODTRAN format
  ; /micron keyword means that the wavelength in the instrument_band-file is in microns and needs to be converted to nanometers
  ; /header: if set, reads file assuming a header with names, and writes the output file with header
  ;
  ;Read ASD file, 1st column wvl

  if keyword_set(header) eq 1 then begin
    temp=read_ascii(asd_spec_file, data_start=1)
  endif else begin
    temp=read_ascii(asd_spec_file)
  endelse

  temp=temp.(0)
  asd_wvl=temp(0,*)
  spec_asd=temp(1:*,*)
  sz=size(spec_asd)
  n_spec=sz(1)

  ;Read instrument band file
  temp=read_ascii(instrument_band_file, data_start=2)
  rsr_read=temp.(0)
  if keyword_set(micron) eq 1 then rsr_read(0,*)=rsr_read(0,*)*1000.
  temp=where(finite(rsr_read(0,*)) ne 1)
  n_bands=size(temp, /n_el) + 1
  rsr_startstop=intarr(2,n_bands)
  for i=0,n_bands-1 do begin
    if i eq 0 then rsr_startstop(0,i) = 0 else begin
      rsr_startstop(0,i) = temp(i-1)+1
    endelse
    if i eq n_bands-1 then rsr_startstop(1,i) = size(rsr_read(0,*), /n_el)-1 else begin
      rsr_startstop(1,i)=temp(i)-1
    endelse
  endfor
  curr_start=rsr_startstop(0,n_bands-1)
  curr_stop=rsr_startstop(1,n_bands-1)
  curr_rsr=rsr_read(*, curr_start:curr_stop)
  rsr_struct=create_struct(strcompress('Band_'+string(n_bands), /rem), curr_rsr)
  for i=1,n_bands-1 do begin
    curr_start=rsr_startstop(0,n_bands-1-i)
    curr_stop=rsr_startstop(1,n_bands-1-i)
    curr_rsr=rsr_read(*, curr_start:curr_stop)
    rsr_struct=create_struct(strcompress('Band_'+string(n_bands-i), /rem), curr_rsr, rsr_struct)
  endfor

  instrument_wvl=fltarr(1,n_bands)
  for j=0,n_bands-1 do begin
    curr_rsr=rsr_struct.(j)
    weights=interpol(curr_rsr(1,*), curr_rsr(0,*), asd_wvl)
    temp=where(asd_wvl lt min(curr_rsr(0,*)) or asd_wvl gt max(curr_rsr(0,*)))
    weights[temp]=0
    weights=weights/total(weights)
    instrument_wvl(j)=total(asd_wvl*weights)
  endfor

  spec_instrument=fltarr(n_spec, n_bands)

  for i=0,n_spec-1 do begin
    curr=spec_asd(i,*)
    for j=0,n_bands-1 do begin
      curr_rsr=rsr_struct.(j)
      weights=interpol(curr_rsr(1,*), curr_rsr(0,*), asd_wvl)
      temp=where(asd_wvl lt min(curr_rsr(0,*)) or asd_wvl gt max(curr_rsr(0,*)))
      weights[temp]=0
      weights=weights/total(weights)
      spec_instrument(i,j)=total(curr*weights)
    endfor
  endfor

  if keyword_set(integer_gain) then begin
    spec_instrument=fix(spec_instrument*integer_gain)
    fmt=strcompress('(f6.1,'+string(n_spec, format='(i)')+'(i,', /rem)+'"  "),(i))'
  endif else begin
    fmt=strcompress('(f6.1,'+string(n_spec, format='(i)')+'(f6.3,', /rem)+'"  "),(f6.3))'
  endelse

  out=[reform(instrument_wvl,1, n_bands), reform(spec_instrument, n_spec,n_bands)]

  ;Get headers if any
  if keyword_set(header) eq 1 then begin
    firstline=''
    get_lun, lun
    openr, lun, asd_spec_file
    readf, lun, firstline, format='(a)'
    close, lun
    free_lun, lun
    firstline=strsplit(firstline, string(byte(9)), /ex) 
    for i=0,n_spec do firstline(i)=' '+firstline(i)
  endif 

  get_lun, wlun
  openw, wlun,instrument_spec_file
  if keyword_set(header) then printf, wlun, firstline, format='('+string(n_spec+1, format='(i)')+'(a))'
  for i=0,n_bands-1 do printf, wlun, strcompress(string(out(*,i), format=fmt))
  close, wlun
  free_lun, wlun

  return
end
