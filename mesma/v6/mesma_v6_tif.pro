pro mesma_v6_tif, case_number, input_file, $
; DEFAULT unmixing mode is to use hard constraints in LLS that the fractions cannot be lower/higher than min_frac/max_frac
    useband_file=useband_file, $
    force_lls=force_lls, $
    lls_sum_to_one=lls_sum_to_one, $    
    sum_to_one=sum_to_one, $
    normalize_fractions_before_rms=normalize_fractions_before_rms, $
    normalize_fractions_after_rms=normalize_fractions_after_rms, $
    guerschman=guerschman, $
;    use_explicit_guerschman_endmembers=use_explicit_guerschman_endmembers, $
    standardize_band=standardize_band, $
    equal_brightness=equal_brightness
;    write_spectral_image=write_spectral_image    

;directory='D:\'
;input_file=dialog_pickfile(path=directory, title='Please select all modis file', filter='*.tif', /multiple_files)
nfiles=n_elements(input_file)
for ss=0, nfiles-1 do begin
    temp_str = strsplit(input_file[ss],'.',/extract)
    ENVI_OPEN_FILE, input_file[ss], R_FID=fid, /INVISIBLE, /NO_INTERACTIVE_QUERY, /NO_REALIZE
    ENVI_File_Query, fid, DIMS=dims, ns=ns, nl=nl, nb=nb, data_type=data_type
    mapinfo=ENVI_GET_MAP_INFO(fid=fid)
    refl=fltarr(ns,nl,nb)
    for i=0,nb-1 do begin
      refl[*,*,i]=envi_get_data(fid=fid,dims=dims,pos=i)
    endfor
       
    ;nargs = N_PARAMS(0)
    ;IF nargs LT 1 THEN BEGIN
    ;        PRINT,'Usage: mesma_v6, in_file[, useband_file=useband_file, /force_lls, lls_sum_to_one=weighting_coefficient, /sum_to_one, /normalize_fractions_before_rms, /normalize_fractions_after_rms, /guerschman, /use_explicit_guerschman_endmembers, standardize_band=band_number, /equal_brightness]
    ;        print, 'DEFAULT unmixing mode is to use hard constraints that the fractions cannot be lower/higher than min_frac/max_frac
    ;        get_lun, lun
    ;        openw, lun, 'mesma_v6_in_template.txt'
    ;        printf, lun, ';MESMA SMA Input File (v6)
    ;        printf, lun, '; name of file containing Endmember Spectra (ASCII) (tab delimited, with header, and band first colum; if hyperspectral, first column must be wavelength in nm)
    ;        printf, lun, 'spectral_lib.txt
    ;        printf, lun, '; Scale value for endmember spectra(what value does 100% reflectance have?)
    ;        printf, lun, '1
    ;        printf, lun, ';Set this to 0 or the name of an instrument band file if a hyperspectral library is to be convolved to instrument wavebands
    ;        printf, lun, 'MODIS_Instrument_Filter_File_For_MODTRAN.flt
    ;        printf, lun, '; Model File - Spectra #s refer to number among the em spectra, not columns in file (since the library file must have a first column with band info, the second column is the first em
    ;        printf, lun, 'model_file.txt
    ;        printf, lun, ';Set this to 0 or the name of a MODIS NBAR .hdf  if that is the input
    ;        printf, lun, '0
    ;        printf, lun, '; name of file containing the reflectance image (16-bit BIL) ;Ignored if modis hdf
    ;        printf, lun, 'image.msb.bil ;;in_file
    ;        printf, lun, '; Swap (1 = yes, any other number = no) ;Ignored if modis hdf
    ;        printf, lun, '0
    ;        printf, lun, '; Image null value
    ;        printf, lun, '32767
    ;        printf, lun, '; Subset image? (1 = yes, any other number = no)
    ;        printf, lun, '0
    ;        printf, lun, ';If subsetting image, the inclusive starting X location (count from 0) ; ignored if subset image ne 1
    ;        printf, lun, '0
    ;        printf, lun, ';If subsetting image, the inclusive ending X location (count from 0) ; ignored if subset image ne 1
    ;        printf, lun, '0
    ;        printf, lun, ';If subsetting image, the inclusive starting Y location (count from 0) ; ignored if subset image ne 1
    ;        printf, lun, '0
    ;        printf, lun, ';If subsetting image, the inclusive ending Y location (count from 0) ; ignored if subset image ne 1
    ;        printf, lun, '0
    ;        printf, lun, '; X dimension of the BIL (samples) ;Ignored if modis hdf
    ;        printf, lun, '141
    ;        printf, lun, '; Y dimension of the BIL (rows) ;Ignored if modis hdf
    ;        printf, lun, '1
    ;        printf, lun, '; number of bands ;Ignored if modis hdf
    ;        printf, lun, '7
    ;        printf, lun, '; Scale value for reflectance image (what value does 100% reflectance have?) ;Ignored if modis hdf
    ;        printf, lun, '10000
    ;        printf, lun, '; Output file name (bil) - will write same byte order as reflectance file
    ;        printf, lun, 'mesma_out.bil
    ;        printf, lun, ';Maximum RMS Error (reflectance units) (needs to be higher for guerschmann transform)
    ;        printf, lun, '1
    ;        printf, lun, ';Maximum allowable fraction value (usually ~ 1.0)
    ;        printf, lun, '1.02
    ;        printf, lun, ';Minimum allowable fraction value (usually ~ 0.0)
    ;        printf, lun, '-0.02
    ;        printf, lun, ';Use an image as the source for one of the endmembers, model list should still have n endmembers (not n-1)
    ;        printf, lun, '0
    ;        printf, lun, ';Which endmember will be replaced with the value from a pixel from the image (1 = first endmember) Duplicate mixtures of em2 and em3 will be run (so watch out!)
    ;        printf, lun, '0
    ;        printf, lun, ';What is the path for the file that will be used as the endmember image (the source of the pixels)? Size must be same as subsetted size, reflectance range = [0,1]
    ;        printf, lun, 'EM_3_spectra_flt_MCD43A4.A2010129.h19v07_mesma_out_515_1230_6_9models_standardize=4.flt.bil
    ;        
    ;        
    ;    close, lun
    ;    Return
    ;endif
    
    ;Get system time
    t0=systime(1)
    x=0
    y=0
    
    ;initialize variables
    temp=''
    lib_file='C:\Users\zhang\Dropbox\mesma_AU\code\mesma\combined_library_isric_aster_thoralf_814x2151.txt'
    em_scale=1
    instrument_band_file='C:\Users\zhang\Dropbox\mesma_AU\code\mesma\MODIS_Instrument_Filter_File_For_MODTRAN.flt'
    model_file='C:\Users\zhang\Dropbox\mesma_AU\code\mesma\models_subset_1000_combined_library_isric_aster_thoralf_814x2151.txt'
    hdf_file='0'
    refl_file=input_file[ss]
    swap=0
    null=!VALUES.F_NAN
    subset_flag=0
    subset_xstart=0
    subset_xstop=0
    subset_ystart=0
    subset_ystop=0
    nx=ns  ;this is how many pixels you have
    ny=nl  ;this is how many rows you have in the pseudo stack image
    nbands_orig=nb ;number of bands
    refl_scale=10000.  
    out_file=temp_str[0]+'_mesma_'
    out_file_ave=temp_str[0]+'_mesma_ave_'
    max_rmse=1
    max_frac=1
    min_frac=0
    use_spectral_img=0
    use_spectral_img_em_number=0
    use_spectral_img_file=''
    
    print, 'Reading reflectance file...'

    ;Convert reflectance to floating point (0 to 1)
    refl=float(refl)/float(refl_scale)
    ;read library file
    print, 'Reading library file...'
    ;convert lib file to modis fit
    if instrument_band_file ne '0' then begin
      asd_to_instrument_v5, instrument_band_file, lib_file, 'temp_lib.txt', /header, /micron
      temp=read_ascii('temp_lib.txt', data_start=1)
    endif else begin  
      temp=read_ascii(lib_file, data_start=1)
    endelse
    temp=temp.(0)
    lib=temp(1:*,*)
    sz=size(lib)
    n_spectra=sz(1)
    spec_bands=sz(2)
    
    ;Convert library to floating point reflectance (0 to 1)
    lib=float(lib)/float(em_scale)
    
    if spec_bands ne nbands_orig then begin
      print, 'The number of bands for the endmembers is not equal to the number of bands'
      print, '  in the image.... Aborting.'
      stop
      return
    endif
    
    if keyword_set(guerschman) eq 1 then begin
    
      nbands_orig=84 ;number of bands in the guerschmann transform  
      new_refl=fltarr(nx,ny,nbands_orig)
      for i=0,nx-1 do for j=0,ny-1 do new_refl[i,j,*]=guerschman_transform_v5(reform(refl[i,j,*]))
      refl=new_refl
      if keyword_set(use_explicit_guerschman_endmembers) ne 1 then $
           lib=guerschman_transform_v5(lib)
    
    endif
    
    ;Keep a copy of the original library used
    orig_lib=lib
    ;Open the spectral image to be used as an EM, if necessary.  These will be from the orig_lib, not the subsetted or standardized libraries.
    if use_spectral_img eq 1 then begin
      readbil_flt, use_spectral_img_file, nx, ny, nbands_orig, em_refl_img
    endif
    
    if keyword_set(useband_file) then begin
     temp=read_ascii(useband_file)
     useband_list=temp.(0)
     useband=where(useband_list eq 1, nbands)
     if size(useband, /n_dim) eq 0 then begin
      print, 'Use Band File contains no bands to use.... returning.'
      return
     endif
     lib=lib[*,useband]
     refl=refl[*,*,useband]
     if use_spectral_img eq 1 then em_refl_img=em_refl_img[*,*,useband]
    endif else begin
     nbands=nbands_orig
    endelse
    
    if keyword_set(standardize_band) then begin
      print, 'Standardizing spectra to brightness of band (0 is 1st band):', standardize_band
      ;image
      if keyword_set(force_lls) then begin ;Force LLS doesn't work if there's a band with EXACTLY the same reflectance
        addme=(2*randomu(seed,nx,ny)-1)*1e-4
        addme2=(2*randomu(seed,nx,ny)-1)*1e-4
      endif else begin
        addme=fltarr(nx,ny)
        addme2=fltarr(nx,ny)
      endelse
      
      for i=0,nx-1 do begin
        for j=0, ny-1 do begin
          refl[i,j,*]=refl[i,j,*]/refl[i,j,standardize_band]+addme[i,j]
        endfor
      endfor
    
      if use_spectral_img eq 1 then begin
        for i=0,nx-1 do begin
          for j=0, ny-1 do begin
            em_refl_img[i,j,*]=em_refl_img[i,j,*]/em_refl_img[i,j,standardize_band]+addme2[i,j]
          endfor
        endfor
      endif
    
      ;library
      if keyword_set(force_lls) then begin ;Force LLS doesn't work if there's a band with EXACTLY the same reflectance
        addme=(2*randomu(seed,n_spectra)-1)*1e-4
      endif else begin
        addme=fltarr(n_spectra)
      endelse
    
      for i=0,n_spectra-1 do lib[i,*]=lib[i,*]/lib[i,standardize_band]+addme[i]
    endif
    
    if keyword_set(equal_brightness) then begin
      print, 'Setting all spectra to have equal brightness...'
      ;image
      for i=0,nx-1 do begin
        for j=0, ny-1 do begin
          refl(i,j,*)=refl(i,j,*)/sqrt(total(refl[i,j,*]^2.))
        endfor
      endfor
      
      ;Do this for the spectral image for endmembers
      if use_spectral_img eq 1 then begin
        for i=0,nx-1 do begin
          for j=0, ny-1 do begin
            em_refl_img(i,j,*)=em_refl_img(i,j,*)/sqrt(total(em_refl_img[i,j,*]^2.))
          endfor
        endfor
      endif
      
      ;library
      for i=0,nbands-1 do lib(i,*)=lib(i,*)/sqrt(total(lib[i,*]^2.))
    endif
    
    ;Add additional columns for linear least squares with sum to 1
    if keyword_set(force_lls) eq 1 and keyword_set(lls_sum_to_one) eq 1 then begin
      print, 'Adding lls_sum_to_one coefficient to library and reflectance image'
      print, strcompress('lls_sum_to_one coefficient set to: '+string(lls_sum_to_one))+' on a 100% reflectance = 1.0 scale'
      refl=[[[refl]],[[replicate(float(lls_sum_to_one), nx, ny, 1)]]]
      lib=[[lib],[replicate(float(lls_sum_to_one),n_spectra,1)]]
      nbands=nbands+1
    endif
    
    ;Read Model File
    Print, 'Reading model file....'
    temp=read_ascii(model_file)
    models=temp.(0)
    sz=size(models)
    if sz(0) eq 1 then begin
     n_em=sz(1)
     n_models=1
    endif else begin
     n_em=sz(1)
     n_models=sz(2)
    endelse
    
    ; Create output array
    ; out is the best selected model [model number, rmse, gv, npv, soil, shade] 
    out=fltarr(nx,ny, n_em+3)
    ; out_ave is the mean of all selected models [the number of models, mean rmse, mean gv, mean npv, mean soil, nothing]
    out_ave=fltarr(nx,ny,n_em+3)
    nx=nx-1 ; set these to nx-1 so that the for loop works
    ny=ny-1
    ;endelse
    ;Set default model to 9999
    out(*,*,0)=9999
    ;Set default RMSE to maximum RMS
    out(*,*,1)=max_rmse
    unmix=''
    ;Start loop over pixels
    for i=0, nx do begin
     if i mod 100 eq 0 then print, strcompress('Running line '+string(i+1)+' of '+string(nx+1)) ; should be column?
     for j=0, ny do begin
       x=i
       y=j 
      ;Write pixel reflectance to a variable
      pixel=float(refl(x,y,*))
      pixel=reform(pixel, 1, n_elements(pixel))
      temp=where(pixel eq null)
      if size(temp, /n_dim) eq 0 then begin
        ;Begin running models for the pixel
        for k=0,n_models-1 do begin
          abort = 1 ; start with abort = 1
      	  ;Create endmember array
      	  ems=fltarr(n_em, nbands)
      	  ;put in endmembers
      	  for l=0,n_em-1 do ems(l,*)=lib(models(l,k),*)
          ;Swap the use_spectral_img_em_number-th endmember  if using a spectral image as one of the endmembers.
          if use_spectral_img eq 1 then ems[use_spectral_img_em_number-1,*]=em_refl_img[x,y,*]
            ; START UNMIXING!      
            if keyword_set(force_lls) then begin  ; use linear least squares
               mod_grm_smdt_v5, ems, q, r
               fracs=invert(R)##transpose(Q)##pixel
            endif else begin 
               if keyword_set(sum_to_one) eq 1 then begin
                 fracs=la_least_square_equality(ems, replicate(1, n_em,1), pixel, [1], residual=residual)
               endif else begin ; The default is the constrained unmixing with hard fraction bounds [min and max fractions]
                 unmix = constrained_unmix_v6(ems, pixel)
                 fracs = reform(unmix(2:*),1,n_em)
               endelse  
            endelse
          	  ; Calculate RMSE    
              ;if the normalize_fractions_before_rms keyword is set,then do the normalization here - BEFORE calculating RMSE
              fracs=(keyword_set(normalize_fractions_before_rms) eq 0) *fracs + $
                    (keyword_set(normalize_fractions_before_rms) eq 1)*(fracs/total(fracs))     
              ; Calculate RMSE  
              rms_curr=sqrt(mean(total((float(ems)##fracs-float(pixel))^2.)))     
              ;if the normalize_fractions_after_rms keyword is set then do normalization here - AFTER calculting RMSE
              fracs=(keyword_set(normalize_fractions_after_rms) eq 0)*(fracs) + $
                    (keyword_set(normalize_fractions_after_rms) eq 1)*(fracs/total(fracs))  
          	  ;if this is the lowest RMS so far, put results in
          	  if (total(((fracs lt min_frac) + (fracs gt max_frac))) eq 0) or (finite(rms_curr) ne 1) then begin  ; Test for GT/LT max_frac/min_frac and NaNs
                ;calculate the model that has min rms for each endmember 
                if (rms_curr lt out(i,j,1)) then begin
                 ; Calculate shade endmember
                 shade=1.-total(fracs) 
            		 ; Insert model number
            		 out(i,j,0)=k
          	  	 ; Insert RMS
          		   out(i,j,1)=rms_curr
          		   ; Insert Fractions  
            		 out(i,j,2:n_em+1)=fracs
                 ; Calculate Shade
          	  	 out(i,j,n_em+2)=shade
          		  endif
          		  ; This does the endmember summing so that averages can be calculated
          		  out_ave[i,j,0]=out_ave[i,j,0]+1 ;counting variable
          		  out_ave[i,j,1]=out_ave[i,j,1]+rms_curr
          		  out_ave[i,j,2:n_em+1]=out_ave[i,j,2:n_em+1]+(fracs)
          	  endif
        endfor ; End loop over k
      endif ; End test for size(temp, /n_dim)
     endfor ; End loop over j
     time_taken=systime(1)-t0
     time_per_line=time_taken/(i+1)
     lines_remaining=nx-i+1
     time_remaining=round(lines_remaining*time_per_line/60.)
     temp="Completed"+string(i+1)+" of "+string(nx+1)+$
    	    " lines.  Estimated Time Remaining:"+string(time_remaining)+" minutes."
    
    endfor ; End loop over i
    print, 'Calculating average RMSE and endmember fractions for all good models (models within [min_frac, max_frac] and with no NANs)'
    for k=0, n_em+1 do out_ave[*,*,1+k]=out_ave[*,*,1+k]/out_ave[*,*,0]
    
    ;Report the total time for the analysis
    time_taken=round((systime(1)-t0)/60.)
    print, strcompress("Elapsed Time: "+ string(time_taken)+" minutes.")
    
    print, 'Write output image files (including average) in the same byteorder as the input file'
    envi_setup_head, data_type=4, interleave=0, fname=strcompress(out_file+string(case_number)), nb=6, ns=ns, nl=nl, BNAMES=['model number','RMSE','gv','npv','soil','shade']
    ENVI_WRITE_ENVI_FILE, out, interleave=0, nb=6, ns=ns, nl=nl, BNAMES=['model number','RMSE','gv','npv','soil','shade'], out_dt=4, map_info=mapinfo, out_name=strcompress(out_file+string(case_number)), /NO_OPEN
    envi_setup_head, data_type=4, interleave=0, fname=strcompress(out_file_ave+string(case_number)), nb=6, ns=ns, nl=nl, BNAMES=['model number','RMSE','gv','npv','soil','shade']
    ENVI_WRITE_ENVI_FILE, out_ave, interleave=0, nb=6, ns=ns, nl=nl, BNAMES=['model number','RMSE','gv','npv','soil','shade'], out_dt=4, map_info=mapinfo, out_name=strcompress(out_file_ave+string(case_number)), /NO_OPEN

;    print, 'Write output files (including average) in the same byteorder as the input file'
;    if swap eq 1 then begin
;      writebil_flt, out, strcompress(out_file+string(case_number)+'.bil', /rem), /swap
;      writebil_flt, out_ave, strcompress(out_file_ave+string(case_number)+'.bil', /rem), /swap
;    endif else begin
;      writebil_flt, out, strcompress(out_file+string(case_number)+'.bil', /rem)
;      writebil_flt, out_ave, strcompress(out_file_ave+string(case_number)+'.bil', /rem)
;    endelse
    
    print, 'Write output csv files (including average) in the same byteorder as the input file'
    write_csv, strcompress(out_file+string(case_number)+'.csv', /rem), out[*,*,0], out[*,*,1], out[*,*,2], out[*,*,3], out[*,*,4], out[*,*,5]
    write_csv, strcompress(out_file_ave+string(case_number)+'.csv',/rem), out_ave[*,*,0], out_ave[*,*,1], out_ave[*,*,2], out_ave[*,*,3], out_ave[*,*,4], out_ave[*,*,5]
    ENVI_FILE_MNG, id=fid, /remove 
endfor ; End mutliple images

return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;; do_write_spectral_image_v6 ;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro do_write_spectral_image_v6, mesma_out, ems_to_use, lib, models, out_file_name, strip_last_band=strip_last_band

print, 'Writing Endmember Spectral Images...'

n_ems_to_use=size(ems_to_use, /n_el)
;ems_to_use=ems_to_use-1

sz=size(mesma_out)
nx=sz[1]
ny=sz[2]
n_em=sz[3]-3

sz=size(lib)
nbands=sz(2)

if keyword_set(strip_last_band) then  use_lib=lib[*,0:nbands-2] else use_lib=lib
sz=size(use_lib)
nbands=sz(2)

  for k=0, n_ems_to_use-1  do begin
    print, strcompress('...EM:'+string(k+1))
    curr_out_file_name=strcompress('EM_'+string(ems_to_use[k])+'_spectra_flt_'+out_file_name, /rem)
    specimgout=fltarr(nx,ny,nbands)
    for i=0,nx-1 do begin
      for j=0,ny-1 do begin
        
        if mesma_out[i,j,0] ne 9999 then $
        specimgout[i,j,*]=use_lib[models[ems_to_use[k]-1, mesma_out[i,j,0]],*]
        
      endfor
    endfor
    
    writebil_flt, specimgout, curr_out_file_name
    
  endfor
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; readbil_16 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO READBIL_16, fname, cells, lines, bands,outarr, swap=swap
  nargs = N_PARAMS(0)
  IF nargs LT 1 THEN BEGIN
    PRINT,'Usage: READBIL_16, fname, cells, lines, bands, outarr, [/swap]'
    RETURN
  ENDIF
  ; Reads a 16-bit Band Interleaved by Line image and writes it to a
  ; integer array that's [cells, lines, bands] in size, named outarr
  ;  initialize outarr
  outarr=intarr(cells,lines,bands)
  i=long(0)   ; band counter
  j=long(0)   ;  y-counter- count number of lines
  k=long(0) ; total counter
  get_lun, lun
  openr, lun, fname
  if bands gt 1 then begin
    line=intarr(cells)
    repeat begin  ; repeat over lines
      repeat begin  ; repeat over bands for line j
        ; temp_line= (line[k] * (swap eq 0)) + (swap_endian(line[k])*(swap ne 0))
        temp_line=assoc(lun, fltarr(cells))
        outarr(*, j, i) = temp_line[i]
        i=i+1
        k=k+1
      endrep until i eq bands
      i=long(0)
      j=j+1
    endrep until j eq lines
  endif else begin
    img=assoc(lun, intarr(cells, lines))
    outarr=img[0]
  endelse
  free_lun, lun
  ; Perform byte swapping if swap is set
  if keyword_set(swap) ne 0 then byteorder, outarr
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; readbil_flt ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO READBIL_flt, fname, cells, lines, bands, outarr
  nargs = N_PARAMS(0)
  IF nargs LT 1 THEN BEGIN
    print, 'Usage:READBIL_flt, fname, cells, lines, bands, outarr'
    return
  endif
  ; Reads a floating point Band Interleaved by Line image and writes it to a
  ; FLT array that's [cells, lines, bands] in size, named outarr
  ;  initialize outarr
  outarr=fltarr(cells,lines,bands)
  i=long(0)   ; band counter
  j=long(0)   ;  y-counter- count number of lines
  k=long(0)   ; total counter
  get_lun, lun
  openr, lun, fname
  if bands gt 1 then begin
    line=assoc(lun, fltarr(cells))
    repeat begin  ; repeat over lines
      repeat begin  ; repeat over bands for line j
        temp_line=line[k]
        outarr(*, j, i) = temp_line
        i=i+1
        k=k+1
      endrep until i eq bands
      i=long(0)
      j=j+1
    endrep until j eq lines
  endif else begin
    img=assoc(lun, fltarr(cells, lines))
    outarr=img[0]
  endelse
  free_lun, lun
  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; writebil_flt ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO WRITEBIL_flt, img, fname
  nargs = N_PARAMS(0)
  IF nargs LT 1 THEN BEGIN
    print, 'Usage:WRITEBIL_flt, img, fname'
    return
  endif
  ;  This procedure writes a floating-point integer 3-d array to file
  ;  in BIL format
  ;
  ;  This procedure assumes the following file structure.
  ;  1st dimension: cells in a row
  ;  2nd dimension: lines in an image
  ;  3rd dimension: bands for each pixel
  sz=size(img)
  dims=sz(0)
  cells=sz(1)
  lines=sz(2)
  if dims eq 3 then bands=sz(3) else bands=1

  i=long(0)   ; band counter
  j=long(0)   ;  y-counter- count number of lines
  k=long(0) ; total counter
  get_lun, lun
  openw, lun, fname
  if bands ne 1 then begin
    line_out=assoc(lun, fltarr(cells))
    repeat begin  ; repeat over lines
      repeat begin  ; repeat over bands for line j
        line_out[k]=img(*,j,i)
        i=i+1
        k=k+1
      endrep until i eq bands
      i=long(0)
      j=j+1
    endrep until j eq lines
  endif else begin
    img_out=assoc(lun, fltarr(cells, lines))
    img_out[0]=img
  endelse
  free_lun, lun
  return
end

