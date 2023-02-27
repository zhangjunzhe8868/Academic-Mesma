pro mesma_v7, case_number, in_file, method, $ 
    useband_file=useband_file, $
    sum_to_one_weight=sum_to_one_weight, $
    normalize_fractions_before_rms=normalize_fractions_before_rms, $
    normalize_fractions_after_rms=normalize_fractions_after_rms, $
    guerschman=guerschman, $
    use_explicit_guerschman_endmembers=use_explicit_guerschman_endmembers, $
    standardize_band=standardize_band, $
    equal_brightness=equal_brightness, $
    write_spectral_image=write_spectral_image
    

;nargs = N_PARAMS(0)
;IF nargs LT 1 THEN BEGIN
;        PRINT,'Usage: mesma_v7, in_file[, method][, useband_file=useband_file, sum_to_one_weight=weighting_coefficient,
;        print,'   /normalize_fractions_before_rms, /normalize_fractions_after_rms, /guerschman, /use_explicit_guerschman_endmembers, 
;        print,'   standardize_band=band_number, /equal_brightness]
;        
;        print,'Methods:
;        print,'          1 = Modified Graham-Schmidt QR decomposition Linear Least Squares
;        print,'          2 = Hard constraints that fractions sum to one using generalized least squres with QR Factorization; IDL la_least_square_equality (SGGLSE from LAPACK users guide 3rd ed)
;        print,'          3 = Constraints on minimum and maximum fractions using sum of squared error minimization through partial 1st Derivatives (GRG algorithm supplied by Windward Technologies, Inc)
;        print,'          4 = Bounded value linear least squares (BVLS: bound constrained linear least-squares minimization; a generalization of NNLS that appeared in SOLVING LEAST SQUARES PROBLEMS, by Lawson and Hanson, Prentice-Hall, 1974.)
;        print,'          5 = Singular value decomposition (SVD) linear least squares
;        print, 'DEFAULT (if method parameter missing) unmixing mode is Method 1: Modified Graham-Schmidt Linear Least Squares
;        get_lun, lun
;        openw, lun, 'mesma_v7_in_template.txt'
;        printf, lun, ';MESMA SMA Input File (v7)
;        printf, lun, '; name of file containing Endmember Spectra (ASCII) (tab delimited, with header, and band first colum; if hyperspectral, first column must be wavelength in nm)
;        printf, lun, 'spectral_lib.txt
;        printf, lun, '; Scale value for endmember spectra(what value does 100% reflectance have?)
;        printf, lun, '100
;        printf, lun, ';Set this to 0 or the name of an instrument band file if a hyperspectral library is to be convolved to instrument wavebands
;        printf, lun, '0
;        printf, lun, '; Model File - Spectra #s refer to number among the em spectra, not columns in file (since the library file must have a first column with band info, the second column is the first em
;        printf, lun, 'model_file.txt
;        printf, lun, ';Set this to 0 or the name of a MODIS NBAR .hdf  if that is the input
;        printf, lun, '0
;        printf, lun, '; name of file containing the reflectance image (16-bit BIL) ;Ignored if modis hdf
;        printf, lun, 'image.msb.bil
;        printf, lun, '; Image null value
;        printf, lun, '32767
;        printf, lun, '; Subset image? (1 = yes, any other number = no)
;        printf, lun, '0
;        printf, lun, ';If subsetting image, the inclusive starting X location (count from 0) ; ignored if subset image ne 1
;        printf, lun, '200
;        printf, lun, ';If subsetting image, the inclusive ending X location (count from 0) ; ignored if subset image ne 1
;        printf, lun, '209
;        printf, lun, ';If subsetting image, the inclusive starting Y location (count from 0) ; ignored if subset image ne 1
;        printf, lun, '200
;        printf, lun, ';If subsetting image, the inclusive ending Y location (count from 0) ; ignored if subset image ne 1
;        printf, lun, '209
;        printf, lun, '; X dimension of the BIL (samples) ;Ignored if modis hdf
;        printf, lun, '1200
;        printf, lun, '; Y dimension of the BIL (rows) ;Ignored if modis hdf
;        printf, lun, '1200
;        printf, lun, '; number of bands ;Ignored if modis hdf
;        printf, lun, '7
;        printf, lun, '; Scale value for reflectance image (what value does 100% reflectance have?) ;Ignored if modis hdf
;        printf, lun, '10000
;        printf, lun, '; Output file name (bil) - will write same byte order as reflectance file
;        printf, lun, 'mesma_out.bil
;        printf, lun, ';Maximum RMS Error (reflectance units) (needs to be higher for guerschmann transform)
;        printf, lun, '0.5
;        printf, lun, ';Maximum allowable fraction value (usually ~ 1.0)
;        printf, lun, '1.02
;        printf, lun, ';Minimum allowable fraction value (usually ~ 0.0)
;        printf, lun, '-0.02
;        printf, lun, ';Use an image as the source for one of the endmembers, model list should still have n endmembers (not n-1)
;        printf, lun, '0
;        printf, lun, ';Which endmember will be replaced with the value from a pixel from the image (1 = first endmember) Duplicate mixtures of em2 and em3 will be run (so watch out!)
;        printf, lun, '3
;        printf, lun, ';What is the path for the file that will be used as the endmember image (the source of the pixels)? Size must be same as subsetted size, reflectance range = [0,1]
;        printf, lun, 'EM_3_spectra_flt_MCD43A4.A2010129.h19v07_mesma_out_515_1230_6_9models_standardize=4.flt.bil
;        
;        
;    close, lun
;    Return
;endif

;Get system time
t0=systime(1)

;;Set the keyword up so that if xy is set or not, code runs correctly
;;if keyword_set(xy) then begin
;; x=xy[0]
;; y=xy[1]
;;endif else begin
; x=0
; y=0
;;endelse

;initialize variables
temp=''
lib_file='C:\Users\zhang\Dropbox\mesma_AU\code\mesma\combined_library_isric_aster_thoralf_814x2151.txt'
em_scale=1
instrument_band_file='C:\Users\zhang\Dropbox\mesma_AU\code\mesma\MODIS_Instrument_Filter_File_For_MODTRAN.flt'
model_file='C:\Users\zhang\Dropbox\mesma_AU\code\mesma\models_subset_1000_combined_library_isric_aster_thoralf_814x2151.txt'
hdf_file=''
refl_file=in_file
swap=0
null=32767    ;!VALUES.F_NAN
subset_flag=0
subset_xstart=0
subset_xstop=0
subset_ystart=0
subset_ystop=0
nx=''  ;this is how many pixels you have
ny=''  ;this is how many rows you have in the pseudo stack image
nbands_orig=7 ;number of bands
refl_scale=10000.  ;for version 5 data, use 10000 and for version 6, use 1
out_file='d:\mesma_out_test_'
out_file_ave='d:\mesma_out_ave_test_'
max_rmse=1.
max_frac=1.
min_frac=0.
use_spectral_img=0
use_spectral_img_em_number=0
use_spectral_img_file=''

; open and read in_file

;get_lun, lun
;openr, lun, in_file
;print, 'Parsing input file...'
;
;readf, lun, temp, format='(a)'
;
;readf, lun, temp, format='(a)'
;readf, lun, lib_file, format='(a)'
;
;readf, lun, temp, format='(a)'
;readf, lun, em_scale, format='(i)'
;
;readf, lun, temp, format='(a)'
;readf, lun, instrument_band_file, format='(a)'
;
;readf, lun, temp, format='(a)'
;readf, lun, model_file, format='(a)'
;
;readf, lun, temp, format='(a)'
;readf, lun, hdf_file, format='(a)'
;
;readf, lun, temp, format='(a)'
;readf, lun, refl_file, format='(a)'
;
;readf, lun, temp, format='(a)'
;readf, lun, null, format='(i)'
;
;readf, lun, temp, format='(a)'
;readf, lun, subset_flag, format='(i)'
;
;readf, lun, temp, format='(a)'
;readf, lun, subset_xstart, format='(i)'
;
;readf, lun, temp, format='(a)'
;readf, lun, subset_xstop, format='(i)'
;
;readf, lun, temp, format='(a)'
;readf, lun, subset_ystart, format='(i)'
;
;readf, lun, temp, format='(a)'
;readf, lun, subset_ystop, format='(i)'
;
;readf, lun, temp, format='(a)'
;readf, lun, nx, format='(i)'
;
;readf, lun, temp, format='(a)'
;readf, lun, ny, format='(i)'
;
;readf, lun, temp, format='(a)'
;readf, lun, nbands_orig, format='(f)'
;
;readf, lun, temp, format='(a)'
;readf, lun, refl_scale, format='(i)'
;
;readf, lun, temp, format='(a)'
;readf, lun, out_file, format='(a)'
;
;readf, lun, temp, format='(a)'
;readf, lun, max_rmse, format='(f)'
;
;readf, lun, temp, format='(a)'
;readf, lun, max_frac, format='(f)'
;
;readf, lun, temp, format='(a)'
;readf, lun, min_frac, format='(f)'
;
;readf, lun, temp, format='(a)'
;readf, lun, use_spectral_img, format='(i)'
;
;readf, lun, temp, format='(a)'
;readf, lun, use_spectral_img_em_number, format='(i)'
;
;readf, lun, temp, format='(a)'
;readf, lun, use_spectral_img_file, format='(a)'
;
;
;close, lun
;free_lun, lun

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;DEFAULT UNMIXING IS LLS - METHOD 1;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
if keyword_set(method) eq 0 then method=1
if (floor(method) ne method) or method lt 1 or method gt 5 then begin
  print, 'Method does not exist.
  return
endif
if method eq 1 then print, 'Using default unmixing method: unconstained linear least squares.'

if keyword_set(use_explicit_guerschman_endmembers) eq 1 and keyword_set(guerschman) eq 0 then guerschman=1

;Read Refletance
Print, 'Reading reflectance file...
;If hdf_file not equal to 0 then create open the MODIS infile and get other info
if hdf_file ne '0' then begin
;  tmpltexists=file_test('modis_nbar_read_template.sav')
;  if tmpltexists eq 0 then begin
;    print, 'Select every MODIS NBAR band...
;    tmplt=hdf_browser(hdf_file)
;    save, tmplt, filename='modis_nbar_read_template.sav'
;  endif else begin
;    restore, 'modis_nbar_read_template.sav'
;  endelse
  
  hdfid = hdf_sd_start(in_file)
  hdf_sd_fileinfo, hdfid, nvars, ngatts
  refl=[]
  for i=0,nvars-1 do begin
    varid = hdf_sd_select(hdfid, i)
    hdf_sd_getdata, varid, data
    refl=[[[refl]],[[data]]]
  endfor

;  help,refl
;  print,refl[1000,1000,*]
;  refl=refl[1000:1009,1000:1009,*]
  
  sz=size(refl)
  nx=sz[1]
  ny=sz[2]
  nbands_orig=sz[3]
endif else begin
  ;Open reflectance file
  readbil_16, refl_file, nx, ny, nbands_orig, refl
endelse

;Subset image if needed
if subset_flag eq 1 then begin
  print, 'Subsetting...
  refl=refl(subset_xstart:subset_xstop,subset_ystart:subset_ystop,*)
  sz=size(refl)
  nx=sz[1]
  ny=sz[2]
endif

;Convert reflectance to floating point (0 to 1)
refl=float(refl)/float(refl_scale)
;read library file
print, 'Reading library file...'
if keyword_set(use_explicit_guerschman_endmembers) eq 1 then begin
  print, '.... will use explicit Guershman endmembers.'
endif else begin
  if instrument_band_file ne '0' then begin
    asd_to_instrument_v7, instrument_band_file, lib_file, 'instrument_lib.csv', /header, /micron
    temp=read_ascii('instrument_lib.csv', data_start=1)
  endif else begin  
    temp=read_ascii(lib_file, data_start=1)
  endelse
  temp=read_ascii('instrument_lib.csv', data_start=1)
  temp=temp.(0)
  lib=temp(1:*,*)
endelse

;Convert library to floating point reflectance (0 to 1)
if keyword_set(use_explicit_guerschman_endmembers) ne 1 then lib=float(lib)/float(em_scale)

;DO Guerschman Conversion
if keyword_set(guerschman) eq 1 then begin
  nbands_orig=84 ;number of bands in the guerschmann transform  
  new_refl=fltarr(nx,ny,nbands_orig)
  for i=0,nx-1 do for j=0,ny-1 do new_refl[i,j,*]=guerschman_transform_v7(reform(refl[i,j,*]))
  refl=new_refl
  
  if keyword_set(use_explicit_guerschman_endmembers) eq 1 then begin
    lib=guerschman_transform_v7(/return_explicit)   
  endif else begin
    lib=guerschman_transform_v7(lib)
  endelse
endif

;Check Size of Library
sz=size(lib)
n_spectra=sz(1)
spec_bands=sz(2)
if spec_bands ne nbands_orig then begin
  refl=refl[*,*,7:13]
  nbands_orig=7
endif

if spec_bands ne nbands_orig then begin
  print, 'The number of bands for the endmembers is not equal to the number of bands
  print, '  in the image.... Aborting.
  return
endif

;Keep a copy of the original library used
orig_lib=lib

;Open the spectral image to be used as an EM, if necessary. These will be from the orig_lib, not the subsetted or standardized libraries.
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

  if keyword_set(force_lls) then begin ;Force LLS doesn't work if there's a band with EXACTLY the same reflectance
    addme=(2*randomu(seed,nx,ny)-1)*1e-4
    addme2=(2*randomu(seed,nx,ny)-1)*1e-4
  endif else begin
    addme=fltarr(nx,ny)
    addme2=fltarr(nx,ny)
  endelse

  ;Actually none of the linear unmixing likely work if they're all they same - so add a little random for all
  addme=(2*randomu(seed,nx,ny)-1)*1e-4
  addme2=(2*randomu(seed,nx,ny)-1)*1e-4  
  
  ;Image
  for i=0,nx-1 do begin
    for j=0, ny-1 do begin
      refl[i,j,*]=refl[i,j,*]/refl[i,j,standardize_band] + addme[i,j]
    endfor
  endfor

  ;Spectral Image
  if use_spectral_img eq 1 then begin
    for i=0,nx-1 do begin
      for j=0, ny-1 do begin
        em_refl_img[i,j,*]=em_refl_img[i,j,*]/em_refl_img[i,j,standardize_band] + addme2[i,j]
      endfor
    endfor
  endif

  ;Library
  addme3=(2*randomu(seed,n_spectra)-1)*1e-4
  for i=0,n_spectra-1 do lib[i,*]=lib[i,*]/lib[i,standardize_band] + addme3[i]
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

;Add additional columns for linear least squares approches so that the fractions sum to 1
if keyword_set(sum_to_one_weight) eq 1 then begin
  print, 'Adding sum_to_one_weight coefficient to library and reflectance image
  print, strcompress('sum_to_one_weight coefficient set to: '+string(sum_to_one_weight))+' on a 100% reflectance = 1.0 scale'
  refl=[[[refl]],[[replicate(float(sum_to_one_weight), nx, ny, 1)]]]
  if use_spectral_img eq 1 then em_refl_img=[[[em_refl_img]],[[replicate(float(sum_to_one_weight), nx, ny, 1)]]]
  lib=[[lib],[replicate(float(sum_to_one_weight),n_spectra,1)]]
  nbands=nbands+1
endif
  
;Read Model File
if keyword_set(use_explicit_guerschman_endmembers) eq 1 then begin
  n_em=3
  n_models=1
  models=reform([0,1,2], 3, 1)
endif else begin
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
endelse 

;Create output array
out=fltarr(nx,ny,n_em+3)
out_ave=fltarr(nx,ny,n_em+3)
;Set default model to 9999
out(*,*,0)=9999
;Set default RMSE to maximum RMS
out(*,*,1)=max_rmse
unmix=''

;Start loop over pixels
for i=0, nx-1 do begin
 if i mod 100 eq 0 then print, strcompress('Running line '+string(i+1)+' of '+string(nx+1)) ; should be column?
 for j=0, ny-1 do begin
  x=i
  y=j 
  ;Write pixel reflectance to a variable
  pixel=float(refl(x,y,*))
  pixel=reform(pixel, 1, n_elements(pixel))
  temp=where(pixel eq null)
  if size(temp, /n_dim) eq 0 then begin    
    ;Begin running models for the pixel
    for k=0, n_models-1 do begin
      abort = 1 ; start with abort = 1
  	  ;Create endmember array
  	  ems=fltarr(n_em, nbands)
  	  ;put in endmembers
  	  for l=0,n_em-1 do ems(l,*)=lib(models(l,k),*)
      ;Swap the use_spectral_img_em_number-th endmember  if using a spectral image as one of the endmembers.
      if use_spectral_img eq 1 then ems[use_spectral_img_em_number-1,*]=em_refl_img[x,y,*]
      ;;;;;;;;;;;;;;;;;;;;;;;
      ;;; START UNMIXING! ;;; 
      ;;;;;;;;;;;;;;;;;;;;;;;
 
      ;Method 1 = MGS LLS     
      if method eq 1 then begin
        mod_grm_smdt_v5, ems, q, r
        fracs=invert(R)##transpose(Q)##pixel
      endif

      ;Method 2 = Hard constraint to one (la_least_Square_equality)
      if method eq 2 then begin
        fracs=la_least_square_equality(ems, replicate(1, n_em,1), pixel, [1], residual=residual)
      endif

      ;Method 3 = Constraints on endmember ---?
      if method eq 3 then begin
        unmix=constrained_unmix_v6(ems, pixel, min_frac, max_frac)
        fracs=reform(unmix(2:*),1,n_em)
      endif

      ;Method 4 = BVLS
      if method eq 4 then begin
        bnd = [replicate(float(min_frac), 1, 3), replicate(float(max_frac), 1, 3)]
        a = transpose(ems) ; must b row matrix (e.g. 84 x 3)
        b = [reform(pixel,n_elements(pixel),1)]
        if keyword_set(sum_to_one_weight) eq 0 then begin
           ;BVLS - unweighted
           bvls, a,b,bnd, fracs
        endif else begin
          ;BVLS - Weighted
           bvls, a, b, bnd, fracs, RNORM=rnorm
        endelse
      endif

      ;Method 5 = SVD
      if method eq 5 then begin
        SVDC, ems, W, U, V, /doubl
        inv_em=dblarr(nbands,n_em)
        for s=0,n_em-1 do inv_em=inv_em + (1/w[s])*(v[s,*]##transpose(u[s,*]))
        fracs=inv_em##pixel
      endif
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
        if (rms_curr lt out(x,y,1)) then begin 
          ; Calculate shade endmember
          shade=1.-total(fracs) 
      		; Insert model number
      		out(x,y,0)=k
  	      ; Insert RMS
    		  out(x,y,1)=rms_curr 
    		  ; Insert Fractions  
      		out(x,y,2:n_em+1)=fracs        
          ; Calculate Shade
    	  	out(x,y,n_em+2)=shade 		 
    		endif 		 
     		; This does the endmember summing so that averages can be calculated
    	  out_ave[x,y,0]=out_ave[x,y,0]+1 ;counting variable
    		out_ave[x,y,1]=out_ave[x,y,1]+rms_curr
  		  out_ave[x,y,2:n_em+1]=out_ave[x,y,2:n_em+1]+(fracs)  		 
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
print, 'Write output files (including average) in the same byteorder as the input file

envi_setup_head, data_type=4, interleave=0, fname=strcompress(out_file+string(case_number)), nb=6, ns=nx, nl=ny, BNAMES=['model number','RMSE','gv','npv','soil','shade']
ENVI_WRITE_ENVI_FILE, out, interleave=0, nb=6, ns=nx, nl=ny, BNAMES=['model number','RMSE','gv','npv','soil','shade'], out_dt=4, map_info=mapinfo, out_name=strcompress(out_file+string(case_number)), /NO_OPEN
envi_setup_head, data_type=4, interleave=0, fname=strcompress(out_file_ave+string(case_number)), nb=6, ns=nx, nl=ny, BNAMES=['model number','RMSE','gv','npv','soil','shade']
ENVI_WRITE_ENVI_FILE, out_ave, interleave=0, nb=6, ns=nx, nl=ny, BNAMES=['model number','RMSE','gv','npv','soil','shade'], out_dt=4, map_info=mapinfo, out_name=strcompress(out_file_ave+string(case_number)), /NO_OPEN

;writebil_flt, out, strcompress(out_file+'.bil', /rem)
;writebil_flt, out_ave, strcompress(out_file+'_averages'+'.bil', /rem)

print, 'Write output csv files (including average) in the same byteorder as the input file'
write_csv, strcompress(out_file+string(case_number)+'.csv', /rem), out[*,*,0], out[*,*,1], out[*,*,2], out[*,*,3], out[*,*,4], out[*,*,5]
write_csv, strcompress(out_file_ave+string(case_number)+'.csv',/rem), out_ave[*,*,0], out_ave[*,*,1], out_ave[*,*,2], out_ave[*,*,3], out_ave[*,*,4], out_ave[*,*,5]

HDF_SD_END, hdfid
;if (nx+1)*(ny+1) eq 1 then begin
;  print, strcompress("Unmixing method: "+string(unmix(0), format='(i)'))
;  print, ''
;  print, strcompress("Best model: "+ string(long(out[0])))
;  print, strcompress("RMSE: "+ string(out[1], format='(f10.3)'))
;  for k=0, n_em-1 do print, strcompress("Endmember #"+string(k)+" MESMA fraction: "+string(out[2+k], format='(f10.3)'))
;  print, strcompress("MESMA Shade fraction: "+string(out[n_em+2], format='(f10.3)'))
;  print, ''
;  print, strcompress("Average RMSE: "+ string(out_ave[1], format='(f10.3)'))
;  for k=0, n_em-1 do print, strcompress("Endmember #"+string(k)+" Average fraction: "+string(out_ave[2+k], format='(f10.3)'))
;  print, strcompress("Average Shade fraction: "+string(out_ave[n_em+2], format='(f10.3)'))
;endif
;
;if keyword_set(write_spectral_image) eq 1 then begin
;  if keyword_set(sum_to_one_weight) eq 1 then do_write_spectral_image_v6, out, indgen(n_em)+1, orig_lib, models, out_file, /strip_last_band $
;    else do_write_spectral_image_v6, out, indgen(n_em)+1, orig_lib, models, out_file
;endif 

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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; constrained_unmix_v6 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function sum_and_sse, fracs
  common shared, em_array, target
  sum=total(fracs)
  sse=total((float(em_array)##fracs-float(target))^2.) ;Sum of the square of errors (sse)
  return, [sum,sse]
end

function constrained_unmix_v6, em_array_in, target_in, min_frac, max_frac;, sum_to_one=sum_to_one
  ;em_array is a n_em x n_band array, target is a 1 x n_band.
  ;Output will be a 1xn_em array
  common shared, em_array, target

  em_array=em_array_in
  target=target_in

  n_bands=size(target, /n_elements)
  n_em=size(em_array, /n_elements)/n_bands
  fracs=fltarr(1,n_em)
  ;fracs(*)=1./float(n_em) *0.5 ; starting guess for fracs
  fracs=randomu(seed,1,n_em)     ; alternate starting guess for fracs
  fracs=fracs/total(fracs)*0.5   ; This also works

  ;print, reform(fracs)
;  frac_lower_bnd=fltarr(n_em)
;  frac_upper_bnd=frac_lower_bnd+1.
  frac_lower_bnd=fltarr(n_em)+min_frac
  frac_upper_bnd=fltarr(n_em)+max_frac  
  frac_bnd=[[frac_lower_bnd],[frac_upper_bnd]]
  sum_and_sse_bnd=[[float(keyword_set(sum_to_one)),0.],[1.,!values.f_infinity]]

  nobj=1 ; I think this is the index in what's return by sum_and_sse for what we want to minimize, here sse
  func='sum_and_sse'
  title='IDL constrained min report
  report='constrained_min_report.txt'
 
  constrained_min, fracs, frac_bnd, sum_and_sse_bnd, nobj, func, inform;, report=report, title=title;, epstop=1e-6
  ;print, sum_and_sse(fracs)
  ;print, reform(fracs)
  ;print,inform

  ;return, fracs
  return, [fix(inform),total((float(em_array)##fracs-float(target))^2.),reform(fracs)]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; mesma_unmixing_engine_v6 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mesma_unmixing_engine_v6, ems, pixel, n_em, max_frac, min_frac, sum_to_one=sum_to_one ;,force_mgs=force_mgs, 
  ;out(0) contains the final state flag
  ; value of 0 means constrained_unmix was used
  ; value of 1 means mgs was used
  ; value of 2 means constrained_unmix failed and mgs failed to give a result with fracs within [0,1]
  ;out(1) contains the sum of squared errors (sse), -1 means no fit
  ;out(2:*) contains the fracsions

  out=fltarr(2+n_em)
  count=0
;  if keyword_set(force_mgs) then begin
;    mod_grm_smdt_v5, ems, q, r
;    fracs=invert(R)##transpose(Q)##pixel
;    sse=total((float(ems)##fracs-float(pixel))^2.)
;    temp=where((fracs lt min_frac) or (fracs gt max_frac))
;    if size(temp, /n_dimensions) ne 0 or total(fracs) gt max_frac then begin
;      out(0)=3   ; indicates that both constrained_min and mgs failed
;      fracs(*)=0
;      sse=-1
;    endif
;  endif else begin
    if keyword_set(sum_to_one) eq 1 then begin
      fracs=la_least_square_equality(ems, replicate(1, n_em,1), pixel, [1], residual=residual)
      sse=total((float(ems)##fracs-float(pixel))^2.)
    endif else begin
;      repeat begin
        unmix=constrained_unmix_v6(ems, pixel)
;        count=count+1
;      endrep until unmix(0) ne 6 or count gt 10
      fracs=reform(unmix(2:*),1,n_em)
      sse=total((float(ems)##fracs-float(pixel))^2.)
      if unmix(0) eq 6 then begin
        out(0)=1    ; indicates that mgs was needed for unmixing
        mod_grm_smdt_v5, ems, q, r
        fracs=invert(R)##transpose(Q)##pixel
        sse=total((float(ems)##fracs-float(pixel))^2.)
        temp=where((fracs lt min_frac) or (fracs gt max_frac))
        if size(temp, /n_dimensions) ne 0 or total(fracs) gt max_frac then begin
          out(0)=3   ; indicates that both constrained_min and mgs failed
          fracs(*)=0
          sse=-1
        endif
      endif
    endelse
;  endelse

  out(1)=sse
  out(2:*)=fracs

  return, out
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; asd_to_instrument_v7 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro asd_to_instrument_v7, instrument_band_file, asd_spec_file, instrument_spec_file, integer_gain=integer_gain, micron=micron, header=header

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
    fmt=strcompress('(f6.1,",",'+string(n_spec, format='(i)')+'(i,', /rem)+'","),(i))'
  endif else begin
    fmt=strcompress('(f6.1,",",'+string(n_spec, format='(i)')+'(f6.3,', /rem)+'","),(f6.3))'
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
    firstline=strsplit(firstline, string(byte(9)), /ex) ; string(byte(9)) is a tab
    ; sz=size(out)
    ; out2=strarr(sz(1), sz(2)+1)
    ; out2(*,1:*)=out
    for i=0,n_spec do firstline(i)=' '+firstline(i)
    ; out2(*,0)=firstline
  endif ;else begin
  ; out2=out
  ;endelse

  get_lun, wlun
  openw, wlun,instrument_spec_file
  if keyword_set(header) then printf, wlun, firstline, format=strcompress('('+string(n_spec+1, format='(a)')+'(a,","),a )')
  for i=0,n_bands-1 do printf, wlun, strcompress(string(out(*,i), format=fmt))
  close, wlun
  free_lun, wlun

return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; guerschman_transform_v7 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function guerschman_transform_v7, mcd43, return_explicit=return_explicit
  ;This function converts either a spectrum or an 3-dimensional spectral image
  ; (with band as the 3rd dimension, bsq) to the transform given by guerschman et al 2015, RSE

  if keyword_set(return_explicit) eq 0 then begin

  ;b1=float(mcd43[*,*,0])
  ;b2=float(mcd43[*,*,1])
  ;b3=float(mcd43[*,*,2])
  ;b4=float(mcd43[*,*,3])
  ;b5=float(mcd43[*,*,4])
  ;b6=float(mcd43[*,*,5])
  ;b7=float(mcd43[*,*,6])

  sz=size(mcd43)

  if sz(0) eq 1 then begin
    nx=1
    ny=1
    b1=float(mcd43[0])
    b2=float(mcd43[1])
    b3=float(mcd43[2])
    b4=float(mcd43[3])
    b5=float(mcd43[4])
    b6=float(mcd43[5])
    b7=float(mcd43[6])
  endif

  if sz(0) eq 2 then begin
    nx=sz[1]
    ny=1
    b1=float(mcd43[*,0])
    b2=float(mcd43[*,1])
    b3=float(mcd43[*,2])
    b4=float(mcd43[*,3])
    b5=float(mcd43[*,4])
    b6=float(mcd43[*,5])
    b7=float(mcd43[*,6])
  endif

  if sz(0) eq 3 then begin
    nx=sz[1]
    ny=sz[2]
    b1=float(mcd43[*,*,0])
    b2=float(mcd43[*,*,1])
    b3=float(mcd43[*,*,2])
    b4=float(mcd43[*,*,3])
    b5=float(mcd43[*,*,4])
    b6=float(mcd43[*,*,5])
    b7=float(mcd43[*,*,6])
  endif

  transform=fltarr(nx,ny,84)

  transform[*,*,0]=B1
  transform[*,*,1]=B2
  transform[*,*,2]=B3
  transform[*,*,3]=B4
  transform[*,*,4]=B5
  transform[*,*,5]=B6
  transform[*,*,6]=B7
  transform[*,*,7]=alog(B1)
  transform[*,*,8]=alog(B2)
  transform[*,*,9]=alog(B3)
  transform[*,*,10]=alog(B4)
  transform[*,*,11]=alog(B5)
  transform[*,*,12]=alog(B6)
  transform[*,*,13]=alog(B7)
  transform[*,*,14]=alog(B1)*B1
  transform[*,*,15]=alog(B2)*B2
  transform[*,*,16]=alog(B3)*B3
  transform[*,*,17]=alog(B4)*B4
  transform[*,*,18]=alog(B5)*B5
  transform[*,*,19]=alog(B6)*B6
  transform[*,*,20]=alog(B7)*B7
  transform[*,*,21]=B1*B2
  transform[*,*,22]=B1*B3
  transform[*,*,23]=B1*B4
  transform[*,*,24]=B1*B5
  transform[*,*,25]=B1*B6
  transform[*,*,26]=B1*B7
  transform[*,*,27]=B2*B3
  transform[*,*,28]=B2*B4
  transform[*,*,29]=B2*B5
  transform[*,*,30]=B2*B6
  transform[*,*,31]=B2*B7
  transform[*,*,32]=B3*B4
  transform[*,*,33]=B3*B5
  transform[*,*,34]=B3*B6
  transform[*,*,35]=B3*B7
  transform[*,*,36]=B4*B5
  transform[*,*,37]=B4*B6
  transform[*,*,38]=B4*B7
  transform[*,*,39]=B5*B6
  transform[*,*,40]=B5*B7
  transform[*,*,41]=B6*B7
  transform[*,*,42]=alog(B1)*alog(B2)
  transform[*,*,43]=alog(B1)*alog(B3)
  transform[*,*,44]=alog(B1)*alog(B4)
  transform[*,*,45]=alog(B1)*alog(B5)
  transform[*,*,46]=alog(B1)*alog(B6)
  transform[*,*,47]=alog(B1)*alog(B7)
  transform[*,*,48]=alog(B2)*alog(B3)
  transform[*,*,49]=alog(B2)*alog(B4)
  transform[*,*,50]=alog(B2)*alog(B5)
  transform[*,*,51]=alog(B2)*alog(B6)
  transform[*,*,52]=alog(B2)*alog(B7)
  transform[*,*,53]=alog(B3)*alog(B4)
  transform[*,*,54]=alog(B3)*alog(B5)
  transform[*,*,55]=alog(B3)*alog(B6)
  transform[*,*,56]=alog(B3)*alog(B7)
  transform[*,*,57]=alog(B4)*alog(B5)
  transform[*,*,58]=alog(B4)*alog(B6)
  transform[*,*,59]=alog(B4)*alog(B7)
  transform[*,*,60]=alog(B5)*alog(B6)
  transform[*,*,61]=alog(B5)*alog(B7)
  transform[*,*,62]=alog(B6)*alog(B7)
  transform[*,*,63]=((B2-B1)/(B2+B1))
  transform[*,*,64]=((B3-B1)/(B3+B1))
  transform[*,*,65]=((B4-B1)/(B4+B1))
  transform[*,*,66]=((B5-B1)/(B5+B1))
  transform[*,*,67]=((B6-B1)/(B6+B1))
  transform[*,*,68]=((B7-B1)/(B7+B1))
  transform[*,*,69]=((B3-B2)/(B3+B2))
  transform[*,*,70]=((B4-B2)/(B4+B2))
  transform[*,*,71]=((B5-B2)/(B5+B2))
  transform[*,*,72]=((B6-B2)/(B6+B2))
  transform[*,*,73]=((B7-B2)/(B7+B2))
  transform[*,*,74]=((B4-B3)/(B4+B3))
  transform[*,*,75]=((B5-B3)/(B5+B3))
  transform[*,*,76]=((B6-B3)/(B6+B3))
  transform[*,*,77]=((B7-B3)/(B7+B3))
  transform[*,*,78]=((B5-B4)/(B5+B4))
  transform[*,*,79]=((B6-B4)/(B6+B4))
  transform[*,*,80]=((B7-B4)/(B7+B4))
  transform[*,*,81]=((B6-B5)/(B6+B5))
  transform[*,*,82]=((B7-B5)/(B7+B5))
  transform[*,*,83]=((B7-B6)/(B7+B6))

  case sz[0] of
    1: out=reform(transform)
    2: out=reform(transform)
    3: out=transform
  endcase

  endif else begin
    out=[$
       [0.36820, 0.42852, 0.43367],$
       [0.47454, 0.44349, 0.43271],$
       [0.22305, 0.21879, 0.20737],$
       [0.30857, 0.32828, 0.32776],$
       [0.50152, 0.48021, 0.44110],$
       [0.31322, 0.39141, 0.32365],$
       [0.22137, 0.23207, 0.33113],$
       [-0.12693, -0.06231, -0.07816],$
       [-0.16635, -0.30239, -0.29336],$
       [-0.19812, -0.28319, -0.35424],$
       [-0.34285, -0.32747, -0.29973],$
       [-0.20843, -0.27223, -0.30794],$
       [-0.32392, -0.30179, -0.39755],$
       [-0.24737, -0.33261, -0.25633],$
       [-0.43749, -0.43443, -0.39807],$
       [-0.28852, -0.28049, -0.26584],$
       [-0.38044, -0.37602, -0.35114],$
       [-0.36823, -0.37852, -0.36705],$
       [-0.25861, -0.25817, -0.25428],$
       [-0.42155, -0.35476, -0.34013],$
       [-0.57556, -0.50045, -0.42295],$
       [0.15701, 0.17908, 0.18490],$
       [0.06274, 0.07019, 0.07285],$
       [0.09088, 0.11098, 0.11895],$
       [0.13343, 0.17567, 0.17961],$
       [0.05773, 0.12955, 0.13801],$
       [0.01405, 0.06526, 0.11622],$
       [0.10327, 0.09705, 0.09425],$
       [0.14150, 0.14525, 0.14810],$
       [0.23646, 0.22149, 0.21179],$
       [0.16604, 0.18842, 0.17554],$
       [0.10612, 0.11357, 0.15255],$
       [0.05207, 0.05364, 0.05393],$
       [0.10719, 0.10370, 0.09720],$
       [0.07812, 0.08745, 0.08097],$
       [0.05434, 0.05965, 0.07223],$
       [0.14064, 0.15202, 0.15123],$
       [0.09720, 0.12829, 0.12790],$
       [0.06463, 0.08557, 0.11408],$
       [0.13803, 0.19331, 0.16413],$
       [0.05192, 0.08677, 0.13355],$
       [-0.08493, 0.00065, 0.05643],$
       [-0.04314, -0.05696, -0.09486],$
       [0.02327, 0.06850, 0.10378],$
       [0.46582, 0.48597, 0.43940],$
       [-0.31506, -0.25992, -0.29827],$
       [-0.30482, -0.21392, -0.33182],$
       [0.00736, -0.13326, -0.02477],$
       [0.06638, 0.03111, 0.04326],$
       [0.06467, 0.15029, 0.16924],$
       [0.08402, 0.04567, 0.00446],$
       [0.00363, -0.03340, 0.09239],$
       [-0.11916, -0.18609, -0.23697],$
       [-0.36308, -0.41304, -0.42504],$
       [0.09809, -0.00275, -0.01553],$
       [0.05874, 0.01206, -0.04059],$
       [0.08658, 0.13523, 0.09648],$
       [-0.11104, -0.10148, -0.07190],$
       [-0.10812, -0.09228, -0.05182],$
       [-0.07056, -0.02000, 0.03352],$
       [0.16846, 0.23203, 0.28153],$
       [0.12738, 0.14114, 0.13712],$
       [0.01632, 0.00232, -0.03020],$
       [-0.10766, -0.23018, -0.20297],$
       [-0.09555, -0.17943, -0.20902],$
       [-0.13120, -0.15798, -0.13758],$
       [-0.20369, -0.25596, -0.26286],$
       [-0.24226, -0.25595, -0.30518],$
       [-0.13133, -0.20186, -0.15313],$
       [-0.18292, -0.12126, -0.15061],$
       [-0.07674, 0.00820, 0.00263],$
       [-0.02534, 0.00701, -0.01539],$
       [-0.03543, 0.03055, -0.03298],$
       [0.11534, 0.10277, 0.12895],$
       [-0.04992, -0.00727, 0.03656],$
       [0.12227, 0.11448, 0.11398],$
       [0.06264, 0.09553, 0.06187],$
       [0.07373, 0.10524, 0.13783],$
       [-0.01237, -0.04502, -0.06621],$
       [-0.04390, -0.04394, -0.09199],$
       [0.02114, -0.00566, 0.01459],$
       [-0.03908, -0.00154, -0.04095],$
       [0.06730, 0.00639, 0.07500],$
       [0.06942, -0.00498, 0.10273]]
  endelse
  return, out
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; mod_grm_smdt_v5 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO mod_grm_smdt_v5, A, Q, R
  nargs = N_PARAMS(0)
  IF nargs LT 1 THEN BEGIN
    print, 'Usage: mod_grm_smdt, A, Q, R
    print, '  A is a matrix comprised of N column vectors with M variables
    print, '   and can be interpreted as a matrix of endmembers with columns
    print, '   corresponding to individual endmember spectra
    print, '  Q is an orthonormal basis set of Q which is to be used in
    print, '   linear least squares unmixing of the equation: Ax=b
    print, '     where x is a column vector giving proportions of each endmember
    print, '     in A to match the target vector (spectrum), b, also a column
    print, '     vector
    print, '  R is an upper triangular square matrix such that A=QR
    print, '
    print, '  The linear least squares solution of this system is given by:
    print, '    x-bar = invert(R)##transpose(Q)##b
    return
  endif
  ;  This procedure calculates the orthonormal basis set for A which is called
  ;  Q as well as the factoralization matrix R, such that A=QR
  ;  The method used is the modified graham-schmidt algorithm described
  ;   in Golub and Van Loan , "Matrix Computations", Johns Hopkins Univ. Press,
  ;   Baltimore MD, 1983, p. 152.

  ;  The basic graham-schmidt orthogonalization is given by the equations:
  ;  1)  q(1) (first column of Q) = a(1)
  ;  2)  q( j ) = a( j ) - Sum from i=1 to j-1 of  [qT( i )a(j) ) q (i) ]
  ;           where qT( i ) is the transpose of q( i )
  ;  Last step, normalize each of the vectors (columns) in q by dividing each element
  ;   of each column by the magnitude of the vector  = the sqrt of the sum of the squares
  ;   of the column entries
  temp=size(A)
  cols=temp[1]  ; columns = number of endmembers
  rows=temp[2]  ; rows = number of bands for each em
  ;  Convert A to a fltarr
  q=float(A)
  ;Do the modified Graham-Schmidt method
  n=cols ;  n= the number of endmembers
  m=rows;         m=the number of bands
  r=fltarr(n,n)
  for k = 0, n-1 do begin
    for i=0, m-1 do begin
      r(k,k)=r(k,k)+(q(k,i)*q(k,i))
    endfor
    r(k,k)=sqrt(r(k,k))
    for i = 0, m-1 do q(k,i) = q(k,i)/R(k,k)
    for J = k+1, n-1 do begin
      sum=0.
      i=0L
      repeat begin
        temp=q(k,i)*q(j,i)
        sum=sum+temp
        i=i+1
      endrep until i eq m
      r(j,K)=sum
      for i=0, m-1 do q(j,i)=q(j,i)-(q(k,i)*r(j,k))
    endfor
  endfor
return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; readbil_16 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO READBIL_16, fname, cells, lines, bands,outarr, swap=swap
  nargs = N_PARAMS(0)
  IF nargs LT 1 THEN BEGIN
    PRINT,'Usage: READBIL_16, fname, cells, lines, bands, outarr, [/swap]
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
    line=assoc(lun, intarr(cells))
    repeat begin  ; repeat over lines
      repeat begin  ; repeat over bands for line j
        ;      temp_line= (line[k] * (swap eq 0)) + (swap_endian(line[k])*(swap ne 0))
        temp_line=line[k]
        outarr(*, j, i) = temp_line
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
    print, 'Usage:READBIL_flt, fname, cells, lines, bands, outarr
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
    print, 'Usage:WRITEBIL_flt, img, fname
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;mesma_create_3em_model_file_v6;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mesma_create_3em_model_file_v6, model_file_name, em1_locs, em2_locs, em3_locs, random_subset=random_subset
;For creating all possible (or a random subset of such) models from a set of indices for emdmembers.  
;emX_locs are indices (first is 0, which should be the wvl column) placed in 1-d arrays.
  n_em1=size(em1_locs, /n_elements)
  n_em2=size(em2_locs, /n_elements)
  n_em3=size(em3_locs, /n_elements)
  print, 'Generating all possible models'
  out=[]
  for i=0, n_em1-1 do begin
    for j=0, n_em2-1 do begin
      for k=0, n_em3-1 do begin
        out=[[out],[ fix(em1_locs(i)), fix(em2_locs(j)), fix(em3_locs(k))]]
      endfor
    endfor
  endfor
  sz=size(out)
  nmodels=sz(2)
  if keyword_set(random_subset) then begin
    if nmodels le random_subset then begin
      print, 'Using all models because there are fewer possible combinations than the random subset requested...'
      sout=out
    endif else begin
      print, 'Generating random subset of models'
      rand=randomu(seed, nmodels)
      mlist=findgen(nmodels)
      mlist=mlist[sort(rand)]
      scrambled=[out[0,mlist], out[1,mlist], out[2,mlist]]
      sout=scrambled(*,0:random_subset-1)
    endelse
  endif else begin
    sout=out
    random_subset=1e19
  endelse
  get_lun, lun
  openw, lun, model_file_name
  for i=0.,min([nmodels,random_subset])-1 do printf, lun, sout(*,i), format='(3(i," "))'
  close, lun
  free_lun, lun
  return
end

;Below is code I got from Juan Guerschmann, the kind of unmixing that he does: bvls

;######################################################################
;
; Copyright (C) 1999-2010, Michele Cappellari
; E-mail: cappellari_at_astro.ox.ac.uk
;
; Updated versions of the software are available from my web page
; http://purl.org/cappellari/idl
;
; If you have found this software useful for your research,
; I would appreciate an acknowledgment and a link to the website.
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
;######################################################################
;+
; NAME:
;   BVLS
;
; AUTHOR:
;   Michele Cappellari, Leiden Observatory, The Netherlands
;   Currently at the University of Oxford UK (cappellari_at_astro.ox.ac.uk)
;
; PURPOSE:
;   Perform bound constrained linear least-squares minimization
;
; CATEGORY:
;   Least Squares
;
; CALLING SEQUENCE:
;   BVLS, A, B, BND, X, $
;       EPS=eps, /FASTNORM, IERR=ierr, INDEX=index, ITER=iter, $
;       ITMAX=itmax, NSETP=nsetp, RNORM=rnorm, W=w
;
; DESCRIPTION:
;
;   Given an M by N matrix, A(M,N), and an M-vector, B(M),  compute an
;   N-vector, X(N), that solves the least-squares problem A # X = B
;   subject to X(J) satisfying  BND(0,J) <= X(J) <= BND(1,J)
;
;   The values BND(0,J) = -(MACHAR()).XMAX and BND(1,J) = (MACHAR()).XMAX
;   are suggested choices to designate that there is no constraint in that
;   direction.
;
;   This algorithm is a generalization of  NNLS, that solves
;   the least-squares problem,  A # X = B,  subject to all X(J) >= 0.
;   The subroutine NNLS appeared in 'SOLVING LEAST SQUARES PROBLEMS,'
;   by Lawson and Hanson, Prentice-Hall, 1974.  Work on BVLS was started
;   by C. L. Lawson and R. J. Hanson at Jet Propulsion Laboratory,
;   1973 June 12.  Many modifications were subsequently made.
;   The Fortran 90 code was completed in April, 1995 by R. J. Hanson.
;   The BVLS package is an additional item for the reprinting of the book
;   by SIAM Publications and the distribution of the code package
;   using netlib and Internet or network facilities.
;
;   This IDL version was ported from the original Fortran 90 code
;   by Michele Cappellari, Leiden Observatory, The Netherlands
;
; INPUT PARAMETERS:
;
;   A(M,N)     [INTENT(InOut)]
;       On entry A() contains the M by N matrix, A.
;       On return A() contains the product matrix, Q*A, where
;       Q is an M by M orthogonal matrix generated by this
;       subroutine.  The dimensions are M=size(A,1) and N=size(A,2).
;
;   B(M)     [INTENT(InOut)]
;       On entry B() contains the M-vector, B.
;       On return, B() contains Q*B.  The same Q multiplies A.
;
;   BND(2,N)  [INTENT(In)]
;       BND(0,J) is the lower bound for X(J).
;       BND(1,J) is the upper bound for X(J).
;       Require:  BND(0,J) <= BND(1,J).
;
; OUTPUT PARAMETER:
;
;   X(N)    [INTENT(Out)]
;       On entry X() need not be initialized.  On return,
;       X() will contain the solution N-vector.
;
; KEYWORD PARAMETERS:
;
;   RNORM    [INTENT(Out)]
;       The Euclidean norm of the residual vector, b - A*X.
;
;   NSETP    [INTENT(Out)]
;       Indicates the number of components of the solution
;       vector, X(), that are not at their constraint values.
;
;   W(N)     [INTENT(Out)]
;       An N-array.  On return, W() will contain the dual solution
;       vector.   Using Set definitions below:
;       W(J) = 0 for all j in Set P,
;       W(J) <= 0 for all j in Set Z, such that X(J) is at its
;       lower bound, and
;       W(J) >= 0 for all j in Set Z, such that X(J) is at its
;       upper bound.
;       If BND(1,J) = BND(2,J), so the variable X(J) is fixed,
;       then W(J) will have an arbitrary value.
;
;   INDEX(N)    [INTENT(Out)]
;       An INTEGER working array of size N.  On exit the contents
;       of this array define the sets P, Z, and F as follows:
;
;   INDEX(1)   through INDEX(NSETP) = Set P.
;   INDEX(IZ1) through INDEX(IZ2)   = Set Z.
;   INDEX(IZ2+1) through INDEX(N)   = Set F.
;   IZ1 = NSETP + 1 = NPP1
;       Any of these sets may be empty.  Set F is those components
;       that are constrained to a unique value by the given
;       constraints.   Sets P and Z are those that are allowed a non-
;       zero range of values.  Of these, set Z are those whose final
;       value is a constraint value, while set P are those whose
;       final value is not a constraint.  The value of IZ2 is not returned.
;...    It is computable as the number of bounds constraining a component
;...    of X uniquely.
;
;   IERR    [INTENT(Out)]
;   Indicates status on return.
;   = 0   Solution completed.
;       = 1   M <= 0 or N <= 0
;       = 2   B(:), X(:), BND(:,:), W(:), or INDEX(:) size or shape violation.
;       = 3   Input bounds are inconsistent.
;       = 4   Exceed maximum number of iterations.
;
;   EPS [real(kind(one))]
;       Determines the relative linear dependence of a column vector
;       for a variable moved from its initial value.  This is used in
;       one place with the default value EPS=(MACHAR()).EPS.  Other
;       values, larger or smaller may be needed for some problems.
;
;   ITMAX  [integer]
;       Set to 3*N.  Maximum number of iterations permitted.
;       This is usually larger than required.
;
;   ITER   [integer]
;       Iteration counter.
;
;   /FASTNORM
;       Perform Euclidean Norm computation without checking for over/underflows.
;       It can speed up the program considerably when M is large, but has to
;       be used with care since may lead to instabilities!
;
; MODIFICATION HISTORY:
;   V1.0: Written by Michele Cappellari, Padova, 2000
;   V1.1: Added /FASTNORM keyword, MC, Leiden, 20 September 2001
;   V1.2: Use MAKE_ARRAY to deal with float or double arrays,
;       MC, Leiden, 19 October 2001
;   V1.3: Added compilation options and converted to IDL V5.0,
;       MC, Leiden 20 May 2002
;   V1.4: Define optional parameters using optional keywords.
;       The new calling sequence is not directly compatible with
;       the previous versions. MC, Leiden, 20 March 2004
;   V1.41: Minor updates to the documentation. MC, Oxford, 01 March 2010
;-
;----------------------------------------------------------------------
FUNCTION NRM2, X, fastNorm
  ;
  ;   NRM2 returns the Euclidean norm of a vector so that
  ;
  ;   NRM2 := sqrt( x'*x )
  ;
  ;;compile_opt IDL2, HIDDEN

  IF fastNorm THEN RETURN, SQRT(TOTAL(X^2)) ; brute force approach: use with care!

  ZERO = 0.0
  ONE = 1.0

  N = N_ELEMENTS(X)
  IF( N LT 1)THEN $
    NORM  = ZERO $
  ELSE IF( N EQ 1 )THEN $
    NORM  = ABS( X[0] ) $
  ELSE BEGIN
    SCALE = ZERO
    SSQ   = ONE
    ;
    FOR IX = 0L, N-1 DO BEGIN
      ABSXI = ABS( X[IX] )
      IF(ABSXI GT ZERO )THEN $
        IF( SCALE LT ABSXI )THEN BEGIN
        SSQ   = ONE + SSQ*( SCALE/ABSXI )^2
        SCALE = ABSXI
      ENDIF ELSE $
        SSQ   = SSQ + ( ABSXI/SCALE )^2
    ENDFOR
    NORM  = SCALE * SQRT( SSQ )
  ENDELSE
  ;
  RETURN, NORM
END
;----------------------------------------------------------------------
PRO TERMINATION
  ;compile_opt IDL2, HIDDEN

  COMMON bvls_global, FIND, HITBND, FREE1, FREE2, FREE, IERR, M, N, I, $
    IBOUND, II, IP, ITER, ITMAX, IZ, IZ1, IZ2, J, JJ, JZ, L, LBOUND, $
    NPP1, NSETP, INDEX, ZERO, ONE, TWO, A, B, S, X, W, Z, BND, ALPHA, $
    ASAVE, CC, EPS, RANGE, RNORM, NORM, SM, SS, T, UNORM, UP, ZTEST, FASTNORM
  ;
  ;    IF (IERR LE 0) THEN BEGIN
  ;
  ;   Compute the norm of the residual vector.
  SM = ZERO
  IF (NPP1 LE M) THEN $
    SM = NRM2(B[NPP1-1:M-1],fastNorm) $
  ELSE $
    W[0:N-1] = ZERO
  RNORM = SM
  ;    ENDIF; ( IERR...)
return
END ; ( TERMINATION )
;----------------------------------------------------------------------
PRO MOVE_COEF_I_FROM_SET_P_TO_SET_Z
  ;compile_opt IDL2, HIDDEN

  COMMON bvls_global

  X[I-1] = BND[IBOUND-1,I-1]
  IF (ABS(X[I-1]) GT ZERO AND JJ GT 0) THEN $
    B[0:JJ-1] = B[0:JJ-1]-A[0:JJ-1,I-1]*X[I-1]

  ;   The following loop can be null.
  FOR J=JJ+1,NSETP DO BEGIN
    II = INDEX[J-1]
    INDEX[J-2] = II

    SCALE = TOTAL(ABS(A[J-2:J-1,II-1]))
    IF (SCALE GT ZERO) THEN BEGIN
      R = SCALE * SQRT(TOTAL( (A[J-2:J-1,II-1]/SCALE)^2 ))
      IF (ABS(A[J-2,II-1]) GT ABS(A[J-1,II-1])) THEN $
        ROE = A[J-2,II-1] $
      ELSE $
        ROE = A[J-1,II-1]
      IF (ROE LT ZERO) THEN R = -R
      CC = A[J-2,II-1]/R
      SS = A[J-1,II-1]/R
      A[J-2,II-1] = R
    ENDIF ELSE BEGIN
      CC = ONE
      SS = ZERO
    ENDELSE

    SM = A[J-2,II-1]
    ;
    ;   The plane rotation is applied to two rows of A and the right-hand
    ;   side.  One row is moved to the scratch array S and THEN the updates
    ;   are computed.  The intent is for array operations to be performed
    ;   and minimal extra data movement.  One extra rotation is applied
    ;   to column II in this approach.
    S = A[J-2,0:N-1]
    A[J-2,0:N-1] = CC*S+SS*A[J-1,0:N-1]
    A[J-1,0:N-1] = CC*A[J-1,0:N-1]-SS*S
    A[J-2,II-1] = SM
    A[J-1,II-1] = ZERO
    SM = B[J-2]
    B[J-2] = CC*SM+SS*B[J-1]
    B[J-1] = CC*B[J-1]-SS*SM
  ENDFOR
  ;
  NPP1 = NSETP
  NSETP = NSETP-1
  IZ1 = IZ1-1
  INDEX[IZ1-1] = I
return
END ; ( MOVE COEF I FROM SET P TO SET Z )
;----------------------------------------------------------------------
PRO SEE_IF_ALL_CONSTRAINED_COEFFS_ARE_FEASIBLE
  ;compile_opt IDL2, HIDDEN

  COMMON bvls_global
  ;
  ;   See if each coefficient in set P is strictly interior to its constraint region.
  ;   If so, set HITBND = false.
  ;   If not, set HITBND = true, and also set ALPHA, JJ, and IBOUND.
  ;   Then ALPHA will satisfy  0.  < ALPHA  <=  1.
  ;
  ALPHA=TWO
  FOR IP=1L,NSETP DO BEGIN
    L = INDEX[IP-1]
    IF  (Z[IP-1]  LE  BND[0,L-1]) THEN $
      ;   Z(IP) HITS LOWER BOUND
      LBOUND=1 $
    ELSE  IF  (Z[IP-1]  GE  BND[1,L-1]) THEN $
      ;   Z(IP) HITS UPPER BOUND
      LBOUND=2 $
    ELSE $
      LBOUND = 0
    ;
    IF  ( LBOUND   NE   0 ) THEN BEGIN
      T = (BND[LBOUND-1,L-1]-X[L-1])/(Z[IP-1]-X[L-1])
      IF  (ALPHA   GT  T) THEN BEGIN
        ALPHA = T
        JJ = IP
        IBOUND = LBOUND
      ENDIF; ( LBOUND )
    ENDIF; ( ALPHA   >  T )
  ENDFOR
  HITBND = ABS(ALPHA - TWO) GT ZERO
return
END ;( SEE IF ALL CONSTRAINED COEFFS ARE FEASIBLE )
;----------------------------------------------------------------------
PRO TEST_SET_P_AGAINST_CONSTRAINTS
  ;compile_opt IDL2, HIDDEN

  COMMON bvls_global

  WHILE 1 DO BEGIN
    ;   The solution obtained by solving the current set P is in the array Z().
    ;
    ITER = ITER+1
    IF (ITER GT ITMAX) THEN BEGIN
      IERR = 4
      GOTO, fine_LOOPB
    ENDIF
    ;
    SEE_IF_ALL_CONSTRAINED_COEFFS_ARE_FEASIBLE
    ;
    ;   The above call sets HITBND.  If HITBND = true THEN it also sets
    ;   ALPHA, JJ, and IBOUND.
    IF (NOT HITBND) THEN GOTO, fine_LOOPB
    ;
    ;   Here ALPHA will be between 0 and 1 for interpolation
    ;   between the old X() and the new Z().
    FOR IP=1L,NSETP DO BEGIN
      L = INDEX[IP-1]
      X[L-1] = X[L-1]+ALPHA*(Z[IP-1]-X[L-1])
    ENDFOR
    ;
    I = INDEX[JJ-1]
    ;   Note:  The exit test is done at the end of the loop, so the loop
    ;   will always be executed at least once.
    WHILE 1 DO BEGIN
      ;
      ;   Modify A(*,*), B(*) and the index arrays to move coefficient I
      ;   from set P to set Z.
      ;
      MOVE_COEF_I_FROM_SET_P_TO_SET_Z
      ;
      IF (NSETP LE 0) THEN GOTO, fine_LOOPB
      ;
      ;   See if the remaining coefficients in set P are feasible.  They should
      ;   be because of the way ALPHA was determined.  If any are infeasible
      ;   it is due to round-off error.  Any that are infeasible or on a boundary
      ;   will be set to the boundary value and moved from set P to set Z.
      ;
      IBOUND = 0
      FOR JJ=1L,NSETP DO BEGIN
        I = INDEX[JJ-1]
        IF  (X[I-1] LE BND[0,I-1]) THEN BEGIN
          IBOUND = 1
          GOTO, fine_ciclo1
        ENDIF ELSE IF (X[I-1] GE BND[1,I-1]) THEN BEGIN
          IBOUND = 2
          GOTO, fine_ciclo1
        ENDIF
      ENDFOR
      fine_ciclo1:
      IF (IBOUND LE 0) THEN GOTO, fine_ciclo2
    ENDWHILE
    fine_ciclo2:
    ;
    ;   Copy B( ) into Z( ).  Then solve again and loop back.
    Z[0:M-1] = B[0:M-1]
    ;
    FOR I=NSETP,1,-1 DO BEGIN
      IF (I NE NSETP) THEN $
        Z[0:I-1] = Z[0:I-1]-A[0:I-1,II-1]*Z[I]
      II = INDEX[I-1]
      Z[I-1] = Z[I-1]/A[I-1,II-1]
    ENDFOR
  ENDWHILE
  fine_LOOPB:

  ;   The following loop can be null.
  FOR IP=1L,NSETP DO BEGIN
    I = INDEX[IP-1]
    X[I-1] = Z[IP-1]
  ENDFOR
return
END ; ( TEST SET P AGAINST CONSTRAINTS)
;----------------------------------------------------------------------
PRO MOVE_J_FROM_SET_Z_TO_SET_P
  ;compile_opt IDL2, HIDDEN

  COMMON bvls_global
  ;
  ;   The index  J=index(IZ)  has been selected to be moved from
  ;   set Z to set P.  Z() contains the old B() adjusted as though X(J) = 0.
  ;   A(*,J) contains the new Householder transformation vector.
  B[0:M-1] = Z[0:M-1]
  ;
  INDEX[IZ-1] = INDEX[IZ1-1]
  INDEX[IZ1-1] = J
  IZ1 = IZ1+1
  NSETP = NPP1
  NPP1 = NPP1+1
  ;   The following loop can be null or not required.
  NORM = A[NSETP-1,J-1]
  A[NSETP-1,J-1] = UP
  IF(ABS(NORM) GT ZERO) THEN BEGIN
    FOR JZ=IZ1,IZ2 DO BEGIN
      JJ = INDEX[JZ-1]
      SM = TOTAL(A[NSETP-1:M-1,J-1]/NORM * A[NSETP-1:M-1,JJ-1])/UP
      A[NSETP-1:M-1,JJ-1] = A[NSETP-1:M-1,JJ-1]+SM*A[NSETP-1:M-1,J-1]
    ENDFOR
    A[NSETP-1,J-1] = NORM
  ENDIF
  ;   The following loop can be null.
  FOR L=NPP1,M DO A[L-1,J-1] = ZERO
  ;
  W[J-1] = ZERO
  ;
  ;   Solve the triangular system.  Store this solution temporarily in Z().
  FOR I=NSETP,1,-1 DO BEGIN
    IF (I NE NSETP) THEN Z[0:I-1] = Z[0:I-1]-A[0:I-1,II-1]*Z[I]
    II = INDEX[I-1]
    Z[I-1] = Z[I-1]/A[I-1,II-1]
  ENDFOR
return
END ; ( MOVE J FROM SET Z TO SET P )
;----------------------------------------------------------------------
PRO TEST_COEF_J_FOR_DIAG_ELT_AND_DIRECTION_OF_CHANGE
  ;compile_opt IDL2, HIDDEN

  COMMON bvls_global
  ;
  ;   The sign of W(J) is OK for J to be moved to set P.
  ;   Begin the transformation and check new diagonal element to avoid
  ;   near linear dependence.
  ;
  ASAVE = A[NPP1-1,J-1]
  ;
  ;   Construct a Householder transformation.

  VNORM = NRM2(A[NPP1-1:M-1,J-1],fastNorm)
  IF (A[NPP1-1,J-1] GT ZERO) THEN VNORM = -VNORM
  UP = A[NPP1-1,J-1] - VNORM
  A[NPP1-1,J-1] = VNORM

  IF (NSETP LT 1) THEN UNORM = 0.0 ELSE UNORM = NRM2(A[0:NSETP-1,J-1],fastNorm)
  IF (ABS(A[NPP1-1,J-1]) GT EPS*UNORM) THEN BEGIN
    ;
    ;   Column J is sufficiently independent.  Copy b into Z, update Z.
    Z[0:M-1] = B[0:M-1]
    ; Compute product of transormation and updated right-hand side.
    NORM = A[NPP1-1,J-1]
    A[NPP1-1,J-1] = UP
    IF (ABS(NORM) GT ZERO) THEN BEGIN
      SM = TOTAL(A[NPP1-1:M-1,J-1]/NORM * Z[NPP1-1:M-1])/UP
      Z[NPP1-1:M-1] = Z[NPP1-1:M-1]+SM*A[NPP1-1:M-1,J-1]
      A[NPP1-1,J-1] = NORM
    ENDIF

    IF (ABS(X[J-1]) GT ZERO) THEN $
      Z[0:NPP1-1] = Z[0:NPP1-1]+A[0:NPP1-1,J-1]*X[J-1]
    ;   Adjust Z() as though X(J) had been reset to zero.
    IF ( FREE ) THEN $
      FIND = 1 $
    ELSE BEGIN
      ;
      ;   Solve for ZTEST ( proposed new value for X(J) ).
      ;   Then set FIND to indicate if ZTEST has moved away from X(J) in
      ;   the expected direction indicated by the sign of W(J).
      ZTEST = Z[NPP1-1]/A[NPP1-1,J-1]
      FIND = ( W[J-1] LT ZERO AND ZTEST LT X[J-1] ) OR $
        ( W[J-1] GT ZERO AND ZTEST GT X[J-1] )
    ENDELSE
  ENDIF
  ;
  ;   If J was not accepted to be moved from set Z to set P,
  ;   restore A(NNP1,J).  Failing these tests may mean the computed
  ;   sign of W(J) is suspect, so here we set W(J) = 0.  This will
  ;   not affect subsequent computation, but cleans up the W() array.
  IF  ( NOT FIND ) THEN BEGIN
    A[NPP1-1,J-1] = ASAVE
    W[J-1] = ZERO
  ENDIF; ( .not. FIND )
return
END ;TEST_COEF_J_FOR_DIAG_ELT_AND_DIRECTION_OF_CHANGE
;----------------------------------------------------------------------
PRO  SELECT_ANOTHER_COEFF_TO_SOLVE_FOR
  ;compile_opt IDL2, HIDDEN

  COMMON bvls_global
  ;
  ;   1. Search through set z for a new coefficient to solve for.
  ;   First select a candidate that is either an unconstrained
  ;   coefficient or ELSE a constrained coefficient that has room
  ;   to move in the direction consistent with the sign of its dual
  ;   vector component.  Components of the dual (negative gradient)
  ;   vector will be computed as needed.
  ;   2. For each candidate start the transformation to bring this
  ;   candidate into the triangle, and THEN do two tests:  Test size
  ;   of new diagonal value to avoid extreme ill-conditioning, and
  ;   the value of this new coefficient to be sure it moved in the
  ;   expected direction.
  ;   3. If some coefficient passes all these conditions, set FIND = true,
  ;   The index of the selected coefficient is J = INDEX(IZ).
  ;   4. If no coefficient is selected, set FIND = false.
  ;
  FIND = 0
  FOR IZ=IZ1,IZ2 DO BEGIN
    J = INDEX[IZ-1]
    ;
    ;   Set FREE1 = true if X(J) is not at the left end-point of its
    ;   constraint region.
    ;   Set FREE2 = true if X(J) is not at the right end-point of its
    ;   constraint region.
    ;   Set FREE = true if X(J) is not at either end-point of its
    ;   constraint region.
    ;
    FREE1 = X[J-1] GT BND[0,J-1]
    FREE2 = X[J-1] LT BND[1,J-1]
    FREE = FREE1 AND FREE2

    IF ( FREE ) THEN $
      TEST_COEF_J_FOR_DIAG_ELT_AND_DIRECTION_OF_CHANGE $
    ELSE BEGIN
      ;   Compute dual coefficient W(J).
      W[J-1] = TOTAL(A[NPP1-1:M-1,J-1] * B[NPP1-1:M-1])
      ;
      ;   Can X(J) move in the direction indicated by the sign of W(J)?
      ;
      IF ( W[J-1] LT ZERO ) THEN BEGIN
        IF ( FREE1 ) THEN $
          TEST_COEF_J_FOR_DIAG_ELT_AND_DIRECTION_OF_CHANGE
      ENDIF ELSE  IF ( W[J-1]  GT ZERO ) THEN BEGIN
        IF ( FREE2 ) THEN $
          TEST_COEF_J_FOR_DIAG_ELT_AND_DIRECTION_OF_CHANGE
      ENDIF
    ENDELSE
    IF ( FIND ) THEN RETURN
  ENDFOR;  IZ
return
END ; ( SELECT ANOTHER COEF TO SOLVE FOR )
;----------------------------------------------------------------------
PRO INITIALIZE, eps1, itmax1
  ;compile_opt IDL2, HIDDEN

  COMMON bvls_global

  IF (N LT 2 OR M LT 2 OR N_ELEMENTS(B) NE M $
    OR (SIZE(BND))[1] NE 2 OR (SIZE(BND))[2] NE N) THEN BEGIN
    IERR = 2
    MESSAGE, 'Wrong input arrays size'
  ENDIF

  IERR = 0
  mch = MACHAR()
  IF N_ELEMENTS(eps1) EQ 0 THEN EPS = mch.EPS ELSE EPS = eps1
  IF N_ELEMENTS(itmax1) EQ 0 THEN ITMAX = 3L*N ELSE ITMAX = itmax1
  ITER = 0L
  ;
  IZ2 = N
  IZ1 = 1L
  NSETP = 0L
  NPP1 = 1L
  ;
  ;   Begin:  Loop on IZ to initialize  X().
  IZ = IZ1
  WHILE 1 DO BEGIN
    IF (IZ GT IZ2) THEN GOTO, fine1
    J = INDEX[IZ-1]
    IF (BND[0,J-1] LE -mch.XMAX) THEN $
      IF (BND[1,J-1] GE mch.XMAX) THEN $
      X[J-1] = ZERO $
    ELSE $
      X[J-1] = ZERO < BND[1,J-1] $
    ELSE IF (BND[1,J-1] GE mch.XMAX) THEN $
      X[J-1] = ZERO > BND[0,J-1] $
    ELSE BEGIN
      RANGE = BND[1,J-1] - BND[0,J-1]
      IF (RANGE LE ZERO) THEN BEGIN
        ;
        ;   Here X(J) is constrained to a single value.
        INDEX[IZ-1] = INDEX[IZ2-1]
        INDEX[IZ2-1] = J
        IZ = IZ-1
        IZ2 = IZ2-1
        X[J-1] = BND[0,J-1]
        W[J-1] = ZERO
      ENDIF ELSE IF (RANGE GT ZERO) THEN $
        ;
        ;   The following statement sets X(J) to 0 if the constraint interval
        ;   includes 0, and otherwise sets X(J) to the endpoint of the
        ;   constraint interval that is closest to 0.
        ;
        X[J-1] = BND[0,J-1] > (BND[1,J-1] < ZERO) $
      ELSE BEGIN
        IERR = 3
        RETURN
      ENDELSE ; ( RANGE:.)
    ENDELSE
    ;
    ;   Change B() to reflect a nonzero starting value for X(J).
    ;
    IF (ABS(X[J-1]) GT ZERO) THEN $
      B[0:M-1] = B[0:M-1]-A[0:M-1,J-1]*X[J-1]
    IZ = IZ+1
  ENDWHILE; ( IZ <= IZ2 )
  fine1:
return
END ; ( INITIALIZE )
;----------------------------------------------------------------------
PRO BVLS, A1, B1, BND1, X1, $
  EPS=eps1, FASTNORM=fastNorm1, IERR=ierr1, INDEX=index1, $
  ITER=iter1, ITMAX=itmax1, NSETP=nsetp1, RNORM=rnorm1, W=w1

  ;compile_opt IDL2
  ;ON_ERROR, 2

  COMMON bvls_global
  ;
  ; Load needed input parameters into the COMMON block
  ;
  A = TEMPORARY(A1)
  B = TEMPORARY(B1)
  BND = TEMPORARY(BND1)

  siz = SIZE(A)
  M = siz[1]
  N = siz[2]
  X = MAKE_ARRAY(N,TYPE=siz[3])
  W = MAKE_ARRAY(N,TYPE=siz[3])
  INDEX = LINDGEN(N)+1

  S = MAKE_ARRAY(N,TYPE=siz[3])
  Z = MAKE_ARRAY(M,TYPE=siz[3])

  ; Load some constants into the COMMON block
  ;
  ZERO = 0.0
  ONE = 1.0
  TWO = 2.0
  IF KEYWORD_SET(fastNorm1) THEN fastNorm = 1 ELSE fastNorm = 0

  INITIALIZE, eps1, itmax1
  ;
  ;   The above call will set IERR.
  ;
  WHILE 1 DO BEGIN
    ;
    ;   Quit on error flag, or if all coefficients are already in the
    ;   solution, .or. if M columns of A have been triangularized.
    IF (IERR NE 0 OR IZ1 GT IZ2 OR NSETP GE M) THEN GOTO, fine
    ;
    SELECT_ANOTHER_COEFF_TO_SOLVE_FOR
    ;
    ;   See if no index was found to be moved from set Z to set P.
    ;   Then go to termination.
    IF  ( NOT FIND ) THEN GOTO, fine
    ;
    MOVE_J_FROM_SET_Z_TO_SET_P
    ;
    TEST_SET_P_AGAINST_CONSTRAINTS
    ;
    ;   The above call may set IERR.
    ;   All coefficients in set P are strictly feasible.  Loop back.
  ENDWHILE
  fine:
  ;
  TERMINATION

  A1 = TEMPORARY(A)
  B1 = TEMPORARY(B)
  BND1 = TEMPORARY(BND)
  X1 = TEMPORARY(X)
  RNORM1 = RNORM
  NSETP1 = NSETP
  W1 = TEMPORARY(W)
  INDEX1 = TEMPORARY(INDEX)
  IERR1 = IERR
  ITER1 = ITER
return
END
;----------------------------------------------------------------------


