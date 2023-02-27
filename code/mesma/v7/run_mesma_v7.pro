pro run_mesma_v7

n=120
ns=10  ;also need to change in mesma_v7,initialize variables
nl=10    ;also need to change in mesma_v7,initialize variables
nb=7    ;also need to change in mesma_v7,initialize variables
in_file='D:\MCD43A4.A2008153.h08v05.005.2008176053238.hdf'


for i=0,1 do begin
;  get_lun, wlun
;  openw, wlun, input_file
;  printf, lun, ';MESMA SMA Input File (v6)
;  printf, lun, '; name of file containing Endmember Spectra (ASCII) (tab delimited, with header, and band first colum; if hyperspectral, first column must be wavelength in nm)
;  printf, lun, 'combined_library_isric_aster_thoralf_814x2151.txt
;  printf, lun, '; Scale value for endmember spectra(what value does 100% reflectance have?)
;  printf, lun, '1.
;  printf, lun, ';Set this to 0 or the name of an instrument band file if a hyperspectral library is to be convolved to instrument wavebands
;  printf, lun, 'MODIS_Instrument_Filter_File_For_MODTRAN.flt
;  printf, lun, '; Model File - Spectra #s refer to number among the em spectra, not columns in file (since the library file must have a first column with band info, the second column is the first em
;  printf, lun, 'models_subset_1000_combined_library_isric_aster_thoralf_814x2151.txt
;  printf, lun, ';Set this to 0 or the name of a MODIS NBAR .hdf  if that is the input
;  printf, lun, '0
;  printf, lun, '; name of file containing the reflectance image (16-bit BIL) ;Ignored if modis hdf
;  printf, lun, 'Australia_mainland.bil
;  printf, lun, '; Swap (1 = yes, any other number = no) ;Ignored if modis hdf
;  printf, lun, '0
;  printf, lun, '; Image null value
;  printf, lun, '32767
;  printf, lun, '; Subset image? (1 = yes, any other number = no)
;  printf, lun, '0
;  printf, lun, ';If subsetting image, the inclusive starting X location (count from 0) ; ignored if subset image ne 1
;  printf, lun, '200
;  printf, lun, ';If subsetting image, the inclusive ending X location (count from 0) ; ignored if subset image ne 1
;  printf, lun, '209
;  printf, lun, ';If subsetting image, the inclusive starting Y location (count from 0) ; ignored if subset image ne 1
;  printf, lun, '200
;  printf, lun, ';If subsetting image, the inclusive ending Y location (count from 0) ; ignored if subset image ne 1
;  printf, lun, '209
;  printf, lun, '; X dimension of the BIL (samples) ;Ignored if modis hdf
;  printf, lun, '240
;  printf, lun, '; Y dimension of the BIL (rows) ;Ignored if modis hdf
;  printf, lun, '2
;  printf, lun, '; number of bands ;Ignored if modis hdf
;  printf, lun, '7
;  printf, lun, '; Scale value for reflectance image (what value does 100% reflectance have?) ;Ignored if modis hdf
;  printf, lun, '10000
;  printf, lun, '; Output file name (bil) - will write same byte order as reflectance file
;  printf, lun, strcompress('mesma_out_'+string(i)+'.bil', /rem)
;  printf, lun, ';Maximum RMS Error (reflectance units) (needs to be higher for guerschmann transform)
;  printf, lun, '100
;  printf, lun, ';Maximum allowable fraction value (usually ~ 1.0)
;  printf, lun, '1.0
;  printf, lun, ';Minimum allowable fraction value (usually ~ 0.0)
;  printf, lun, '-0.0
;  printf, lun, ';Use an image as the source for one of the endmembers, model list should still have n endmembers (not n-1)
;  printf, lun, '0
;  printf, lun, ';Which endmember will be replaced with the value from a pixel from the image (1 = first endmember) Duplicate mixtures of em2 and em3 will be run (so watch out!)
;  printf, lun, '3
;  printf, lun, ';What is the path for the file that will be used as the endmember image (the source of the pixels)? Size must be same as subsetted size, reflectance range = [0,1]
;  printf, lun, 'EM_3_spectra_flt_MCD43A4.A2010129.h19v07_mesma_out_515_1230_6_9models_standardize=4.flt.bil'  
;  
;  close, wlun
;  free_lun, wlun
  print,'$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'
  print,i
  case i of
    0:  mesma_v7,i, in_file, 1
    1:  mesma_v7,i, in_file, 1, sum_to_one_weight=1.0
    2:  mesma_v7,i, in_file, 1, sum_to_one_weight=1.0, standardize_band=5
    3:  mesma_v7,i, in_file, 1, sum_to_one_weight=1.0, /equal_brightness
    4:  mesma_v7,i, in_file, 1, /guerschman
    5:  mesma_v7,i, in_file, 1, /guerschman, sum_to_one_weight=0.02
    6:  mesma_v7,i, in_file, 1, /guerschman, /use_explicit
    7:  mesma_v7,i, in_file, 1, /guerschman, /use_explicit, sum_to_one_weight=0.02
    8:  mesma_v7,i, in_file, 1, /normalize_fractions_before_rms
    9:  mesma_v7,i, in_file, 1, sum_to_one_weight=1.0, /normalize_fractions_before_rms
    10:  mesma_v7,i, in_file, 1, sum_to_one_weight=1.0, standardize_band=5, /normalize_fractions_before_rms
    11:  mesma_v7,i, in_file, 1, sum_to_one_weight=1.0, /equal_brightness, /normalize_fractions_before_rms
    12:  mesma_v7,i, in_file, 1, /guerschman, /normalize_fractions_before_rms
    13:  mesma_v7,i, in_file, 1, /guerschman, sum_to_one_weight=0.02, /normalize_fractions_before_rms
    14:  mesma_v7,i, in_file, 1, /guerschman, /use_explicit, /normalize_fractions_before_rms
    15:  mesma_v7,i, in_file, 1, /guerschman, /use_explicit, sum_to_one_weight=0.02, /normalize_fractions_before_rms
    16:  mesma_v7,i, in_file, 1, /normalize_fractions_after_rms
    17:  mesma_v7,i, in_file, 1, sum_to_one_weight=1.0, /normalize_fractions_after_rms
    18:  mesma_v7,i, in_file, 1, sum_to_one_weight=1.0, standardize_band=5, /normalize_fractions_after_rms
    19:  mesma_v7,i, in_file, 1, sum_to_one_weight=1.0, /equal_brightness, /normalize_fractions_after_rms
    20:  mesma_v7,i, in_file, 1, /guerschman, /normalize_fractions_after_rms
    21:  mesma_v7,i, in_file, 1, /guerschman, sum_to_one_weight=0.02, /normalize_fractions_after_rms
    22:  mesma_v7,i, in_file, 1, /guerschman, /use_explicit, /normalize_fractions_after_rms
    23:  mesma_v7,i, in_file, 1, /guerschman, /use_explicit, sum_to_one_weight=0.02, /normalize_fractions_after_rms
    24:  mesma_v7,i, in_file, 2
    25:  mesma_v7,i, in_file, 2, sum_to_one_weight=1.0
    26:  mesma_v7,i, in_file, 2, sum_to_one_weight=1.0, standardize_band=5
    27:  mesma_v7,i, in_file, 2, sum_to_one_weight=1.0, /equal_brightness
    28:  mesma_v7,i, in_file, 2, /guerschman
    29:  mesma_v7,i, in_file, 2, /guerschman, sum_to_one_weight=0.02
    30:  mesma_v7,i, in_file, 2, /guerschman, /use_explicit
    31:  mesma_v7,i, in_file, 2, /guerschman, /use_explicit, sum_to_one_weight=0.02
    32:  mesma_v7,i, in_file, 2, /normalize_fractions_before_rms
    33:  mesma_v7,i, in_file, 2, sum_to_one_weight=1.0, /normalize_fractions_before_rms
    34:  mesma_v7,i, in_file, 2, sum_to_one_weight=1.0, standardize_band=5, /normalize_fractions_before_rms
    35:  mesma_v7,i, in_file, 2, sum_to_one_weight=1.0, /equal_brightness, /normalize_fractions_before_rms
    36:  mesma_v7,i, in_file, 2, /guerschman, /normalize_fractions_before_rms
    37:  mesma_v7,i, in_file, 2, /guerschman, sum_to_one_weight=0.02, /normalize_fractions_before_rms
    38:  mesma_v7,i, in_file, 2, /guerschman, /use_explicit, /normalize_fractions_before_rms
    39:  mesma_v7,i, in_file, 2, /guerschman, /use_explicit, sum_to_one_weight=0.02, /normalize_fractions_before_rms
    40:  mesma_v7,i, in_file, 2, /normalize_fractions_after_rms
    41:  mesma_v7,i, in_file, 2, sum_to_one_weight=1.0, /normalize_fractions_after_rms
    42:  mesma_v7,i, in_file, 2, sum_to_one_weight=1.0, standardize_band=5, /normalize_fractions_after_rms
    43:  mesma_v7,i, in_file, 2, sum_to_one_weight=1.0, /equal_brightness, /normalize_fractions_after_rms
    44:  mesma_v7,i, in_file, 2, /guerschman, /normalize_fractions_after_rms
    45:  mesma_v7,i, in_file, 2, /guerschman, sum_to_one_weight=0.02, /normalize_fractions_after_rms
    46:  mesma_v7,i, in_file, 2, /guerschman, /use_explicit, /normalize_fractions_after_rms
    47:  mesma_v7,i, in_file, 2, /guerschman, /use_explicit, sum_to_one_weight=0.02, /normalize_fractions_after_rms
    48:  mesma_v7,i, in_file, 3
    49:  mesma_v7,i, in_file, 3, sum_to_one_weight=1.0
    50:  mesma_v7,i, in_file, 3, sum_to_one_weight=1.0, standardize_band=5
    51:  mesma_v7,i, in_file, 3, sum_to_one_weight=1.0, /equal_brightness
    52:  mesma_v7,i, in_file, 3, /guerschman
    53:  mesma_v7,i, in_file, 3, /guerschman, sum_to_one_weight=0.02
    54:  mesma_v7,i, in_file, 3, /guerschman, /use_explicit
    55:  mesma_v7,i, in_file, 3, /guerschman, /use_explicit, sum_to_one_weight=0.02
    56:  mesma_v7,i, in_file, 3, /normalize_fractions_before_rms
    57:  mesma_v7,i, in_file, 3, sum_to_one_weight=1.0, /normalize_fractions_before_rms
    58:  mesma_v7,i, in_file, 3, sum_to_one_weight=1.0, standardize_band=5, /normalize_fractions_before_rms
    59:  mesma_v7,i, in_file, 3, sum_to_one_weight=1.0, /equal_brightness, /normalize_fractions_before_rms
    60:  mesma_v7,i, in_file, 3, /guerschman, /normalize_fractions_before_rms
    61:  mesma_v7,i, in_file, 3, /guerschman, sum_to_one_weight=0.02, /normalize_fractions_before_rms
    62:  mesma_v7,i, in_file, 3, /guerschman, /use_explicit, /normalize_fractions_before_rms
    63:  mesma_v7,i, in_file, 3, /guerschman, /use_explicit, sum_to_one_weight=0.02, /normalize_fractions_before_rms
    64:  mesma_v7,i, in_file, 3, /normalize_fractions_after_rms
    65:  mesma_v7,i, in_file, 3, sum_to_one_weight=1.0, /normalize_fractions_after_rms
    66:  mesma_v7,i, in_file, 3, sum_to_one_weight=1.0, standardize_band=5, /normalize_fractions_after_rms
    67:  mesma_v7,i, in_file, 3, sum_to_one_weight=1.0, /equal_brightness, /normalize_fractions_after_rms
    68:  mesma_v7,i, in_file, 3, /guerschman, /normalize_fractions_after_rms
    69:  mesma_v7,i, in_file, 3, /guerschman, sum_to_one_weight=0.02, /normalize_fractions_after_rms
    70:  mesma_v7,i, in_file, 3, /guerschman, /use_explicit, /normalize_fractions_after_rms
    71:  mesma_v7,i, in_file, 3, /guerschman, /use_explicit, sum_to_one_weight=0.02, /normalize_fractions_after_rms
    72:  mesma_v7,i, in_file, 4
    73:  mesma_v7,i, in_file, 4, sum_to_one_weight=1.0
    74:  mesma_v7,i, in_file, 4, sum_to_one_weight=1.0, standardize_band=5
    75:  mesma_v7,i, in_file, 4, sum_to_one_weight=1.0, /equal_brightness
    76:  mesma_v7,i, in_file, 4, /guerschman
    77:  mesma_v7,i, in_file, 4, /guerschman, sum_to_one_weight=0.02
    78:  mesma_v7,i, in_file, 4, /guerschman, /use_explicit
    79:  mesma_v7,i, in_file, 4, /guerschman, /use_explicit, sum_to_one_weight=0.02
    80:  mesma_v7,i, in_file, 4, /normalize_fractions_before_rms
    81:  mesma_v7,i, in_file, 4, sum_to_one_weight=1.0, /normalize_fractions_before_rms
    82:  mesma_v7,i, in_file, 4, sum_to_one_weight=1.0, standardize_band=5, /normalize_fractions_before_rms
    83:  mesma_v7,i, in_file, 4, sum_to_one_weight=1.0, /equal_brightness, /normalize_fractions_before_rms
    84:  mesma_v7,i, in_file, 4, /guerschman, /normalize_fractions_before_rms
    85:  mesma_v7,i, in_file, 4, /guerschman, sum_to_one_weight=0.02, /normalize_fractions_before_rms
    86:  mesma_v7,i, in_file, 4, /guerschman, /use_explicit, /normalize_fractions_before_rms
    87:  mesma_v7,i, in_file, 4, /guerschman, /use_explicit, sum_to_one_weight=0.02, /normalize_fractions_before_rms
    88:  mesma_v7,i, in_file, 4, /normalize_fractions_after_rms
    89:  mesma_v7,i, in_file, 4, sum_to_one_weight=1.0, /normalize_fractions_after_rms
    90:  mesma_v7,i, in_file, 4, sum_to_one_weight=1.0, standardize_band=5, /normalize_fractions_after_rms
    91:  mesma_v7,i, in_file, 4, sum_to_one_weight=1.0, /equal_brightness, /normalize_fractions_after_rms
    92:  mesma_v7,i, in_file, 4, /guerschman, /normalize_fractions_after_rms
    93:  mesma_v7,i, in_file, 4, /guerschman, sum_to_one_weight=0.02, /normalize_fractions_after_rms
    94:  mesma_v7,i, in_file, 4, /guerschman, /use_explicit, /normalize_fractions_after_rms
    95:  mesma_v7,i, in_file, 4, /guerschman, /use_explicit, sum_to_one_weight=0.02, /normalize_fractions_after_rms
    96:  mesma_v7,i, in_file, 5
    97:  mesma_v7,i, in_file, 5, sum_to_one_weight=1.0
    98:  mesma_v7,i, in_file, 5, sum_to_one_weight=1.0, standardize_band=5
    99:  mesma_v7,i, in_file, 5, sum_to_one_weight=1.0, /equal_brightness
    100:  mesma_v7,i, in_file, 5, /guerschman
    101:  mesma_v7,i, in_file, 5, /guerschman, sum_to_one_weight=0.02
    102:  mesma_v7,i, in_file, 5, /guerschman, /use_explicit
    103:  mesma_v7,i, in_file, 5, /guerschman, /use_explicit, sum_to_one_weight=0.02
    104:  mesma_v7,i, in_file, 5, /normalize_fractions_before_rms
    105:  mesma_v7,i, in_file, 5, sum_to_one_weight=1.0, /normalize_fractions_before_rms
    106:  mesma_v7,i, in_file, 5, sum_to_one_weight=1.0, standardize_band=5, /normalize_fractions_before_rms
    107:  mesma_v7,i, in_file, 5, sum_to_one_weight=1.0, /equal_brightness, /normalize_fractions_before_rms
    108:  mesma_v7,i, in_file, 5, /guerschman, /normalize_fractions_before_rms
    109:  mesma_v7,i, in_file, 5, /guerschman, sum_to_one_weight=0.02, /normalize_fractions_before_rms
    110:  mesma_v7,i, in_file, 5, /guerschman, /use_explicit, /normalize_fractions_before_rms
    111:  mesma_v7,i, in_file, 5, /guerschman, /use_explicit, sum_to_one_weight=0.02, /normalize_fractions_before_rms
    112:  mesma_v7,i, in_file, 5, /normalize_fractions_after_rms
    113:  mesma_v7,i, in_file, 5, sum_to_one_weight=1.0, /normalize_fractions_after_rms
    114:  mesma_v7,i, in_file, 5, sum_to_one_weight=1.0, standardize_band=5, /normalize_fractions_after_rms
    115:  mesma_v7,i, in_file, 5, sum_to_one_weight=1.0, /equal_brightness, /normalize_fractions_after_rms
    116:  mesma_v7,i, in_file, 5, /guerschman, /normalize_fractions_after_rms
    117:  mesma_v7,i, in_file, 5, /guerschman, sum_to_one_weight=0.02, /normalize_fractions_after_rms
    118:  mesma_v7,i, in_file, 5, /guerschman, /use_explicit, /normalize_fractions_after_rms
    119:  mesma_v7,i, in_file, 5, /guerschman, /use_explicit, sum_to_one_weight=0.02, /normalize_fractions_after_rms
  endcase
  
;  out=fltarr(ns,nl,nb)
;  readbil_flt, file_name, ns, nl, nb, img
;  for j=0,ns-1 do begin
;      out[j,*,*]=img[j,*,*]
;  endfor
;  write_csv, strcompress('D:\mesma_out_'+string(i)+'.csv', /rem),  out[*,*,0], out[*,*,1], out[*,*,2], out[*,*,3], out[*,*,4], out[*,*,5]
;  writebil_flt, out, strcompress('D:\mesma_out_'+string(i)+'.bil',/rem)
;
;  out_ave=fltarr(ns,nl,nb)
;  readbil_flt, file_name_ave, ns, nl, nb, img_ave
;  for j=0,ns-1 do begin
;      out_ave[j,*,*]=img_ave[j,*,*]
;  endfor
;  write_csv, strcompress('D:\mesma_out_averages_'+string(i)+'.csv', /rem), out_ave[*,*,0], out_ave[*,*,1], out_ave[*,*,2], out_ave[*,*,3], out_ave[*,*,4], out_ave[*,*,5]
;  writebil_flt, out_ave, strcompress('D:\mesma_out_averages_'+string(i)+'.bil',/rem)   
  print,'$$$$$$$$$$$$$$$$$$$$$$$$$'
endfor

return
end
