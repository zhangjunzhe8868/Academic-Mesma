pro run_mesma_v6

;modis v5 is 8-day mean, so there is only mean image
;modis v6 is daily, so there are +4 to -4 images and calculated mean
input_file='D:\qinghai_2016_201.dat'
;temp_str = strsplit(input_file,'.',/extract)
ns=10
nl=10
nb=7
;file_name=temp_str[0]+'_mesma1.bil'
;file_name_ave=temp_str[0]+'_mesma_ave1.bil'

for i=0,1 do begin
  print,'#############################'
  print,i
  case i of
    0: mesma_v6_tif, i, input_file, /force_lls
    1: mesma_v6_tif, i, input_file, /force_lls, lls_sum_to_one=1.
    2: mesma_v6_tif, i, input_file, /force_lls, /normalize_fractions_before_rms
    3: mesma_v6_tif, i, input_file, /force_lls, /normalize_fractions_after_rms
    4: mesma_v6_tif, i, input_file, /force_lls, /normalize_fractions_before_rms, lls_sum_to_one=1.
    5: mesma_v6_tif, i, input_file, /force_lls, /normalize_fractions_after_rms, lls_sum_to_one=1.
    6: mesma_v6_tif, i, input_file, /force_lls, lls_sum_to_one=1., standardize_band=5
    7: mesma_v6_tif, i, input_file, /force_lls, /normalize_fractions_before_rms, standardize_band=5
    8: mesma_v6_tif, i, input_file, /force_lls, /normalize_fractions_after_rms, standardize_band=5
    9: mesma_v6_tif, i, input_file, /force_lls, /normalize_fractions_before_rms, lls_sum_to_one=1., standardize_band=5
    10: mesma_v6_tif, i, input_file, /force_lls, /normalize_fractions_after_rms, lls_sum_to_one=1., standardize_band=5
    11: mesma_v6_tif, i, input_file, /force_lls, lls_sum_to_one=1., /equal_brightness
    12: mesma_v6_tif, i, input_file, /force_lls, /normalize_fractions_before_rms, /equal_brightness
    13: mesma_v6_tif, i, input_file, /force_lls, /normalize_fractions_after_rms, /equal_brightness
    14: mesma_v6_tif, i, input_file, /force_lls, /normalize_fractions_before_rms, lls_sum_to_one=1., /equal_brightness
    15: mesma_v6_tif, i, input_file, /force_lls, /normalize_fractions_after_rms, lls_sum_to_one=1., /equal_brightness
    16: mesma_v6_tif, i, input_file, /force_lls, lls_sum_to_one=0.2, /guerschman
    17: mesma_v6_tif, i, input_file, /force_lls, /normalize_fractions_before_rms, /guerschman
    18: mesma_v6_tif, i, input_file, /force_lls, /normalize_fractions_after_rms, /guerschman
    19: mesma_v6_tif, i, input_file, /force_lls, /normalize_fractions_before_rms, lls_sum_to_one=0.2, /guerschman
    20: mesma_v6_tif, i, input_file, /force_lls, /normalize_fractions_after_rms, lls_sum_to_one=0.2, /guerschman
    21: mesma_v6_tif, i, input_file, /force_lls, standardize_band=5
    22: mesma_v6_tif, i, input_file, /force_lls, /equal_brightness
    23: mesma_v6_tif, i, input_file, /force_lls, /guerschman
    24: mesma_v6_tif, i, input_file
    25: mesma_v6_tif, i, input_file, lls_sum_to_one=1.
    26: mesma_v6_tif, i, input_file, /normalize_fractions_before_rms
    27: mesma_v6_tif, i, input_file, /normalize_fractions_after_rms
    28: mesma_v6_tif, i, input_file, /normalize_fractions_before_rms, lls_sum_to_one=1.
    29: mesma_v6_tif, i, input_file, /normalize_fractions_after_rms, lls_sum_to_one=1.
    30: mesma_v6_tif, i, input_file, lls_sum_to_one=1., standardize_band=5
    31: mesma_v6_tif, i, input_file, /normalize_fractions_before_rms, standardize_band=5
    32: mesma_v6_tif, i, input_file, /normalize_fractions_after_rms, standardize_band=5
    33: mesma_v6_tif, i, input_file, /normalize_fractions_before_rms, lls_sum_to_one=1., standardize_band=5
    34: mesma_v6_tif, i, input_file, /normalize_fractions_after_rms, lls_sum_to_one=1., standardize_band=5
    35: mesma_v6_tif, i, input_file, lls_sum_to_one=1., /equal_brightness
    36: mesma_v6_tif, i, input_file, /normalize_fractions_before_rms, /equal_brightness
    37: mesma_v6_tif, i, input_file, /normalize_fractions_after_rms, /equal_brightness
    38: mesma_v6_tif, i, input_file, /normalize_fractions_before_rms, lls_sum_to_one=1., /equal_brightness
    39: mesma_v6_tif, i, input_file, /normalize_fractions_after_rms, lls_sum_to_one=1., /equal_brightness
    40: mesma_v6_tif, i, input_file, lls_sum_to_one=0.2, /guerschman
    41: mesma_v6_tif, i, input_file, /normalize_fractions_before_rms, /guerschman
    42: mesma_v6_tif, i, input_file, /normalize_fractions_after_rms, /guerschman
    43: mesma_v6_tif, i, input_file, /normalize_fractions_before_rms, lls_sum_to_one=0.2, /guerschman
    44: mesma_v6_tif, i, input_file, /normalize_fractions_after_rms, lls_sum_to_one=0.2, /guerschman
    45: mesma_v6_tif, i, input_file, standardize_band=5
    46: mesma_v6_tif, i, input_file, /equal_brightness
    47: mesma_v6_tif, i, input_file, /guerschman
    48: mesma_v6_tif, i, input_file, /sum_to_one
    49: mesma_v6_tif, i, input_file, /sum_to_one, lls_sum_to_one=1.
    50: mesma_v6_tif, i, input_file, /sum_to_one, /normalize_fractions_before_rms
    51: mesma_v6_tif, i, input_file, /sum_to_one, /normalize_fractions_after_rms
    52: mesma_v6_tif, i, input_file, /sum_to_one, /normalize_fractions_before_rms, lls_sum_to_one=1.
    53: mesma_v6_tif, i, input_file, /sum_to_one, /normalize_fractions_after_rms, lls_sum_to_one=1.
    54: mesma_v6_tif, i, input_file, /sum_to_one, lls_sum_to_one=1., standardize_band=5
    55: mesma_v6_tif, i, input_file, /sum_to_one, /normalize_fractions_before_rms, standardize_band=5
    56: mesma_v6_tif, i, input_file, /sum_to_one, /normalize_fractions_after_rms, standardize_band=5
    57: mesma_v6_tif, i, input_file, /sum_to_one, /normalize_fractions_before_rms, lls_sum_to_one=1., standardize_band=5
    58: mesma_v6_tif, i, input_file, /sum_to_one, /normalize_fractions_after_rms, lls_sum_to_one=1., standardize_band=5
    59: mesma_v6_tif, i, input_file, /sum_to_one, lls_sum_to_one=1., /equal_brightness
    60: mesma_v6_tif, i, input_file, /sum_to_one, /normalize_fractions_before_rms, /equal_brightness
    61: mesma_v6_tif, i, input_file, /sum_to_one, /normalize_fractions_after_rms, /equal_brightness
    62: mesma_v6_tif, i, input_file, /sum_to_one, /normalize_fractions_before_rms, lls_sum_to_one=1., /equal_brightness
    63: mesma_v6_tif, i, input_file, /sum_to_one, /normalize_fractions_after_rms, lls_sum_to_one=1., /equal_brightness
    64: mesma_v6_tif, i, input_file, /sum_to_one, lls_sum_to_one=0.2, /guerschman
    65: mesma_v6_tif, i, input_file, /sum_to_one, /normalize_fractions_before_rms, /guerschman
    66: mesma_v6_tif, i, input_file, /sum_to_one, /normalize_fractions_after_rms, /guerschman
    67: mesma_v6_tif, i, input_file, /sum_to_one, /normalize_fractions_before_rms, lls_sum_to_one=0.2, /guerschman
    68: mesma_v6_tif, i, input_file, /sum_to_one, /normalize_fractions_after_rms, lls_sum_to_one=0.2, /guerschman
    69: mesma_v6_tif, i, input_file, /sum_to_one, standardize_band=5
    70: mesma_v6_tif, i, input_file, /sum_to_one, /equal_brightness
    71: mesma_v6_tif, i, input_file, /sum_to_one, /guerschman

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
  print,'#############################'  
endfor

return
end
