;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; constrained_unmix_v6 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function sum_and_sse, fracs
  common shared, em_array, target
  sum=total(fracs)
  sse=total((float(em_array)##fracs-float(target))^2.) ;Sum of the square of errors (sse)
  return, [sum,sse]
end

function constrained_unmix_v6, em_array_in, target_in ;sum_to_one=sum_to_one
  ;em_array is a n_em x n_band array, target is a 1 x n_band
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

  frac_lower_bnd=fltarr(n_em)
  frac_upper_bnd=frac_lower_bnd+1.
  frac_bnd=[[frac_lower_bnd],[frac_upper_bnd]]
  sum_and_sse_bnd=[[float(keyword_set(sum_to_one)),0.],[1.,!values.f_infinity]]

  nobj=1 ;I think this is the index in what's return by sum_and_sse for what we want to minimize, here sse
  func='sum_and_sse'
  title='IDL constrained min report
  report='constrained_min_report.txt'

  constrained_min, fracs, frac_bnd, sum_and_sse_bnd, nobj, func, inform; report=report, title=title;, epstop=1e-6

  return, [fix(inform),total((float(em_array)##fracs-float(target))^2.),reform(fracs)]
end
