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

  out(1)=sse
  out(2:*)=fracs

  return, out
end