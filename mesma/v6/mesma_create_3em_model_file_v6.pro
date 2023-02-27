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


