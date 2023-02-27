;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; guerschman_transform_v5 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function guerschman_transform_v5, mcd43
  ;This function converts either a spectrum or an 3-dimensional spectral image
  ; (with band as the 3rd dimension, bsq) to the transform given by guerschman et al 2015, RSE

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

  return, out
end
