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
