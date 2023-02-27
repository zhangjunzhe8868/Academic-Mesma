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
