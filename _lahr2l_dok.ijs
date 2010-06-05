lahr2l=: 3 : 0
  Y=. V=. 0 {. y
  T=. H=. 0 0 $ 0
  for_j. i. HRDBS do.
    b=. IOSFR { y
    y=. 1 0 }. y
    b=. b - (+ (< a: ; <: j) { V) mp Y
    b=. b - (((0 (_1) } b) mp (ct V)) mp (ct T)) mp V  NB. matrix-by-vector ops only
    z1=. 1 (0) } z=. larfgf j }. b
    'beta tau'=. 0 _1 { z
    u=. z1 * - tau
    w=. + (((0 , j) }. V) mp (+ 0 (_1) } u))
    T=. T (0 append) ((w mp T) , tau)
    Y=. Y , ((w mp Y) - (u mp y))
    H=. H (0 append) ((j {. b) , beta)
    V=. V (_1 append) z1
  end.
  Y ; V ; H ; T
)
