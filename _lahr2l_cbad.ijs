lahr2l=: 3 : 0
  Y=. V=. 0 {. y
  T=. H=. 0 0 $ 0
  for_j. i. HRDBS do.
    b=. IOSFR { y
    y=. 1 0 }. y
    b=. b (- dbg 'L1 -') (+ (< a: ; <: j) { V) (mp dbg 'L1 mp') Y
    b=. b (- dbg 'L2 -') (((0 (_1) } b) (mp dbg 'L2 mp1') (ct V)) (mp dbg 'L2 mp2') T) (mp dbg 'L2 mp3') V  NB. matrix-by-vector ops only
    z1=. 1 (0) } z=. (larfgfc dbg 'L3 larfgfc') j }. b
    'beta tau'=. 0 _1 { z
    u=. z1 (* dbg 'L5 *') + - tau
    w=. (0 (_1) } u) (mp dbg 'L6 mp') (ct ((0 , j) }. V))
    T=. T ((0 append) dbg 'L7 0append T') ((w (mp dbg 'L7 mp') (ct T)) , tau)
    Y=. Y (, dbg 'L8 , Y') ((w (mp dbg 'L8 mp2') Y) - (u (mp dbg 'L8 mp1') y))
    H=. H (0 append) ((j {. b) , beta)
    V=. V (_1 append) z1
  end.
  Y ; V ; H ; T
)
