NB. pure_j_lapack.ijs
NB. Implement LAPACK routines on pure J

script_z_ '~system/packages/math/matutil.ijs'  NB. diag

dump=: (3 : 'smoutput ''y'' ; <y') : (4 : 'smoutput 2 2 $ ''x'' ; (<x) ; ''y'' ; (<y)')
prc=: [ (C. " 1) C.                NB. apply permutation x to both rows and columns of table y

NB.                                         ( A11 A12 A13 )
NB. Aperm = P*A*inv(P) = P*A*transpose(P) = ( 0   A22 A23 )
NB.                                         ( 0   0   A33 )
NB. 'start_size permutation'=. gebalpi mat

gebalpi=: 3 : 0
  nz=. 0 & ~:
  nzd=. nz @ diag
  cnzcwod=. (0:`[`(nz @ ({ " 1))) }  NB. count non-zero elms in col x of table y without elm on diagonal
  cnzrwod=. (0:`[`(nz @  {     )) }  NB. count non-zero elms in row x of table y without elm on diagonal
  IOlz=. i: & 0                      NB. IO last zero
  IOfz=. i. & 0                      NB. IO first zero

NB. --- 'orig. A' dump y

  n=. # y
  kl=. 0 , n                         NB. left (up) corner and width (height) of A22
  rcp=. a:                           NB. cyclic permutations stack to separate A33
  ccp=. a:                           NB. cyclic permutations stack to separate A11

  nzcc=. ((+/ - diag) @ nz) y        NB. count non-zeros in cols of square table y without diagonal
  zc=. n $ 0                         NB. zeros correction
NB. --- 'nzcc' dump nzcc
NB. --- 'zc' dump zc
NB. --- 'ccp' dump ccp
  for_c. ({. + (i. @ {:)) kl do.     NB. traverse columns from left to right
NB. --- 'c' dump c
    nzcc=. ccp C. (nzcc - zc)        NB. correct zeros and shift them to end
NB. --- 'new nzcc' dump nzcc
NB. --- 'limited nzcc' dump (,. kl) (] ;. 0) nzcc
    zi=. ({. kl) + (,. kl) (IOfz ;. 0) nzcc    NB. find leftmost zero's index
NB. --- 'zi' dump zi
    if. zi >: {: kl do.              NB. no zeros left?
NB. --- 'break, kl' dump kl
      break.
    end.
    zc=. zc + zi cnzrwod y           NB. correct zeros count ин shifted out
NB. --- 'new zc' dump zc
    ccp=. (< zi , c) , ccp           NB. save cyclic permutation
NB. --- 'new ccp' dump ccp
NB. --- 'new A' dump (rcp , ccp) prc y
    kl=. kl + 1 _1                   NB. exclude leading column
NB. --- 'new kl' dump kl
  end.
  y=. ccp prc y
NB. --- 'curr. A' dump y

  nzcr=. ((+/ " 1 - diag) @ nz) y    NB. count non-zeros in rows of square table y without diagonal
  zc=. n $ 0                         NB. zeros correction
NB. --- 'nzcr' dump nzcr
NB. --- 'zc' dump zc
NB. --- 'rcp' dump rcp
  for_r. |. ({. + (i. @ {:)) kl do.  NB. traverse rows from bottom to up
NB. --- 'r' dump r
    nzcr=. rcp C. (nzcr - zc)        NB. correct zeros and shift them to end
NB. --- 'new nzcr' dump nzcr
NB. --- 'limited nzcr' dump (,. kl) (] ;. 0) nzcr
    zi=. ({. kl) + (,. kl) (IOlz ;. 0) nzcr    NB. find rightmost zero's index
NB. --- 'zi' dump zi
    if. zi >: {: kl do.              NB. no zeros left?
NB. --- 'break, kl' dump kl
      break.
    end.
    zc=. zc + zi cnzcwod y           NB. correct zeros count by shifted out
NB. --- 'new zc' dump zc
    rcp=. (< zi , r) , rcp           NB. save cyclic permutation
NB. --- 'new rcp' dump rcp
NB. --- 'new A' dump (rcp , ccp) prc y
    kl=. kl + 0 _1                   NB. exclude last row
NB. --- 'new kl' dump kl
  end.
NB. xxx   y=. rcp prc y
NB. --- 'curr. A' dump y

  kl ; (C. rcp , ccp)
)

NB. mat_permuted=. gebalp mat
NB. Permute rows and columns to isolate eogenvalues

NB. B -: P mp A mp invP
NB. invP -: |: P
NB. P=. ??? p ilo ihi
NB. perm has not LAPACK, but J standart permutation format

gebalp=: prc~ ((1 & {::) @ gebalpi)

NB. scale=. [start_size] gebalsi mat
gebalsi=: ($:~ (0 , #)) : gebalsidyad

NB. mat_scaled=. [start_size] gebals mat
gebals=: ] * ((*/ %) @ gebalsi)


NB. 'start_size permutation'=. gebalpi mat

NB. 'start_size permutation scale'=. gebali mat
gebali=. 3 : 0
  'ss p'=. gebalpi y
  ss ; p ; ss gebalsi y <<<<<<<<<<<<<<<<<<<<<,
  ( gebalpi) y
)

((0 & {::)) @ gebalpi


NB. mat_balanced=. gebal mat
NB. Balance square matrix A

NB.    load '/home/jip/j602-user/projects/lapack/gebal.ijs'  NB. permute only
NB.    ts=: 6!:2, 7!:2@]
NB.    a1000f=. dzero_jlapack_ + ? 1000 1000 $ 7
NB.    a1000c=. zzero_jlapack_ + j./ ? 2 1000 1000 $ 7
NB.    5 ts 'gebal_jlapack_ a1000f'
NB. 0.239994 2.51897e7
NB.    5 ts 'gebalp a1000f'
NB. 0.151756 1.68269e7
NB.    (gebalp a1000f) -: (2b1000 gebal_jlapack_ a1000f)
NB. 1
NB.    5 ts 'gebal_jlapack_ a1000c'
NB. 0.570236 4.19459e7
NB.    5 ts 'gebalp a1000c'
NB. 0.245782 3.36287e7
NB.    (gebalp a1000c) -: (2b1000 gebal_jlapack_ a1000c)
NB. 1

gebal=: gebals @ gebalp

