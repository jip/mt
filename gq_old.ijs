coclass 'mt'

gqi=: (0 >. <.@(GQNB %~ (_1 + GQNB - GQNX)&+))M.
gqb=: (<. (GQNB * gqi))M.

ungl2step_old=: ((1 _ ,:~ (0 {:: [) <:@-  #@]) (,;.0) 1 {:: [) (((>:@]  0} *) +@-@{:)@[ ,   larfrcfr) 0 ,.  ]
ung2lstep_old=: ((_ 1 ,:~ (0 {:: [)    -~ c@]) (,;.0) 1 {:: [) (((>:@] _1} *)   -@{.)@[ ,.~ larflnbc) 0 , ~ ]
ung2rstep_old=: ((_ 1 ,:~ (0 {:: [) <:@-  c@]) (,;.0) 1 {:: [) (((>:@]  0} *)   -@{:)@[ ,.  larflnfc) 0 ,   ]
ungr2step_old=: ((1 _ ,:~ (0 {:: [)    -~ #@]) (,;.0) 1 {:: [) (((>:@] _1} *) +@-@{.)@[ , ~ larfrcbr) 0 ,.~ ]

ungl2_old=: ungl2step_old^:(;`(#@])`( idmat      @((,  c)- #@])))
ung2l_old=: ung2lstep_old^:(;`(c@])`((idmat~ -~/)@((,~ #)- c@])))
ung2r_old=: ung2rstep_old^:(;`(c@])`( idmat      @((,~ #)- c@])))
ungr2_old=: ungr2step_old^:(;`(#@])`((idmat~ -~/)@((,  c)- #@])))

unglqstep_old=: ((((GQNB ,  _) ,:~ -&c                           ) (] ;. 0)       [) (((ungl2_old~ #)@[) ,   larfbrcfr) ]) ({."1~ ((-GQNB) - c))
ungqrstep_old=: ((((GQNB ,~ _) ,:~ -&#                           ) (] ;. 0)       [) (((ung2r_old~ c)@[) ,.  larfblnfc) ]) ({.  ~ ((-GQNB) - #))

ungqlstep_old=: ((((GQNB ,~ _) ,:~ (0 {:: [) ((<: GQNB) + -~) c@]) (] ;. 0) 1 {:: [) (((ung2l_old~ c)@[) ,.~ larfblnbc) ]) ({.  ~ (  GQNB  + #))
ungrqstep_old=: ((((GQNB ,  _) ,:~ (0 {:: [) ((<: GQNB) + -~) #@]) (] ;. 0) 1 {:: [) (((ungr2_old~ #)@[) , ~ larfbrcbr) ]) ({."1~ (  GQNB  + c))

unglq_old=: ($:~  0 _1 <./@:+ $) : (}:"1@([ (unglqstep_old^:(]`(gqi@#@])`((- gqb@#) ungl2_old (}.~ 2 #   gqb@#)@])))  tru1        @((<.  0 _1    <./ @:+ $) {.   ])))
ungql_old=: ($:~ _1  0 <./@:+ $) : (}.  @([ (ungqlstep_old^:(;`(gqi@c@])`((- gqb@c) ung2l_old (}.~ 2 # -@gqb@c)@]))) (tru1~ -~/@$)@((<. _1  0 -@(<./)@:+ $) {."1 ])))
ungqr_old=: ($:~ _1  0 <./@:+ $) : (}:  @([ (ungqrstep_old^:(]`(gqi@c@])`((- gqb@c) ung2r_old (}.~ 2 #   gqb@c)@])))  trl1        @((<. _1  0    <./ @:+ $) {."1 ])))
ungrq_old=: ($:~  0 _1 <./@:+ $) : (}."1@([ (ungrqstep_old^:(;`(gqi@#@])`((- gqb@#) ungr2_old (}.~ 2 # -@gqb@#)@]))) (trl1~ -~/@$)@((<.  0 _1 -@(<./)@:+ $) {.   ])))
