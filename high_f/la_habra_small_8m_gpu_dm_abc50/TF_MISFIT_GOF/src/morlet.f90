FUNCTION MORLET ( S, FA )

  USE GLOBAL, ONLY: WP, PI, W0

  IMPLICIT NONE 

  REAL(WP), PARAMETER  :: PI14 =  0.7511255444649
  REAL(WP), INTENT(IN) :: S, FA
  COMPLEX              :: MORLET
 
  IF (FA == 0.) THEN 
    MORLET = (0.,0.)
  ELSE
    MORLET = CMPLX(PI14*EXP( (-(S*2.*PI*FA-W0)**2) / 2.),0.)
  ENDIF

END FUNCTION MORLET
