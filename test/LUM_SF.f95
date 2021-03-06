PROGRAM LUM_SF
  ! This program generates a known bin file to READ by the OSOAA.OUTPUTS.BIN
  ! program and test it

IMPLICIT NONE

INTEGER I,K
INTEGER NT_TOT, NBMU
PARAMETER(NT_TOT=107)
PARAMETER(NBMU=51)

DOUBLE PRECISION I3(0:NT_TOT,-NBMU:NBMU)
DOUBLE PRECISION Q3(0:NT_TOT,-NBMU:NBMU)
DOUBLE PRECISION U3(0:NT_TOT,-NBMU:NBMU)

DO I=0,NT_TOT,1
  DO K=-NBMU,NBMU,1
    I3(I,K) = 0+K+I/108.
    Q3(I,K) = 1000+K+I/108.
    U3(I,K) = 2000+K+I/108.
  END DO
END DO

OPEN(UNIT=101,FILE="OSOAA_RESULTS/Advanced_outputs/LUM_SF.bin",FORM='UNFORMATTED')

WRITE(UNIT=101)((I3(I,K),I=0,NT_TOT),K=-NBMU,NBMU),((Q3(I,K),I=0,NT_TOT),K=-NBMU,NBMU),((U3(I,K),I=0,NT_TOT),K=-NBMU,NBMU)

CLOSE(UNIT=101)

END PROGRAM
