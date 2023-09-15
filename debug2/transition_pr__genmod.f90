        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 15 10:49:48 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TRANSITION_PR__genmod
          INTERFACE 
            SUBROUTINE TRANSITION_PR(RHO_AV,S2_NU_AV,NZZ,MAG_1,MAG_2,   &
     &PR_P_NEW)
              INTEGER(KIND=4), INTENT(IN) :: NZZ
              REAL(KIND=8), INTENT(IN) :: RHO_AV
              REAL(KIND=8), INTENT(IN) :: S2_NU_AV
              REAL(KIND=8), INTENT(IN) :: MAG_1(NZZ,1)
              REAL(KIND=8), INTENT(IN) :: MAG_2(NZZ,1)
              REAL(KIND=8), INTENT(OUT) :: PR_P_NEW(NZZ,NZZ)
            END SUBROUTINE TRANSITION_PR
          END INTERFACE 
        END MODULE TRANSITION_PR__genmod
