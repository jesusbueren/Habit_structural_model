        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 15 10:49:48 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DISCRETIZE_IID_SHOCKS__genmod
          INTERFACE 
            SUBROUTINE DISCRETIZE_IID_SHOCKS(S2_NU_AV,NZZ,MAG_P,PR0_P)
              INTEGER(KIND=4), INTENT(IN) :: NZZ
              REAL(KIND=8), INTENT(IN) :: S2_NU_AV
              REAL(KIND=8), INTENT(OUT) :: MAG_P(NZZ,1)
              REAL(KIND=8), INTENT(OUT) :: PR0_P(NZZ,1)
            END SUBROUTINE DISCRETIZE_IID_SHOCKS
          END INTERFACE 
        END MODULE DISCRETIZE_IID_SHOCKS__genmod
