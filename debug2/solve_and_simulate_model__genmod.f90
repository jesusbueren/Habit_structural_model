        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 15 10:49:51 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_AND_SIMULATE_MODEL__genmod
          INTERFACE 
            SUBROUTINE SOLVE_AND_SIMULATE_MODEL(ASSET_DISTRIBUTION,     &
     &AV_VSL,AV_V_INI,P50_DELTA)
              REAL(KIND=8), INTENT(OUT) :: ASSET_DISTRIBUTION(36,3,2,4)
              REAL(KIND=8), INTENT(OUT) :: AV_VSL(3)
              REAL(KIND=8), INTENT(OUT) :: AV_V_INI(1,3,2)
              REAL(KIND=8), INTENT(OUT) :: P50_DELTA(3,2,36)
            END SUBROUTINE SOLVE_AND_SIMULATE_MODEL
          END INTERFACE 
        END MODULE SOLVE_AND_SIMULATE_MODEL__genmod
