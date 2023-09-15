        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 15 10:49:51 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SIMULATE_MODEL__genmod
          INTERFACE 
            SUBROUTINE SIMULATE_MODEL(A_POLICY,VSL,ASSET_DISTRIBUTION,  &
     &AV_VSL,AV_V_INI,P50_DELTA,LAMBDA_REF)
              REAL(KIND=8), INTENT(IN) :: A_POLICY(200,36,2,10,3,2,1)
              REAL(KIND=8), INTENT(IN) :: VSL(200,2,10,3,2,1)
              REAL(KIND=8), INTENT(OUT) :: ASSET_DISTRIBUTION(36,3,2,4)
              REAL(KIND=8), INTENT(OUT) :: AV_VSL(3)
              REAL(KIND=8), INTENT(OUT) :: AV_V_INI(1,3,2)
              REAL(KIND=8), INTENT(OUT) :: P50_DELTA(3,2,36)
              REAL(KIND=8) ,OPTIONAL :: LAMBDA_REF(3,2)
            END SUBROUTINE SIMULATE_MODEL
          END INTERFACE 
        END MODULE SIMULATE_MODEL__genmod
