        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 15 10:49:49 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE COMPUTE_LABOR_INCOME_REVENUE__genmod
          INTERFACE 
            SUBROUTINE COMPUTE_LABOR_INCOME_REVENUE(JOINT_PR,G,AV_V_INI,&
     &COST_EY,COUNTERFACTUALS)
              REAL(KIND=8), INTENT(IN) :: JOINT_PR(1,3,2)
              REAL(KIND=8), INTENT(OUT) :: G
              REAL(KIND=8), INTENT(IN) :: AV_V_INI(1,3,2)
              REAL(KIND=8), INTENT(IN) :: COST_EY(1,3,2)
              REAL(KIND=8) :: COUNTERFACTUALS(6)
            END SUBROUTINE COMPUTE_LABOR_INCOME_REVENUE
          END INTERFACE 
        END MODULE COMPUTE_LABOR_INCOME_REVENUE__genmod
