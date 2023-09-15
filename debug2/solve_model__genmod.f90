        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 15 10:49:50 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_MODEL__genmod
          INTERFACE 
            SUBROUTINE SOLVE_MODEL(A_POLICY,VSL)
              REAL(KIND=8), INTENT(OUT) :: A_POLICY(200,36,2,10,3,2,1)
              REAL(KIND=8), INTENT(OUT) :: VSL(200,2,10,3,2,1)
            END SUBROUTINE SOLVE_MODEL
          END INTERFACE 
        END MODULE SOLVE_MODEL__genmod
