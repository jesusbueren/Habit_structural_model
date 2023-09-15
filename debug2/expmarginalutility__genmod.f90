        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 15 10:49:50 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EXPMARGINALUTILITY__genmod
          INTERFACE 
            FUNCTION EXPMARGINALUTILITY(A,FSAVINGS,T_L,H_L,PI_L,E_L,Y_L,&
     &DF_L)
              REAL(KIND=8), INTENT(IN) :: A
              REAL(KIND=8), INTENT(IN) :: FSAVINGS(200,2,10)
              INTEGER(KIND=4), INTENT(IN) :: T_L
              INTEGER(KIND=4), INTENT(IN) :: H_L
              INTEGER(KIND=4), INTENT(IN) :: PI_L
              INTEGER(KIND=4), INTENT(IN) :: E_L
              INTEGER(KIND=4), INTENT(IN) :: Y_L
              INTEGER(KIND=4), INTENT(IN) :: DF_L
              REAL(KIND=8) :: EXPMARGINALUTILITY
            END FUNCTION EXPMARGINALUTILITY
          END INTERFACE 
        END MODULE EXPMARGINALUTILITY__genmod
