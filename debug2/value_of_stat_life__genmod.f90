        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 15 10:49:50 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE VALUE_OF_STAT_LIFE__genmod
          INTERFACE 
            FUNCTION VALUE_OF_STAT_LIFE(A_TODAY,A_PRIME,FV,T_L,H_L,PI_L,&
     &E_L,Y_L,DF_L)
              REAL(KIND=8), INTENT(IN) :: A_TODAY
              REAL(KIND=8), INTENT(IN) :: A_PRIME
              REAL(KIND=8), INTENT(IN) :: FV(200,3,10)
              INTEGER(KIND=4), INTENT(IN) :: T_L
              INTEGER(KIND=4), INTENT(IN) :: H_L
              INTEGER(KIND=4), INTENT(IN) :: PI_L
              INTEGER(KIND=4), INTENT(IN) :: E_L
              INTEGER(KIND=4), INTENT(IN) :: Y_L
              INTEGER(KIND=4), INTENT(IN) :: DF_L
              REAL(KIND=8) :: VALUE_OF_STAT_LIFE
            END FUNCTION VALUE_OF_STAT_LIFE
          END INTERFACE 
        END MODULE VALUE_OF_STAT_LIFE__genmod
