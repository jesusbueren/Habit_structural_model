        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 15 10:49:49 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SORT__genmod
          INTERFACE 
            SUBROUTINE SORT(ARR,LENGTH)
              INTEGER(KIND=4), INTENT(IN) :: LENGTH
              REAL(KIND=8), INTENT(INOUT) :: ARR(LENGTH)
            END SUBROUTINE SORT
          END INTERFACE 
        END MODULE SORT__genmod
