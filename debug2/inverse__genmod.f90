        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 15 10:49:49 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INVERSE__genmod
          INTERFACE 
            SUBROUTINE INVERSE(A,C,N)
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: A(N,N)
              REAL(KIND=8) :: C(N,N)
            END SUBROUTINE INVERSE
          END INTERFACE 
        END MODULE INVERSE__genmod
