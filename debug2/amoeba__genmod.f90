        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 15 10:49:50 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE AMOEBA__genmod
          INTERFACE 
            SUBROUTINE AMOEBA(P,Y,FTOL,FUNC,ITER)
              REAL(KIND=8), INTENT(INOUT) :: P(:,:)
              REAL(KIND=8), INTENT(INOUT) :: Y(:)
              REAL(KIND=8), INTENT(IN) :: FTOL
              INTERFACE 
                FUNCTION FUNC(X)
                  REAL(KIND=8), INTENT(IN) :: X(:)
                  REAL(KIND=8) :: FUNC
                END FUNCTION FUNC
              END INTERFACE 
              INTEGER(KIND=4), INTENT(OUT) :: ITER
            END SUBROUTINE AMOEBA
          END INTERFACE 
        END MODULE AMOEBA__genmod
