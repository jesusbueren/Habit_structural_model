        !COMPILER-GENERATED INTERFACE MODULE: Fri Sep 15 10:49:51 2023
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TRANSITIONS__genmod
          INTERFACE 
            SUBROUTINE TRANSITIONS(BETA_H,BETA_D,H,LE,JOINT_YH)
              REAL(KIND=8), INTENT(IN) :: BETA_H(4,2,2,3)
              REAL(KIND=8), INTENT(IN) :: BETA_D(4,2,2,3)
              REAL(KIND=8), INTENT(OUT) :: H(3,3,37,2,2,3)
              REAL(KIND=8), INTENT(OUT) :: LE(2,2,3,3)
              REAL(KIND=8), INTENT(IN) :: JOINT_YH(37,2,2,3,2,5)
            END SUBROUTINE TRANSITIONS
          END INTERFACE 
        END MODULE TRANSITIONS__genmod
