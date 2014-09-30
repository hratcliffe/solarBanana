MODULE APPLYBOUNDS

IMPLICIT NONE

CONTAINS

SUBROUTINE electron_bounds(f_new)

 USE CONSTANTS

  USE PARAMETERS

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(-n_v:, :), INTENT(INOUT) :: f_new
  INTEGER(KIND=4):: j

!   open boundary condition

  DO j=3, n_xp+2

    f_new(-n_v, j) = f_new(-n_v+1, j)
    f_new(n_v, j) = f_new(n_v-1, j)
!       infinite velocity boundary condition

    f_new(1, j) = f_new(2, j)
    f_new(-1, j) = f_new(-2, j)
    f_new(0, j) = f_new(1,j)
!  low velocity conditions, enforce zero gradients except at zero
  
  ENDDO

END SUBROUTINE electron_bounds


SUBROUTINE langmuir_bounds(w_new)

 USE CONSTANTS

  USE PARAMETERS

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(-n_kv:, :), INTENT(INOUT) :: w_new
  INTEGER(KIND=4):: j

!   open boundary condition

  DO j=3, n_xp+2

    w_new(-n_kv, j) = -w_new(-n_kv+1, j)
    w_new(n_kv, j) = -w_new(n_kv-1, j)
    w_new(0, j) = 0.d0
!       k=0 bounds

    w_new(1, j) = -w_new(2, j)
    w_new(-1, j) = -w_new(-2, j)

!   k=k_de bounds
  
  ENDDO

END SUBROUTINE langmuir_bounds


SUBROUTINE harmonic_bounds(w_new)

 USE CONSTANTS

  USE PARAMETERS

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(-n_kv:, :), INTENT(INOUT) :: w_new
  INTEGER(KIND=4):: j

!   open boundary condition

  DO j=3, n_xp+2

    w_new(-n_kv, j) = -w_new(-n_kv+1, j)
    w_new(n_kv, j) = -w_new(n_kv-1, j)
!     w_new(0, j) = 0.d0
!       k=0 bounds

!     w_new(1, j) = -w_new(2, j)
!     w_new(-1, j) = -w_new(-2, j)

!   k=k_de bounds
  
  ENDDO

END SUBROUTINE harmonic_bounds

SUBROUTINE fundamental_bounds(w_new)

 USE CONSTANTS

  USE PARAMETERS

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(-n_kv:, :), INTENT(INOUT) :: w_new
  INTEGER(KIND=4):: j

!   open boundary condition

  DO j=3, n_xp+2

    w_new(-n_kv, j) = -w_new(-n_kv+1, j)
    w_new(n_kv, j) = -w_new(n_kv-1, j)
    w_new(0, j) = 0.d0
!       k=0 bounds

    w_new(1, j) = -w_new(2, j)
    w_new(-1, j) = -w_new(-2, j)

!   k=k_de bounds
  
  ENDDO

END SUBROUTINE fundamental_bounds


! SUBROUTINE box_bounds(f_newL, w_newL)
! 
!  USE CONSTANTS
! 
!   USE PARAMETERS
! 
!   IMPLICIT NONE
! 
!   REAL(KIND=8), DIMENSION(-n_kv:, :), INTENT(INOUT) :: w_new
!   INTEGER(KIND=4):: j
! 
! !   open boundary condition
! 
!   DO j=3, n_xp+2
! 
!     w_new(-n_kv, j) = -w_new(-n_kv+1, j)
!     w_new(n_kv, j) = -w_new(n_kv-1, j)
!     w_new(0, j) = 0.d0
! !       k=0 bounds
! 
!     w_new(1, j) = -w_new(2, j)
!     w_new(-1, j) = -w_new(-2, j)
! 
! !   k=k_de bounds
!   
!   ENDDO
! 
! END SUBROUTINE box_bounds
! 









END MODULE APPLYBOUNDS