MODULE DISTRIBSINIT
 
! define initial electron and wave distributions

USE CONSTANTS
USE PARAMETERS
 
IMPLICIT NONE

CONTAINS


SUBROUTINE distrib_init(fvx, wvx, w_s, w_t, w_em, w_em_f, omega_x, r_x)
!   sets up initial spectral energy density of Langmuir waves and IS waves
!   also calculates total initial wave energy

USE CONSTANTS
USE PARAMETERS
USE MPI

IMPLICIT NONE

REAL(KIND=8), DIMENSION(-n_v:,:), INTENT(OUT) :: fvx
REAL(KIND=8), DIMENSION(-n_kv:,:), INTENT(OUT) :: wvx, w_s, w_t, w_em, w_em_f
REAL(KIND=8), DIMENSION(:), INTENT(IN) :: omega_x, r_x
INTEGER(KIND=4):: i,j


 DO i=-n_v, n_v
    fvx(i,:)=0.
    v = velocity(i)
    wvx(i,:) = (k_b*T_e*omega_x*omega_x/(4.*pi*pi*v*v*k2))*log(abs(v)/v_t)
    w_s(i,:) = 1d-20*abs(k_x(i))
 ENDDO

 DO i=n_v+1,n_kv
	wvx(i,:) = (k_b*T_e*omega_x*omega_x/(4.*pi*pi*v_t*v_t*k2))*k_x(i)*k_x(i)*log(1./abs(k_x(i)))
	wvx(-i,:) = wvx(i,:)
	w_s(i,:) = 1d-20*abs(k_x(i))
	w_s(-i,:) = 1d-20*abs(k_x(i))
 ENDDO


!----------- Boundary conditions -------------------------------------

  wvx(0,:) = 0.
  w_s(0,:) = 0.
! Makes sure the infinite velocity term is initialised to zero

  w_t=wvx

  ew_0 = 0.

DO j=3, n_xp+2
  DO i = 1, n_v
    ew_0 = omega_x(j) * dv/(velocity(i)*velocity(i)) * w_t(i,j) *k2 + ew_0	
  ENDDO

  DO i = n_v+1, n_kv
    ew_0 = omega_x(j) / v_t *k2 *dk * w_t(i,j) + ew_0
  ENDDO
ENDDO

  ew_0 = ew_0*2.
! Take into account the negative side which is equal to the positive side
   
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, ew_0, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ecode)

!---------------------------------------------------------------------



  DO i=-n_kv, n_kv
    w_em(i,:) = kbTe_1d *k_h(i)*k_h(i)/(v_t*v_t)*omega_x*omega_x
    w_em_f(i,:) = kbTe_1d *k_f(i)*k_f(i)/(v_t*v_t)*omega_x*omega_x
  ENDDO



END SUBROUTINE distrib_init

END MODULE DISTRIBSINIT
