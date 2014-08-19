MODULE EMPROP

CONTAINS


SUBROUTINE upwind_fund(w_em, omega, l, delta_x)

  USE CONSTANTS
  USE PARAMETERS

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(-n_kv:, :), INTENT(INOUT) :: w_em
  REAL(KIND=8), DIMENSION(:):: omega, l, delta_x
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE:: w_em_new
  REAL(KIND=8), DIMENSION(-n_kv:n_kv):: omega_emf_om, v_grp, coeff2
  INTEGER(KIND=4)::i, j

  ALLOCATE(w_em_new(-n_kv:n_kv, n_xp+3))

  omega_emF_om = sqrt(1.+v_c**2/v_t**2 *k_f**2)
  
  v_grp = v_c*v_c*k_f/v_t/omega_emF_om
  coeff2 = v_grp*k_f


  w_em_new = w_em

  DO j=3, n_xp+2

    DO i = -n_kv,-1
     w_em_new(i,j) = w_em_new(i,j) - dt*v_grp(i)*(w_em(i,j+1) - w_em(i,j))/delta_x(j)
    ENDDO

    DO i = 1, n_kv
     w_em_new(i,j) = w_em_new(i,j) - dt*v_grp(i)*(w_em(i,j) - w_em(i,j-1))/delta_x(j-1)
    ENDDO

    DO i = -n_kv, -1

      w_em_new(i,j) = w_em_new(i,j) + dt *(1./omega_emf_om(i) + coeff2(i))* &
      max(-l(j), 0.d0)*( w_em(i,j) - &
      w_em(i+1,j))/(k_f(i)-k_f(i+1))
    ENDDO


    DO i = 1, n_kv

      w_em_new(i,j) = w_em_new(i,j) - dt *(1./omega_emf_om(i) + coeff2(i))* &
      max(-l(j), 0.d0)*( w_em(i,j) - &
      w_em(i-1,j))/(k_f(i)-k_f(i-1))
    ENDDO

    
    DO i = -n_kv+1, -1

      w_em_new(i,j) = w_em_new(i,j) - dt *(1./omega_emf_om(i) + coeff2(i))* &
      max(l(j), 0.d0)*( w_em(i,j) - &
      w_em(i-1,j))/(k_f(i)-k_f(i-1))
    ENDDO

      i= -n_kv
!       w_em_new(i,j) = w_em_new(i,j) + dt *(1./omega_emf_om(i) + coeff2(i))* &
! 					    max(l(j), 0.d0)*( w_em(i,j))/(k_f(i+1)-k_f(i))
  


    DO i = 1, n_kv-1

      w_em_new(i,j) = w_em_new(i,j) + dt *(1./omega_emf_om(i) + coeff2(i))* &
					    max(l(j), 0.d0)*( w_em(i+1,j) - &
							  w_em(i,j))/(k_f(i+1)-k_f(i))
    ENDDO
    i=n_kv
!       w_em_new(i,j) = w_em_new(i,j) + dt *(1./omega_emf_om(i) + coeff2(i))* &
! 					    max(l(j), 0.d0)*(- w_em(i,j))/(k_f(i)-k_f(i-1))

  

  ENDDO


  w_em = w_em_new

  DEALLOCATE(w_em_new)

END SUBROUTINE upwind_fund



SUBROUTINE upwind_harm(w_em, omega, l, delta_x)

  USE CONSTANTS
  USE PARAMETERS

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(-n_kv:, :), INTENT(INOUT) :: w_em
  REAL(KIND=8), DIMENSION(:):: omega, l, delta_x
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE:: w_em_new
  REAL(KIND=8), DIMENSION(-n_kv:n_kv):: omega_emH_om, v_grp, coeff2
  INTEGER(KIND=4)::i, j

  ALLOCATE(w_em_new(-n_kv:n_kv, n_xp+3))

  omega_emH_om = sqrt(1.+v_c**2/v_t**2 *k_h**2)
  
  v_grp = v_c*v_c*k_h/v_t/omega_emH_om
  coeff2 = v_grp*k_h

  w_em_new = w_em

  DO j=3, n_xp+2
     w_em_new(:,j) = w_em_new(:,j) - dt*v_grp*(w_em(:,j) - w_em(:,j-1))/delta_x(j-1)


    DO i = -n_kv, n_kv-1

    w_em_new(i,j) = w_em_new(i,j) + dt * (1./omega_emH_om(i) + coeff2(i))* &
					  max(l(j), 0.d0)*( w_em(i+1,j) -  w_em(i,j))/(k_h(i+1)-k_h(i))
    ENDDO

    i=n_kv
    w_em_new(i,j) = w_em_new(i,j) + dt * (1./omega_emH_om(i) + coeff2(i))* &
					  max(l(j), 0.d0)*( -1.*w_em(i,j))/(k_h(i)-k_h(i-1))



    DO i = -n_kv+1, n_kv

      w_em_new(i,j) = w_em_new(i,j) - dt *(1./omega_emH_om(i) + coeff2(i))* &
					    max(-l(j), 0.d0)*( w_em(i,j) - &
							  w_em(i-1,j))/(k_h(i)-k_h(i-1))
    ENDDO

    i=-n_kv
    w_em_new(i,j) = w_em_new(i,j) - dt * (1./omega_emH_om(i) + coeff2(i))* &
					    max(-l(j), 0.d0)*( w_em(i,j))/(k_h(i+1)-k_h(i))

  ENDDO


  w_em = w_em_new

  DEALLOCATE(w_em_new)

END SUBROUTINE upwind_harm



END MODULE EMPROP