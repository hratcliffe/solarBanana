MODULE ELECTRON_SOURCE

! This function is the source function which will inject the electrons over time into
! the simulation.  It uses the characteristic length of d and time of tau with a velocity
! spectral index of alpha.  The velocity bounds are between v_min and v_beam.

  USE CONSTANTS
  USE PARAMETERS
 
   IMPLICIT NONE

   CONTAINS


SUBROUTINE fsource(f_new, x)

 ! Electron source term
 
   IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), INTENT(IN) :: x
  REAL(KIND=8), DIMENSION(-n_v:,:) :: f_new
  REAL(KIND=8) :: a, alpha,tau_c, tau_t, prefac
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: fconst
  INTEGER(KIND=4)::i, j
  
  ALLOCATE(fconst(n_xp+3))
  fconst = 0.d0

  tau_c = 1./(sqrt(pi)*0.5*(tau+tau2))
  tau_t = tau
  if (t > t_0) tau_t = tau2

   alpha = 6.  
   a = n_beam*(alpha-1)/(v_min**(-(alpha-1)) - v_beam**(-(alpha-1)))
   prefac = a * tau_c * exp(-(t-t_0)*(t-t_0)/(tau_t*tau_t)) * dt / k1
  
  DO j=3, n_xp+2
   fconst(j) = prefac*exp(-x(j)*x(j)/(d*d)) 
  ENDDO

  DO j=3, n_xp+2

    DO i=1, n_v
      f_new(i,j) = f_new(i,j) + (fconst(j)/ (velocity(i)**alpha))
    ENDDO

!     DO i=v_beam_i+1, n_v-2
!       f_new(i,j) = f_new(i,j) + (fconst(j) / ((v_beam-dv)**alpha)) *(i - v_beam_i)*1d-5
!     ENDDO
  ENDDO

  DEALLOCATE(fconst)

END SUBROUTINE fsource

END MODULE electron_source
