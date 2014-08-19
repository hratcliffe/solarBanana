MODULE QUASILINEAR_TERMS

!subroutines for quasilinear, collisional and transport evolution of electrons and Langmuir waves

IMPLICIT NONE

CONTAINS

SUBROUTINE quasilinear(f, w, f_new, w_new, coeff2)
! Quasilinear alteration of the electron distribution function 

USE CONSTANTS

USE PARAMETERS

IMPLICIT NONE


REAL(KIND=8), DIMENSION(-n_v:, :) :: f_new,f
REAL(KIND=8), DIMENSION(-n_kv:,:) :: w_new,w
REAL(KIND=8), DIMENSION(:) :: coeff2
! REAL(KIND=8), DIMENSION(-n_v:n_v, n_xp+3) :: m1, m2
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: m1, m2
INTEGER(KIND=4)::i

ALLOCATE(m1(-n_v:n_v,n_xp+3))
ALLOCATE(m2(-n_v:n_v,n_xp+3))


DO i = -n_v+1,-1
  m1(i,:) = w(i,:)*(f(i+1,:)-f(i,:)) / (velocity(i)*dvdv)
  m2(i,:) = w(i-1,:)*(f(i,:)-f(i-1,:)) / (velocity(i-1)*dvdv)
ENDDO

DO i = 1,n_v-1
  m1(i,:) = w(i,:)*(f(i+1,:)-f(i,:)) / (velocity(i)*dvdv)   
  m2(i,:) = w(i-1,:)*(f(i,:)-f(i-1,:)) / (velocity(i-1)*dvdv)	
ENDDO


DO i = -n_v+1,-1
  f_new(i,:) = f_new(i,:) - a1*dt_ql*(m1(i,:)-m2(i,:))
  w_new(i,:) = w_new(i,:) + velocity(i)*velocity(i)*dt_ql/dv*coeff2 * W(i,:) * (f(i,:)-f(i-1,:))
ENDDO
! 
DO i = 1,n_v-1
  f_new(i,:) = f_new(i,:) + a1*dt_ql*(m1(i,:)-m2(i,:))
  w_new(i,:) = w_new(i,:) + velocity(i)*velocity(i)*dt_ql/dv*coeff2 * W(i,:)* (f(i+1,:)-f(i,:))
ENDDO


DEALLOCATE(m1, m2)


END SUBROUTINE quasilinear


SUBROUTINE landau_background(w, w_new, landau_arr)
!   resonant interaction of L waves with thermal background electrons
USE CONSTANTS
USE PARAMETERS

IMPLICIT NONE

REAL(KIND=8), DIMENSION(-n_kv:,:) :: w_new,w
REAL(KIND=8), DIMENSION(-n_v:,:), INTENT(IN) ::landau_arr


w_new(-n_v:n_v,:) = w_new(-n_v:n_v,:) - w(-n_v:n_v,:)*dt_ql*Landau_arr(-n_v:n_v,:)


END SUBROUTINE landau_background


SUBROUTINE inhomogeneity(w, w_new, l)
! calculation of the velocity diffusion term 

USE CONSTANTS

USE PARAMETERS

IMPLICIT NONE

REAL(KIND=8), DIMENSION(-n_kv:, :) :: w
REAL(KIND=8), DIMENSION(-n_kv:, :) :: w_new
REAL(KIND=8), DIMENSION(n_xp+3) :: l  
INTEGER(KIND=4)::i
REAL(KIND=8):: const

  const=v_t*dt_ql

  DO i= 1, n_kv-1
    w_new(i,:) = w_new(i,:) - const*max(l, 0.d0) * (W(i,:) - W(i-1,:))/dk_x(i-1)
    w_new(i,:) = w_new(i,:) + const*max(-l, 0.d0) * (W(i+1,:) - W(i,:))/dk_x(i)

  ENDDO

  DO i= -n_kv +1, -1
    w_new(i,:) = w_new(i,:) + const*max(l, 0.d0) * (W(i+1,:) - W(i,:))/dk_x(i)
    w_new(i,:) = w_new(i,:) - const*max(-l, 0.d0) * (W(i,:) - W(i-1,:))/dk_x(i-1)

  ENDDO


END SUBROUTINE inhomogeneity

 
SUBROUTINE collision(f,f_new,w_new,gamma_xp,gamma_xw,w,w_t)
! Collisional Term of the electron distribution function 

USE CONSTANTS

USE PARAMETERS

IMPLICIT NONE

REAL(KIND=8), DIMENSION(-n_v:,:), INTENT(IN) :: f
REAL(KIND=8), DIMENSION(-n_v:,:), INTENT(OUT) :: f_new
REAL(KIND=8), DIMENSION(-n_kv:,:), INTENT(IN) :: w, w_t
REAL(KIND=8), DIMENSION(-n_kv:,:), INTENT(OUT) :: w_new
REAL(KIND=8), DIMENSION(:), INTENT(IN) :: gamma_xp,gamma_xw
INTEGER(KIND=4)::i

DO i=1, n_v-1
  f_new(i,:) = f_new(i,:) + dt_ql*gamma_xp(:)/dv * &
		      (f(i+1,:)/(velocity(i+1)*velocity(i+1)) &
			- f(i,:)/(velocity(i)*velocity(i)))
ENDDO

DO i=-n_v+1, -1
  f_new(i,:) = f_new(i,:) - dt_ql*gamma_xp(:)/dv * (f(i,:)/(velocity(i)*velocity(i)) &
			- f(i-1,:)/(velocity(i-1)*velocity(i-1)))
ENDDO

DO i=-n_kv,n_kv
  w_new(i,:) = w_new(i,:) - dt_ql*gamma_xw(:)*(w(i,:)-w_t(i,:))
ENDDO


END SUBROUTINE collision



SUBROUTINE spontaneous(f, w_new, log_v_v_t ,coeff)
! Spontaneous Term of the spectral energy density 

USE CONSTANTS

USE PARAMETERS

IMPLICIT NONE

REAL(KIND=8), DIMENSION(-n_v:,:), INTENT(IN) :: f
REAL(KIND=8), DIMENSION(-n_v:), INTENT(IN) :: log_v_v_t
REAL(KIND=8), DIMENSION(-n_kv:,:), INTENT(OUT) :: w_new
REAL(KIND=8), DIMENSION(:) :: coeff
INTEGER(KIND=4)::i

DO i = -n_v,n_v
	w_new(i,:) = w_new(i,:) + dt_ql*coeff(:)*velocity(i)*f(i,:)*log_v_v_t(i)
ENDDO


END SUBROUTINE spontaneous



SUBROUTINE radial(f,f_new,r_x)
! Radial Term of the electron distribution function 

USE CONSTANTS

USE PARAMETERS

IMPLICIT NONE

REAL(KIND=8), DIMENSION(-n_v:, :), INTENT(IN) :: f
REAL(KIND=8), DIMENSION(-n_v:, :), INTENT(OUT) :: f_new
REAL(KIND=8), DIMENSION(:) :: r_x
INTEGER(KIND=4) :: i

DO i = 3, n_xp+2

  f_new(:,i) = f_new(:,i) - dt *2.*velocity(:)*f(:,i)/( r_x(i) + 3.4d9 )

ENDDO

END SUBROUTINE radial


SUBROUTINE vanleer(f,f_new,delta_x)
! Van Leer Transport Term 

USE CONSTANTS

USE PARAMETERS

IMPLICIT NONE

REAL(KIND=8), DIMENSION(-n_v:,:), INTENT(IN) :: f
REAL(KIND=8), DIMENSION(-n_v:n_v, n_xp+3), INTENT(OUT) :: f_new
REAL(KIND=8), DIMENSION(:), INTENT(IN) :: delta_x
REAL(KIND=8), DIMENSION(-n_v:n_v) :: d1 = 0.d0, DF1 = 0.d0, d2 = 0.d0, DF2 = 0.d0
REAL(KIND=8), DIMENSION(-n_v:n_v) :: aa_minus = 0.0, aa_plus = 0.0, denom
INTEGER(KIND=4)::i,j


DO j = 3, n_xp+2

    d1 = ( f(:,j) - f(:,j-1) )*( f(:,j+1) - f(:,j) )
    denom(:)=max(abs(f(:, j+1) - f(:,j-1)), tiny(0.d0))
    DF1 = (delta_x(j+1) + delta_x(j))*( f(:,j+1) - f(:,j-1) )/denom(:)/(delta_x(j)*denom(:))* MAX(0.d0,d1(:))

 ! faster method to check no zero rather than using if statement... Max is intrinsic function
  !added protection from divide by zero error 

! IF ( d1 > 0.) THEN 	
! 		DF1 = d1*(delta_x1+delta_x2)/(delta_x1*(f1(i)-f3(i)))
! 	ELSE 			
! 		DF1 = 0.
! 	END IF

    d2 = ( f(:,j-1) - f(:,j-2) )*( f(:,j) - f(:,j-1) )
    denom(:)=max(abs(f(:,j)-f(:,j-2)), tiny(0.d0))
    DF2 = (delta_x(j) + delta_x(j-1))*( f(:,j) - f(:,j-2) )/denom(:)/(delta_x(j)*denom(:)) * MAX(0.d0,d2(:))

! IF ( d2 > 0.) THEN	
! 		DF2 = d2*(delta_x2+delta_x3)/(delta_x2*(f2(i)-f4(i)))
! 	ELSE  			
! 		DF2 = 0.
! 	END IF

  aa_minus = velocity*dt/delta_x(j)
  aa_plus = velocity*dt/delta_x(j+1)

  f_new(:,j) = f_new(:,j) - aa_minus(:)*(f(:,j) - f(:,j-1) + (1.-aa_plus(:))*Df1(:)/2. - (1.-aa_minus(:))*Df2(:)/2.)


ENDDO

END SUBROUTINE vanleer

SUBROUTINE upwind(w, w_new, delta_x)
! Upwind Transport Term 

USE CONSTANTS
USE PARAMETERS

IMPLICIT NONE

REAL(KIND=8), DIMENSION(-n_kv:, :), INTENT(IN) :: w
REAL(KIND=8), DIMENSION(-n_kv:, :), INTENT(OUT) :: w_new
REAL(KIND=8), DIMENSION(:), INTENT(IN) :: delta_x
INTEGER(KIND=4)::i,j

DO j=3, n_xp+2

  DO i=1, n_kv
    w_new(i,j) = w_new(i,j) - 3.*v_t*k_x(i)*dt/max(delta_x(j),tiny(0.d0)) * ( w(i,j) - w(i,j-1) )
  ENDDO
  
  DO i=-n_kv, -1
    w_new(i,j) = w_new(i,j) - 3.*v_t*k_x(i)*dt/max(delta_x(j+1),tiny(0.d0)) * ( w(i,j+1) - w(i,j) )
  ENDDO

ENDDO

END SUBROUTINE upwind


END MODULE QUASILINEAR_TERMS
