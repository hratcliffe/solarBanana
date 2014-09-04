
MODULE INITIALISE

!setup code. define space and velocity grids, coronal density profile

IMPLICIT NONE

CONTAINS

SUBROUTINE num_points_check

  USE CONSTANTS
  USE PARAMETERS

  IMPLICIT NONE

  INTEGER:: n_x_old

  IF(n_xp .LE. 25) print*, 'WARNING, too few points per core. Reduce n_p'
  IF(n_xp*n_p .NE. n_x-3) THEN
!    print*, 'WARNING, n_x does not divide n_p, terminating'
!    stop
    !terminate code OR next line instead increases n_x
  n_x_old = n_x
  n_x = (n_xp+1)*n_p + 3
  print*, 'WARNING, n_x does not divide n_p.............................................'
  print*, 'changing n_x from', n_x_old, 'to', n_x

  ENDIF

END SUBROUTINE


SUBROUTINE space_grid_define(xmin, r_x, delta_x)

  USE CONSTANTS
  USE PARAMETERS

  IMPLICIT NONE

  REAL(KIND=8), INTENT(IN) ::xmin
  REAL(KIND=8), DIMENSION(:),INTENT(OUT)::delta_x, r_x
  REAL(KIND=8):: x_local, d_x, dtemp, travel
  INTEGER(KIND=4)::i

  d_x = d/7.
  dtemp = d
  x_local = xmin;

  DO i=1, n_x
    
    travel = d_x/v_beam
    dtemp = dtemp + dv*travel*0.225
    d_x = dtemp/10.
    r_x(i) = x_local
      
    IF (i .LE. 15040) then
      d_x = d/6.
      dtemp = d
    ENDIF

    x_local = x_local+d_x

  ENDDO
  
  DO i=2, n_x

    delta_x(i) = r_x(i) - r_x(i-1)

  ENDDO

  delta_x(1) = delta_x(2)

END SUBROUTINE space_grid_define

SUBROUTINE velocity_define

  USE CONSTANTS
  USE PARAMETERS

  IMPLICIT NONE

  INTEGER, DIMENSION(1):: min_tmp
  REAL(KIND=8), DIMENSION(40)::tmp
  INTEGER(KIND=4)::i


  v_ti = sqrt (1.*k_b*t_i/m_i)
  v_t = sqrt (1.*k_b*t_e/m_e)
  three_v_t_sq = 3.0*v_t*v_t
  v_s = sqrt(1.*k_b*t_e/m_i * (1. + 3.*t_i/t_e))
  ! Thermal Velocities and Timesaving constant

! Define velocity and wavenumber grids and steps ........................................................

  dv = (v_max - v_min)/(n_v-1);
  dvdv = dv*dv

  DO i= 1, n_v  
    velocity(i) = v_min + (i-1)*dv
    velocity(-i) = -velocity(i)
  ENDDO


  k_x(-n_v:n_v) = v_t/velocity(-n_v:n_v)

  dk_x(1:n_v-1) = abs(k_x(2:n_v) -  k_x(1:n_v-1))
  dk_x(-n_v+1:-1) = abs(k_x(-n_v:-2) - k_x(-n_v+1:-1)) 

  dk=dk_x(n_v-1)

  
  DO i= n_v+1, n_kv
    k_x(i) = (v_t/(v_max+dv))*float(n_kv-i+1)/(n_k)
    k_x(-i) = -k_x(i)
  ENDDO

!   DO i= n_v+1, n_kv
!     k_x(i) = k_x(n_v) - dk*0.98*float(i-n_v)
!     k_x(-i) = -k_x(i)
!   ENDDO
! 
!   k_x(n_v+1)=(k_x(n_v)+k_x(n_v+2))/2.
! 
! print*, 'WARNING k_x is BAD!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

  velocity(0) = 1.
  k_x(0) = 0.

  
!   wavenumber step in equally spaced part

  dk_x(1:n_kv-1) = abs(k_x(2:n_kv) -  k_x(1:n_kv-1))
  dk_x(-n_kv+1:-1) = abs(k_x(-n_kv:-2) - k_x(-n_kv+1:-1)) 
  dk_x(n_kv) = dk_x(n_kv-1)

  dk_x(-n_kv) = dk_x(-n_kv+1)
  dk_x(0) = dk_x(1)
!   dk = dk_x(n_kv)

! tmp(0:20)=
! if(rank ==0) print*,  'dk',n_kv, dk, dk_x, k_x


! stop

! ----------------------------------------------------------------------------------------------

!   find index of velocity=v_beam
  min_tmp=minloc(abs(velocity - v_beam))
  v_beam_i=min_tmp(1) - n_v

END SUBROUTINE velocity_define

SUBROUTINE em_wavevec_define

  USE CONSTANTS
  USE PARAMETERS

  IMPLICIT NONE

  INTEGER, DIMENSION(1):: min_tmp
  REAL(KIND=8) :: k_h_rng, k_f_max, dk_f,k_h_max,k_h_min, k_f_min
  INTEGER(KIND=4)::i

!   pick max wanted frequency of fund emm treatment 
!   make it just beyond resonant region w. L waves... say 10% further.

  k_f_max = 1.05*sqrt(3.)*v_t*maxval(k_x)/v_c
  dk_f = k_f_max/n_kv
  k_f_min = dk_f*0.

  DO i=1, n_kv
    k_f(i) = k_f_min + dk_f*float(i)
    k_f(-i) = -k_f(i)
  ENDDO

  k_f(0)=0.


  !   pick max wanted frequency of harm emm treatment 
!   make it in terms of omega_pe say 10 percent total range

  k_h_min = sqrt(1.9*1.9 - 1.)*v_t/v_c
  k_h_max = sqrt(2.3*2.3 - 1.)*v_t/v_c
  dk_h = (k_h_max - k_h_min)/(2.*float(n_kv))

  DO i=-n_kv, -1
!     k_h(i)=sqrt(3.)/v_c*v_t  + float(i+1)/float(n_kv)/200.      
     k_h(i) = k_h_min + float(i+1+n_kv)*dk_h
  ENDDO

  DO i=0, n_kv
!     k_h(i)=sqrt(3.)/v_c*v_t  + float(i+1)/float(n_kv)/200.
    k_h(i) = k_h_min + float(i+1+n_kv)*dk_h
  ENDDO


END SUBROUTINE em_wavevec_define


SUBROUTINE density_define(Density_x, density_smooth_x, Omega_x, Wind_x, R_x, dn_dx, delta_x, l)

  USE CONSTANTS
  USE PARAMETERS
  USE MPI

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:),INTENT(OUT)::Density_x,Density_smooth_x,Omega_x,Wind_x,dn_dx,l;
  REAL(KIND=8), DIMENSION(:),INTENT(IN):: R_x, delta_x
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE::density2_x, r_corr, r_corr2
  REAL(KIND=8):: si, dn0_dx, au_pos,total, temp;
  REAL(KIND=8):: log_rho_min, log_rho_max, phi_min, phi_max, beta, density_1AU
  REAL(KIND=8), DIMENSION(N):: alpha, rho, phi
  INTEGER(KIND=4)::i,j


  ALLOCATE(density2_x(n_xp+3))
  ALLOCATE(r_corr(n_xp+3))
  ALLOCATE(r_corr2(n_xp+3))

  density_x=0.
  density_smooth_x=0.
  omega_x=0.
  wind_x=0.
  dn_dx=0.
  l=0.

! define smooth density profile ...................................................................
 
   DO i=1, n_xp+3
  
     CALL corona(density_x(i), wind_x(i), r_x(i) + r_s + x_0)
   
   ENDDO

   density_smooth_x=density_x

! define randomly fluctuating density and add to smooth profile
!  number of components is N

! parameters of flcutuations
  si = 0.8
  log_rho_min = 7.
  log_rho_max = 10.
  phi_min = 0.
  phi_max = 2.*pi
  beta = -0.25

  density_1AU = 4.

  r_corr = (density_x/density_1AU)**beta
  r_corr2 = (density_x/density_1AU)**(beta-1.)

  DO i=1,N 
    CALL random_number(temp)
    CALL MPI_BCAST(temp, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ecode)
!     make sure random number is same across processors...

    rho(i) = 10.**(log_rho_min + temp*(log_rho_max-log_rho_min))
    alpha(i) = (rho(i)/2d12)**(si)*5.
    CALL random_number(temp)
    CALL MPI_BCAST(temp, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ecode)
    phi(i) = phi_min + temp*(phi_max-phi_min)
  ENDDO

  density2_x = density_x

  DO i=1,N
 	density2_x(1) = density2_x(1) + 1*density_x(1)*alpha(i)*sin(2.*pi*r_x(1)/rho(i) + phi(i))
  ENDDO


  DO i=2,n_xp+3


    dn_dx(i-1) = (density_x(i) - density_x(i-1))/delta_x(i) 
    dn0_dx = dn_dx(i-1)
    
    DO j=1,N
      density2_x(i) = density2_x(i) + density_x(i)*r_corr(i)*alpha(j)*sin(2.*pi*r_x(i)/rho(j) + phi(j))
      dn_dx(i-1) = dn_dx(i-1) + dn0_dx*r_corr(i)*alpha(j)*sin(2.*pi*r_x(i)/rho(j) + phi(j))
      dn_dx(i-1) = dn_dx(i-1) + density_x(i)*r_corr(i)*alpha(j)*2.*pi/rho(j)*cos(2.*pi*r_x(i)/rho(j) + phi(j))
      dn_dx(i-1) = dn_dx(i-1) + density_x(i)*dn_dx(i-1)*beta/density_1AU*r_corr2(i)*alpha(j)*sin(2.*pi*r_x(i)/rho(j) + phi(j))
!       analytic dn/dx rather than straight numerical derivative

!       total = total + alpha(j)*sin(2.*pi*r_x(i)/rho(j) + phi(j)) 	
    ENDDO
	total = 0.

  ENDDO

  IF (rank == n_p-1) dn_dx(n_xp+3) = dn_dx(n_xp+2)
!     absolute last cell in job is given a value...

  
  density_x = density2_x

  CALL MPI_SENDRECV(dn_dx(n_xp+1:n_xp+2),2,MPI_DOUBLE_PRECISION, right, tag, &
  dn_dx(1:2),2,MPI_DOUBLE_PRECISION, left, tag,MPI_COMM_WORLD, status, ecode)

  CALL MPI_SENDRECV(dn_dx(3),1,MPI_DOUBLE_PRECISION, left, tag, &
  dn_dx(n_xp+3),1,MPI_DOUBLE_PRECISION, right, tag,MPI_COMM_WORLD, status, ecode)

  Omega_x = sqrt (4.0*Pi*e*e*density_x/m_e)
 !defines background plasma frequency

  DEALLOCATE(density2_x, r_corr, r_corr2)

  CALL MPI_SENDRECV(density_x(n_xp+1:n_xp+2),2,MPI_DOUBLE_PRECISION, right, tag, &
  density_x(1:2),2,MPI_DOUBLE_PRECISION, left, tag,MPI_COMM_WORLD, status, ecode)

  CALL MPI_SENDRECV(density_x(3),1,MPI_DOUBLE_PRECISION, left, tag, &
  density_x(n_xp+3),1,MPI_DOUBLE_PRECISION, right, tag,MPI_COMM_WORLD, status, ecode)

   l  = 0.5*sqrt(4*pi*e*e/m_e)/sqrt(density_x)*dn_dx/omega_x
 ! Local Plasma Inhomogeneity

  CALL MPI_SENDRECV(l(n_xp+1:n_xp+2),2,MPI_DOUBLE_PRECISION, right, tag, &
  l(1:2),2,MPI_DOUBLE_PRECISION, left, tag,MPI_COMM_WORLD, status, ecode)

  CALL MPI_SENDRECV(l(3),1,MPI_DOUBLE_PRECISION, left, tag, &
  l(n_xp+3),1,MPI_DOUBLE_PRECISION, right, tag,MPI_COMM_WORLD, status, ecode)


END SUBROUTINE density_define

SUBROUTINE corona(dens,sw,x)

USE CONSTANTS

IMPLICIT NONE

  REAL(KIND=8), PARAMETER:: V_1AU = 6.25e+7
    ! wind velocity at the distance 1AU =625km/s
  REAL(KIND=8), PARAMETER:: N_1AU = 6.59
  ! electron number density at 1AU =6.59cm^(-3)
  ! the constant defined experimentally
  ! Mann, G et al A&A, 348, 614-620 (1999)
  REAL(KIND=8), PARAMETER:: C=6.3e+34
  ! see ref. : Mann, G et al A&A, 348, 614-620 (1999)
  REAL(KIND=8), INTENT(IN)::x
  REAL(KIND=8), INTENT(OUT)::dens, sw
  REAL(KIND=8):: v_j, v_k, v_0, v_old, x_c, v_crit
 
  v_crit = sqrt(1d6*k_b/(mu*m_i))
! critical velocity

  x_c = G*m_s/(2.0*v_crit*v_crit) 
! critical distance
  
  IF (x > x_c) then
    v_j = 5e+9
    v_k = v_crit

  ELSE
    v_k = 1e-12
    v_j = v_crit

  ENDIF
! velocity of the solar wind at the distance of 1AU 

  v_0=1e+10
  v_old = 100.

  DO WHILE ((abs(v_0-v_old)/v_0 > abs(1e-12)).and. (abs(PARKER(x,v_0))> 1e-12) )

    v_old=v_0
    v_0 = v_j - PARKER(x,v_j)*(v_j-v_k)/(PARKER(x,v_j)-PARKER(x,v_k))


    IF (Parker(x,v_0)*Parker(x,v_k)< 0.) THEN
      v_j=v_0
    ELSE
      v_k=v_0
    ENDIF

  ENDDO

  sw =v_0
  !final solar wind velocity

  Dens= C/(x*x*v_0) 
  ! final density

  CONTAINS

  DOUBLE PRECISION FUNCTION parker(R, v_wind)
    ! Function of a density and solar wind velocity based
    ! on Parker's model of corona
    USE CONSTANTS
 
    REAL(KIND=8), INTENT(IN)::R, v_wind;
  
    Parker = (v_wind/v_crit)**2 - log(v_wind*v_wind/(v_crit*v_crit))-4.0*log(R/x_c)- 4.0*x_c/R + 3.0

  END FUNCTION parker

END SUBROUTINE corona

END MODULE INITIALISE


