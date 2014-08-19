
MODULE CALC_INDICES

! routines to calculate coefficients and indices for wave-wave processes

IMPLICIT NONE

CONTAINS

SUBROUTINE coeff(Omega_x,landau_arr,gamma_cw,gamma_cp,density_x,log_v_v_t, a2, b2)

  USE CONSTANTS
  USE PARAMETERS

IMPLICIT NONE

REAL(KIND=8), DIMENSION(:),INTENT(IN) :: density_x,Omega_x
REAL(KIND=8), DIMENSION(:):: gamma_cw, gamma_cp, a2, b2
REAL(KIND=8), DIMENSION(-n_v:n_v),INTENT(OUT):: log_v_v_t
REAL(KIND=8), DIMENSION(-n_v:,:):: landau_arr
REAL(KIND=8), DIMENSION(n_xp+3) :: Tql,debye_l,c_log
INTEGER(KIND=4)::i

  Tql = density_x / (omega_x*n_beam)
  ! The quasilinear time

  k1 = n_beam/v_beam;
  ! electron distribution dimensionl parameter

  k2 = m_e*n_beam*v_beam**3/Omega_0;
  ! spectral energy density dimensional parameter

  a1 = v_beam**3*4.*PI*PI*n_beam*e*e/(Omega_0*m_e) 
  ! First equation Coefficient

  a2 = Pi/(Tql*v_beam);
  ! Second equation Coefficient

  b2 = omega_x*e*e*omega_0/(m_e*v_beam**4)
  ! Second spontaneous equation Coefficient

  t_0 = 4.*tau
  ! Characteristic Time

  debye_l = v_t/omega_x
  ! Debye Length of the background plasma

  C_Log = log(debye_l*k_b*t_e/(e*e))
  !  Coulomb Logarith for the Background plasma

  log_v_v_t = log(abs(velocity)/v_t)
  log_v_v_t(0) = 1.

  DO i = -n_v,n_v
	landau_arr(i,:) = sqrt(PI/2.)*omega_x*exp(-(velocity(i)**2)/(2.*v_t**2))*(abs(velocity(i)/v_t)**3)
  ENDDO
	
  landau_arr(0,:) = 0.

  dt_ql_const = tql(20)/16./pi/(n_v-1)**2

  dt = 2d-3
!   fixed maximum main timestep

  gamma_cw = PI*e**4*density_x*(15.89+log(T_e)-0.5*log(density_x))/(m_e*m_e*V_T*V_T*V_T)

  gamma_cp =4.*PI*e**4*density_x*(15.89+log(T_e)-0.5*log(density_x))/(m_e*m_e)



END SUBROUTINE COEFF

SUBROUTINE em_prop_dt

  USE CONSTANTS
  USE PARAMETERS

  IMPLICIT NONE

  REAL(KIND=8):: dt_x, dt_k
  REAL(KIND=8), DIMENSION(-n_kv:n_kv):: omega_emH_om, v_grp, coeff2

 

  omega_emH_om = sqrt(1.+v_c**2/v_t**2 *k_h**2)
  
  v_grp = v_c*v_c*k_h/v_t/omega_emH_om
  coeff2 = v_grp*k_h
  
  dt_x = min_dx/maxval(abs(v_grp))
  dt_k = dk_h/maxval(1./omega_emH_om + coeff2)/max_l

!   print*, dt_x, dt_k
  dt = min(dt_x, dt_k)

!   EM propagation timescale. need to do comms on this. no change in time. based on harmonic as v_grp faster...


END SUBROUTINE em_prop_dt


SUBROUTINE coeff_sound
! define k indices for wave-wave coupling processes

  USE CONSTANTS
  USE PARAMETERS

  IMPLICIT NONE


  REAL(KIND=8) :: k_new, k_max_max, k_max_s, k_prime
  INTEGER(KIND=4)::i

  k_s_stretch =1.
!   stretch the k_s grid by factor, usually either 1 or 2...

  alpha_is = Pi*(4.*pi*e**2/m_e)**2*v_s*k2/(6.*k_b*T_e*v_t**2)
  ! coefficient for decay minus the density change

  gamma_s = sqrt(Pi/2.)*v_s/v_t*(V_s/V_t+(v_s/v_Ti)**3*exp(-v_s*v_s/(2.*v_Ti*v_Ti)))
  ! coefficient for S wave damping minus the k and omega change

  k_prime = v_s/(3.*v_t)
  k_max_max = k_x(4)
  k_max_s = k_s_stretch*k_x(4)


  ind_2pl = 0
  ind_pls = 0

  DO i=-n_kv,-1

    k_new = -k_x(i) - 2.*k_prime

    IF (k_new .le. k_max_max .AND. k_new .ge. -k_max_max)	ind_2pl(i) = minloc(abs(k_x - k_new),1)-n_kv-1
    IF (abs(k_x(ind_2pl(i))) .GE. abs(k_new) ) ind_2pl(i) = ind_2pl(i) + 1
    IF (abs(k_new) .lt. k_prime) ind_2pl(i) = 0
 
!   find k_x nearest to AND LESS THAN k_new. zero index if k_x less than k_prime (no decay occurs)

    k_new = 2.*k_x(i) + 2.*k_prime

    IF (k_new .lt. k_max_s .AND. k_new .gt. -k_max_s)	ind_pls(i) = minloc(abs(k_s_stretch*k_x - k_new),1)-n_kv-1
    IF (abs(k_s_stretch*k_x(ind_pls(i))) .GE. abs(k_new) ) ind_pls(i) = ind_pls(i) - 1
    IF(abs(ind_pls(i)) .GE. n_kv-1) ind_pls(i) = 0


  ENDDO

  DO i=1,n_kv

    k_new = -k_x(i) + 2.*k_prime

    IF (k_new .le. k_max_max .AND. k_new .ge. -k_max_max)	ind_2pl(i) = minloc(abs(k_x - k_new),1)-n_kv-1
    IF (abs(k_x(ind_2pl(i))) .GE. abs(k_new) ) ind_2pl(i) = ind_2pl(i) - 1
    IF (abs(k_new) .lt. k_prime) ind_2pl(i) = 0
        
!   find k_x nearest to AND LESS THAN k_new. zero index if k_x less than k_prime (no decay occurs)

    
    k_new = 2.*k_x(i) - 2.*k_prime

! yes this is a hard coded 2. see momentum-frequency conservation relations

    IF (k_new .lt. k_max_s .AND. k_new .gt. -k_max_s)	ind_pls(i) = minloc(abs(k_s_stretch*k_x - k_new),1)-n_kv-1
    IF (k_s_stretch*k_x(ind_pls(i)) .GE. k_new ) ind_pls(i) = ind_pls(i) + 1
   !     IF (abs(k_new) .le. k_prime) ind_pls(i)=0
    IF(abs(ind_pls(i)) .GE. n_kv-1) ind_pls(i) = 0


  ENDDO


  DO i=-n_kv, -1
    up_2pl(i) = ind_2pl(i+1) - ind_2pl(i) 
!         up_2pl(i) = ind_2pl(i-1) - ind_2pl(i)

!    UP TO NEXT INDEX NB next may be same number...
    frac_2pl(i) = 1./(abs(float(up_2pl(i))))
!   fraction from each resulting point
    frac_2pl(i) = min(frac_2pl(i), 1.)
!   fix in case up is zero
    if(abs(up_2pl(i)) .GE. 1) up_2pl(i) = up_2pl(i) + 1
!   ditto
    if(abs(up_2pl(i)) .GE. 20) up_2pl(i) = 0
!   here indexing is wonky
!
    up_pls(i) = ind_pls(i+1) - ind_pls(i) 
    frac_pls(i) = 1./(abs(float(up_pls(i))))
    frac_pls(i) = min(frac_pls(i), 1.)
    if(abs(up_pls(i)) .GE. 1) up_pls(i) = up_pls(i) - 1
    if(abs(up_pls(i)) .GE. 20) up_pls(i) = 0

!      up_pls(i) = 1./float(abs(ind_pls(i) - ind_pls(i+1)))
!     up_pls(i) = min(up_pls(i), 1)
  
  ENDDO

  DO i=1, n_kv
!     up_2pl(i) = ind_2pl(i+1) - ind_2pl(i) 
    up_2pl(i) = ind_2pl(i-1) - ind_2pl(i)
    frac_2pl(i) = 1./float(abs(up_2pl(i)))
    frac_2pl(i) = min(frac_2pl(i), 1.)
    if(abs(up_2pl(i)) .GE. 1) up_2pl(i) = up_2pl(i) - 1
    if(abs(up_2pl(i)) .GE. 20) up_2pl(i) = 0

   up_pls(i) = ind_pls(i-1) - ind_pls(i)
    frac_pls(i) = 1./float(abs(up_pls(i)))
    frac_pls(i) = min(frac_pls(i), 1.)
    if(abs(up_pls(i)) .GE. 1) up_pls(i) = up_pls(i) + 1
    if(abs(up_pls(i)) .GE. 20) up_pls(i) = 0


  
  ENDDO

  lpc=minloc(abs(k_x - k_prime),1)-n_kv-1

  up_2pl(-1:1) = 0.
  up_pls(-1:1) = 0.

  frac_2pl(-1:1)=1.
  frac_pls(-1:1)=1.



END SUBROUTINE coeff_sound

SUBROUTINE coeff_emf
! define k indices for wave-wave coupling processes

  USE CONSTANTS
  USE PARAMETERS

  IMPLICIT NONE


  REAL(KIND=8) :: k_new, k_max_max, k_max_s,k_max_f
  INTEGER(KIND=4)::i

  kbTe_1d = 4.*pi*k_b*T_e/(2.*pi)**3
  gamma_d_em = 4.*pi*(e**2/m_e)*20./v_t**3*sqrt(2./pi)/3.

! defined these here and in harmonic routine, to make sure gets done...

 alpha_em_f = pi*(1.+3.*(T_i/T_e))/(4.*k_b*T_e)*v_s/(3.*v_t**2) ! *omega**2*omega**2/density*

  ind_fs = 0
  ind_fe = 0


  k_max_max = k_x(4)
  k_max_s = k_s_stretch*k_x(4)
  k_max_f = maxval(k_f)

DO i=-n_kv, -1

  k_new = k_s_stretch*k_x(i)
!   sound wave index for same k_x as Langmuir wave, using stretch factor defined in coeff_sound

  IF (k_new .le. k_max_s .AND. k_new .ge. -k_max_s)	ind_fs(i) = minloc(abs(k_x - k_new),1)-n_kv-1

  k_new = -sqrt(max((3.*k_x(i)**2 - 2.*abs(k_x(i))*v_s/v_t), 0.d0))*v_t/sqrt(v_c**2-3.*v_t**2) 
  
  IF (k_new .le. k_max_f .AND. k_new .ge. -k_max_f)	ind_fe(i) = minloc(abs(k_f - k_new),1)-n_kv-1

  IF (abs(k_f(ind_fe(i))) .GE. abs(k_new) ) ind_fe(i) = ind_fe(i) + 1
  IF (abs(k_new) .LE. minval(abs(k_f))) ind_fe(i)=0
  up_fe(i)=0
  frac_fe(i)=1.


  k_new = -sqrt(max((3.*k_x(i)**2 + 2.*abs(k_x(i))*v_s/v_t), 0.d0))*v_t/sqrt(v_c**2-3.*v_t**2) 
  
  IF (k_new .le. k_max_f .AND. k_new .ge. -k_max_f)	ind_fe2(i) = minloc(abs(k_f - k_new),1)-n_kv-1
  IF (abs(k_f(ind_fe2(i))) .GE. abs(k_new) ) ind_fe2(i) = ind_fe2(i) + 1
  IF (abs(k_new) .LE. minval(abs(k_f))) ind_fe2(i)=0
  up_fe2(i)=0
  frac_fe2(i)=1.


  
ENDDO


DO i=1, n_kv

  k_new = k_s_stretch*k_x(i)
!   sound wave index for same k_x as Langmuir wave, using stretch factor defined in coeff_sound

  IF (k_new .le. k_max_s .AND. k_new .ge. -k_max_s)	ind_fs(i) = minloc(abs(k_x - k_new),1)-n_kv-1

  k_new = sqrt(max((3.*k_x(i)**2 - 2.*k_x(i)*v_s/v_t), 0.d0))*v_t/sqrt(v_c**2-3.*v_t**2)
  
  IF (k_new .le. k_max_f .AND. k_new .ge. -k_max_f)	ind_fe(i) = minloc(abs(k_f - k_new),1)-n_kv-1
  IF (abs(k_f(ind_fe(i))) .GE. abs(k_new) ) ind_fe(i) = ind_fe(i) - 1
  IF (abs(k_new) .LE. minval(abs(k_f))) ind_fe(i)=0

  up_fe(i)=0
  frac_fe(i)=1.


  k_new = sqrt(max((3.*k_x(i)**2 + 2.*k_x(i)*v_s/v_t), 0.d0))*v_t/sqrt(v_c**2-3.*v_t**2)
  
  IF (k_new .le. k_max_f .AND. k_new .ge. -k_max_f)	ind_fe2(i) = minloc(abs(k_f - k_new),1)-n_kv-1
  IF (abs(k_f(ind_fe2(i))) .GE. abs(k_new) ) ind_fe2(i) = ind_fe2(i) - 1
  IF (abs(k_new) .LE. minval(abs(k_f))) ind_fe2(i)=0

  up_fe2(i)=0
  frac_fe2(i)=1.
  
ENDDO


  DO i=-n_kv, -1
    up_fe(i) = ind_fe(i+1) - ind_fe(i) 
!         up_fe(i) = ind_fe(i-1) - ind_fe(i)

!    UP TO NEXT INDEX NB next may be same number...
    frac_fe(i) = 1./(abs(float(up_fe(i))))
!   fraction from each resulting point
    frac_fe(i) = min(frac_fe(i), 1.)
!   fix in case up is zero
    if(abs(up_fe(i)) .GE. 1) up_fe(i) = up_fe(i) + 1
  
    up_fe2(i) = ind_fe2(i+1) - ind_fe2(i) 
!         up_fe(i) = ind_fe(i-1) - ind_fe(i)

!    UP TO NEXT INDEX NB next may be same number...
    frac_fe2(i) = 1./(abs(float(up_fe2(i))))
!   fraction from each resulting point
    frac_fe2(i) = min(frac_fe2(i), 1.)
!   fix in case up is zero
    if(abs(up_fe2(i)) .GE. 1) up_fe2(i) = up_fe2(i) + 1


  ENDDO

  DO i=1, n_kv
    up_fe(i) = ind_fe(i-1) - ind_fe(i)
    frac_fe(i) = 1./float(abs(up_fe(i)))
    frac_fe(i) = min(frac_fe(i), 1.)
    if(abs(up_fe(i)) .GE. 1) up_fe(i) = up_fe(i) - 1

    up_fe2(i) = ind_fe2(i-1) - ind_fe2(i)
    frac_fe2(i) = 1./float(abs(up_fe2(i)))
    frac_fe2(i) = min(frac_fe2(i), 1.)
    if(abs(up_fe2(i)) .GE. 1) up_fe2(i) = up_fe2(i) - 1
  
  ENDDO

  up_fe(-1:1) = 0.
  frac_fe(-1:1)=1.

   up_fe2(-1:1) = 0.
  frac_fe2(-1:1)=1.

! DO i=-n_kv, n_kv
! 
!   print*, i, ind_fs(i), ind_fe(i), ind_fe2(i) 
!   print*, k_x(i), k_x(ind_fs(i)), k_f(ind_fe(i)), k_f(ind_fe2(i))
! ENDDO
! stop

! minimum i to consider for l,s waves
lpe=5


END SUBROUTINE coeff_emf

SUBROUTINE coeff_emh
! define k indices for harmonic emission

  USE CONSTANTS
  USE PARAMETERS

  IMPLICIT NONE


  REAL(KIND=8):: k_temp, omega_emH, cosine, k_1, k_2
  INTEGER(KIND=4)::i

  kbTe_1d = 4.*pi*k_b*T_e/(2.*pi)**3
  gamma_d_em = 4.*pi*(e**2/m_e)*20./v_t**3*sqrt(2./pi)/3.

!   coefficient for harmonic emission, assuming approx angle pi/4
  cosine=sqrt(2.)/2.
  alpha_em_h = pi/(2.*16.*m_e)*(1.-cosine**2) 


  DO i=-n_kv, n_kv
      omega_emH = sqrt(1.+v_c**2/v_t**2 *k_h(i)**2)

      k_temp = (4.*(omega_emH-2.)/3. + (cosine**2-2.)*k_h(i)**2)
      k_1=0.

      IF(k_temp .ge. 0) k_1 =  0.5*k_h(i)*cosine + 0.5*sqrt(k_temp)
    
      k_2=0.
      k_temp= k_h(i)**2 + k_1**2 - 2.*k_h(i)*k_1*cosine

      IF( k_temp .ge.0) k_2 = sqrt(k_temp)
  
      ind_h1(i) = minloc(abs(k_x - k_1), 1)-n_kv-1
      ind_h2(i) = minloc(abs(k_x - k_2), 1)-n_kv-1

   ENDDO


END SUBROUTINE coeff_emh




END MODULE CALC_INDICES