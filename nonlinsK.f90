MODULE NONLIN
! routines for 3-wave interactions and S-wave evolution

CONTAINS

 
SUBROUTINE soundwaves(w_s, w_l_new, omega, density)

  USE CONSTANTS

  USE PARAMETERS

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(-n_kv:n_kv,n_xp+3), INTENT(INOUT) :: W_s, w_l_new
  REAL(KIND=8), DIMENSION(:), INTENT(IN) :: omega, density
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE ::  gamma_s_kd
  REAL(KIND=8), DIMENSION(-n_kv:n_kv) :: omega_s_om, s_terms, damping,s_terms2, a_l2, l_terms
  REAL(KIND=8), DIMENSION(-n_kv:n_kv) :: n_l, n_s, a_l, a_s, a_s2, omega_l_om, a_lnew, a_l2new
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: alpha_is_d
  REAL(KIND=8):: t_is, dt_is, t_is_new, frac1, frac2, NL1, NL2, ns1, ns2
  INTEGER(KIND=4)::i, j
  
  

  ALLOCATE(gamma_s_kd(-n_kv:n_kv,n_xp+3))
  ALLOCATE(alpha_is_d(n_xp+3))


! defining the various constants involved...................................................

  
   
  DO j= 1, n_xp+3
     gamma_s_kd(:,j) = gamma_s*omega(j)*abs(k_x)

  ENDDO

  alpha_is_d = alpha_is*density/omega/omega



! define frequency for L, S waves	..........................................
    omega_l_om = (1.+3.*k_x**2/2.)
    omega_s_om = v_s*abs(k_x)*k_s_stretch/v_t


! subcyling all IS wave processes.......................................................

  DO j= 3, n_xp+2
    t_is=0.
    dt_is = dt_ql/4.

    DO WHILE (abs(t_is - dt_ql) .GE. tiny(1.e0))

! wave occupation numbers.................................................................

      n_l(:) = w_l_new(:,j)/omega_l_om(:)
      n_s(:) = w_s(:,j)/omega_s_om(:)

      n_s(0)=0.
      n_l(0)=0.
!wave wave square bracket terms............................................................
      a_s(:)=0.
      a_s2(:)=0.
      a_l(:)=0.
      a_l2(:)=0.



    DO i=-lpc, -1


 	nl2 = sum(N_L(ind_2pl(i)+up_2pl(i):ind_2pl(i)))
	ns2 = sum(N_s(ind_pls(i):ind_pls(i)+up_pls(i)))


! 	evolution of all wave spectral energy densities in SINGLE L <--> L /pm s encounter
! 	forces conservation of energy in the scattering process
! 	handles both 1-1, 1-many and many-1 situations

	a_l(i) = a_l(i)  + ns2*nl2    !+ Ns1*NL1
	a_l2(i) = a_l2(i)  - nl2 - ns2  !+ nl1 - nl2

	a_l(ind_2pl(i)+up_2pl(i):ind_2pl(i)) = a_l(ind_2pl(i)+up_2pl(i):ind_2pl(i)) + ns2*n_l(i)*frac_2pl(i)
	a_l2(ind_2pl(i)+up_2pl(i):ind_2pl(i)) = a_l2(ind_2pl(i)+up_2pl(i):ind_2pl(i)) - (ns2 -n_l(i))*frac_2pl(i)

	a_s(ind_pls(i):ind_pls(i)+up_pls(i)) = a_s(ind_pls(i):ind_pls(i)+up_pls(i)) + n_l(i)*nl2*frac_pls(i)
	a_s2(ind_pls(i):ind_pls(i)+up_pls(i)) = a_s2(ind_pls(i):ind_pls(i)+up_pls(i)) + (n_l(i) - nl2)*frac_pls(i)


      ENDDO


      DO i=-1, lpc

! wave ocupation numbers involved	

 	nl2 = sum(N_L(ind_2pl(i):ind_2pl(i)+up_2pl(i)))
	ns2 = sum(N_s(ind_pls(i):ind_pls(i)+up_pls(i)))


! 	evolution of all wave spectral energy densities in SINGLE L <--> L /pm s encounter
! 	forces conservation of energy in the scattering process
! 	handles both 1-1, 1-many and many-1 situations

	a_l(i) = a_l(i)  + ns2*nl2    !+ Ns1*NL1
	a_l2(i) = a_l2(i)  - nl2 - ns2  !+ nl1 - nl2

	a_l(ind_2pl(i):ind_2pl(i)+up_2pl(i)) = a_l(ind_2pl(i):ind_2pl(i)+up_2pl(i)) + ns2*n_l(i)*frac_2pl(i)
	a_l2(ind_2pl(i):ind_2pl(i)+up_2pl(i)) = a_l2(ind_2pl(i):ind_2pl(i)+up_2pl(i)) - (ns2 -n_l(i))*frac_2pl(i)

	a_s(ind_pls(i):ind_pls(i)+up_pls(i)) = a_s(ind_pls(i):ind_pls(i)+up_pls(i)) + n_l(i)*nl2*frac_pls(i)
	a_s2(ind_pls(i):ind_pls(i)+up_pls(i)) = a_s2(ind_pls(i):ind_pls(i)+up_pls(i)) + (n_l(i) - nl2)*frac_pls(i)




      ENDDO


      a_L(0)=0.
      A_L2(0)=0.
      a_s(0)=0.
      a_s2(0)=0.

      s_terms= (alpha_is_d(j)*a_s2/2.  -gamma_s_kd(:,j))
      s_terms2 = a_s *alpha_is_d(j)*omega_s_om/2.
      l_terms = a_L2*alpha_is_d(j)
      
!   subcycle timestep determination at given space pt
      dt_is =min(0.1/(maxval(s_terms)+tiny(1.d0)), 1./(maxval(l_terms(-lpc:lpc))+tiny(1.d0)), dt_ql/4.)
!        if(j== 50) print*, j, dt_is
!        dt_is=dt_ql/5.
      t_is_new=min(t_is + dt_is, dt_ql)
      dt_is = t_is_new-t_is
      t_is = t_is_new

!       w_l_new(:,j) = W_L_new(:,j) + alpha_is_d(j)*omega_l_om*dt_is*(a_L + a_l2*N_L)
      w_l_new(:,j) = (W_L_new(:,j) + a_L*alpha_is_d(j)*omega_l_om*dt_is)/(1. - dt_is*l_terms)

      W_s(:,j) = (W_s(:,j) + s_terms2*dt_is)/(1.- dt_is*s_terms)
!       w_s(:,j) = w_s(:,j) + dt_is *(s_terms

!   semi-implicit evolution of wave spectral enery densities

      where(w_s(:,j) .LE. 0) w_s(:,j) = 0.
      w_s(0 , j) = 0.
   

    ENDDO
!     time loop
!   print*, 'done subloop.....................................', j
  ENDDO
!   spatial coordinate loop

  DEALLOCATE( gamma_s_kd, alpha_is_d)

END SUBROUTINE soundwaves


SUBROUTINE soundwaves_with_em(w_s, w_l_new, w_em, omega, density)

  USE CONSTANTS

  USE PARAMETERS

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(-n_kv:n_kv,n_xp+3), INTENT(INOUT) :: W_s, w_l_new, w_em
  REAL(KIND=8), DIMENSION(:), INTENT(IN) :: omega, density
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE ::  gamma_s_kd, gamma_d_f
  REAL(KIND=8), DIMENSION(-n_kv:n_kv) :: omega_s_om, s_terms, damping,s_terms2, a_l2, l_terms, em_terms
  REAL(KIND=8), DIMENSION(-n_kv:n_kv):: omega_emf_om, w_f0_om, w_f0
  REAL(KIND=8), DIMENSION(-n_kv:n_kv) :: n_l, n_s, a_l, a_s, a_s2, omega_l_om, a_eme,a_eml, a_ems, n_em, a_eme2
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: alpha_is_d, alpha_em
  INTEGER(KIND=4)::i, j
  REAL(KIND=8):: t_is, dt_is, t_is_new, frac1, frac2, NL1, NL2, ns1, ns2, nee, nel, nes
  

  ALLOCATE(gamma_s_kd(-n_kv:n_kv,n_xp+3))
  ALLOCATE(gamma_d_f(-n_kv:n_kv,n_xp+3))
  ALLOCATE(alpha_is_d(n_xp+3))
  ALLOCATE(alpha_em(n_xp+3))


! defining the various constants involved...................................................
 
   
  DO j= 1, n_xp+3
     gamma_s_kd(:,j) = gamma_s*omega(j)*abs(k_x)
     gamma_d_f(:,j) = gamma_d_emf*omega(j)*omega(j)

  ENDDO

  alpha_is_d = alpha_is*density/omega/omega
  alpha_em = alpha_em_f*omega*omega/density


  W_f0_om = kbTe_1d *k_f*k_f/(v_t*v_t)


! define frequency for L, S, T waves	..........................................
    omega_l_om = (1.+3.*k_x**2/2.)
    omega_s_om = v_s*abs(k_x)*k_s_stretch/v_t
    omega_emf_om = sqrt(1.+v_c**2/v_t**2 *k_f**2)


! subcyling all IS wave processes.......................................................

  DO j= 3, n_xp+2

    w_f0=w_f0_om*omega(j)*omega(j)

    t_is=0.
    dt_is = dt_ql/4.

    DO WHILE (abs(t_is - dt_ql) .GE. tiny(1.e0))

! wave occupation numbers.................................................................

      n_l = w_l_new(:,j)/omega_l_om
      n_s = w_s(:,j)/omega_s_om
      n_em = w_em(:,j)/omega_emf_om
    
      n_s(0)=0.
      n_l(0)=0.
      n_em(0)=0.
!wave wave square bracket terms............................................................
      a_s=0.
      a_s2=0.
      a_l=0.
      a_l2=0.
      a_eme=0.
      a_eme2=0.
      a_eml=0.
      a_ems=0.

    DO i=-lpc, -1

 	nl2 = sum(N_L(ind_2pl(i)+up_2pl(i):ind_2pl(i)))
	ns2 = sum(N_s(ind_pls(i):ind_pls(i)+up_pls(i)))

! 	evolution of all wave spectral energy densities in SINGLE L <--> L /pm s encounter
! 	forces conservation of energy in the scattering process
! 	handles both 1-1, 1-many and many-1 situations

	a_l(i) = a_l(i)  + ns2*nl2    !+ Ns1*NL1
	a_l2(i) = a_l2(i)  - nl2 - ns2  !+ nl1 - nl2

	a_l(ind_2pl(i)+up_2pl(i):ind_2pl(i)) = a_l(ind_2pl(i)+up_2pl(i):ind_2pl(i)) + ns2*n_l(i)*frac_2pl(i)
	a_l2(ind_2pl(i)+up_2pl(i):ind_2pl(i)) = a_l2(ind_2pl(i)+up_2pl(i):ind_2pl(i)) - (ns2 -n_l(i))*frac_2pl(i)

	a_s(ind_pls(i):ind_pls(i)+up_pls(i)) = a_s(ind_pls(i):ind_pls(i)+up_pls(i)) + n_l(i)*nl2*frac_pls(i)
	a_s2(ind_pls(i):ind_pls(i)+up_pls(i)) = a_s2(ind_pls(i):ind_pls(i)+up_pls(i)) + (n_l(i) - nl2)*frac_pls(i)


      ENDDO

      DO i=1, lpc

! wave ocupation numbers involved	

 	nl2 = sum(N_L(ind_2pl(i):ind_2pl(i)+up_2pl(i)))
	ns2 = sum(N_s(ind_pls(i):ind_pls(i)+up_pls(i)))


! 	evolution of all wave spectral energy densities in SINGLE L <--> L /pm s encounter
! 	forces conservation of energy in the scattering process
! 	handles both 1-1, 1-many and many-1 situations

	a_l(i) = a_l(i)  + ns2*nl2 
	a_l2(i) = a_l2(i)  - nl2 - ns2

	a_l(ind_2pl(i):ind_2pl(i)+up_2pl(i)) = a_l(ind_2pl(i):ind_2pl(i)+up_2pl(i)) + ns2*n_l(i)*frac_2pl(i)
	a_l2(ind_2pl(i):ind_2pl(i)+up_2pl(i)) = a_l2(ind_2pl(i):ind_2pl(i)+up_2pl(i)) - (ns2 -n_l(i))*frac_2pl(i)

	a_s(ind_pls(i):ind_pls(i)+up_pls(i)) = a_s(ind_pls(i):ind_pls(i)+up_pls(i)) + n_l(i)*nl2*frac_pls(i)
	a_s2(ind_pls(i):ind_pls(i)+up_pls(i)) = a_s2(ind_pls(i):ind_pls(i)+up_pls(i)) + (n_l(i) - nl2)*frac_pls(i)


      ENDDO


      DO i=-n_kv, -lpe

	nel = n_l(i)
	nes = n_s(ind_fs(i))
	nee = sum(n_em(ind_fe(i):ind_fe(i)+up_fe(i)))

! 	evolution of all wave spectral energy densities in SINGLE L <--> T + s encounter
! 	forces conservation of energy in the scattering process
! 	WILL BE FIXED TO handle both 1-1, 1-many and many-1 situations

	a_eml(i) = a_eml(i) + nes*nel - nes*nee + nel*nee
 	a_ems(ind_fs(i)) = a_ems(ind_fs(i)) + nes*nel - nes*nee + nel*nee
 	a_eme(ind_fe(i):ind_fe(i)+up_fe(i)) = a_eme(ind_fe(i):ind_fe(i)+up_fe(i)) + (-nes*nel)*frac_fe(i)
	a_eme2(ind_fe(i):ind_fe(i)+up_fe(i)) = a_eme2(ind_fe(i):ind_fe(i)+up_fe(i)) + (nes - nel)*frac_fe(i)

      ENDDO

      DO i=lpe, n_kv

	nel = n_l(i)
	nes = n_s(ind_fs(i))
	nee = sum(n_em(ind_fe(i):ind_fe(i)+up_fe(i)))

! 	evolution of all wave spectral energy densities in SINGLE L <--> T + s encounter
! 	forces conservation of energy in the scattering process
! 	WILL BE FIXED TO handle both 1-1, 1-many and many-1 situations

	a_eml(i) = a_eml(i) + nes*nel - nes*nee + nel*nee
 	a_ems(ind_fs(i)) = a_ems(ind_fs(i)) + nes*nel - nes*nee + nel*nee
!  	a_eme(ind_fe(i)) = a_eme(ind_fe(i)) - nes*nel + nes*nee - nel*nee
	a_eme(ind_fe(i):ind_fe(i)+up_fe(i)) = a_eme(ind_fe(i):ind_fe(i)+up_fe(i)) + (-nes*nel)*frac_fe(i)
	a_eme2(ind_fe(i):ind_fe(i)+up_fe(i)) = a_eme2(ind_fe(i):ind_fe(i)+up_fe(i)) + ( nes - nel)*frac_fe(i)

      ENDDO

        DO i=-n_kv, -lpe

	nel = n_l(i)
	nes = n_s(ind_fs(i))
	nee = sum(n_em(ind_fe2(i):ind_fe2(i)+up_fe2(i)))

! 	evolution of all wave spectral energy densities in SINGLE L <--> T - s encounter
! 	forces conservation of energy in the scattering process
! 	WILL BE FIXED TO handle both 1-1, 1-many and many-1 situations

	a_eml(i) = a_eml(i) - nes*nel + nes*nee + nel*nee
 	a_ems(ind_fs(i)) = a_ems(ind_fs(i)) - nes*nel + nes*nee + nel*nee
 	a_eme(ind_fe2(i):ind_fe2(i)+up_fe2(i)) = a_eme(ind_fe2(i):ind_fe2(i)+up_fe2(i)) + (nes*nel)*frac_fe2(i)
	a_eme2(ind_fe2(i):ind_fe2(i)+up_fe2(i)) = a_eme2(ind_fe2(i):ind_fe2(i)+up_fe2(i)) + (- nes - nel)*frac_fe2(i)


      ENDDO

      DO i=lpe, n_kv

	nel = n_l(i)
	nes = n_s(ind_fs(i))
	nee = sum(n_em(ind_fe2(i):ind_fe2(i)+up_fe2(i)))

! 	evolution of all wave spectral energy densities in SINGLE L <--> T - s encounter
! 	forces conservation of energy in the scattering process
! 	WILL BE FIXED TO handle both 1-1, 1-many and many-1 situations

	a_eml(i) = a_eml(i) - nes*nel + nes*nee + nel*nee
 	a_ems(ind_fs(i)) = a_ems(ind_fs(i)) - nes*nel + nes*nee + nel*nee
!  	a_eme(ind_fe(i)) = a_eme(ind_fe(i)) - nes*nel + nes*nee - nel*nee
	a_eme(ind_fe2(i):ind_fe2(i)+up_fe2(i)) = a_eme(ind_fe2(i):ind_fe2(i)+up_fe2(i)) + (nes*nel)*frac_fe2(i)
	a_eme2(ind_fe2(i):ind_fe2(i)+up_fe2(i)) = a_eme2(ind_fe2(i):ind_fe2(i)+up_fe2(i)) + (- nes - nel)*frac_fe2(i)

      ENDDO


      a_L(0)=0.
      A_L2(0)=0.
      a_s(0)=0.
      a_s2(0)=0.
      a_eml(0)=0.
      a_eme(0)=0.
      a_ems(0)=0.


      s_terms= (alpha_is_d(j)*a_s2/2.  -gamma_s_kd(:,j))
      s_terms2 = a_s *alpha_is_d(j)*omega_s_om/2.
      l_terms = a_L2*alpha_is_d(j)
      em_terms = a_eme2*alpha_em(j)*k2
      
!   subcycle timestep determination at given space pt
      dt_is =min(0.1/(maxval(s_terms)+tiny(1.d0)),0.1/(maxval(em_terms)+tiny(1.d0)), &
  1./(maxval(l_terms(-lpc:lpc))+tiny(1.d0)), dt_ql/8.)
!         if(j== n_xp/2) print*,rank, dt_is, 1./(maxval(l_terms(-lpc:lpc))+tiny(1.d0))
!          dt_is=dt_ql/8.
      t_is_new=min(t_is + dt_is, dt_ql)
      dt_is = t_is_new-t_is
      t_is = t_is_new

!       w_l_new(:,j) = W_L_new(:,j) + alpha_is_d(j)*omega_l_om*dt_is*(a_L + a_l2*N_L)
      w_l_new(:,j) = (W_L_new(:,j) + a_L*alpha_is_d(j)*omega_l_om*dt_is)/(1. - dt_is*l_terms)

      W_s(:,j) = (W_s(:,j) + s_terms2*dt_is)/(1.- dt_is*s_terms)
!       w_s(:,j) = w_s(:,j) + dt_is *(s_terms

!   semi-implicit evolution of wave spectral enery densities

     W_em(:,j) = w_em(:,j)  + (a_eme + a_eme2*n_em)*alpha_em(j)*k2*dt_is   &
 	  - (w_em(:,j) - w_f0)*gamma_d_f(:,j)/omega_emf_om**2 *dt_is
      w_l_new(:,j) = w_l_new(:,j) + a_eml*alpha_em(j)*dt_is*k2
      w_s(:,j) = w_s(:,j)+ a_ems*alpha_em(j)*dt_is*k2

      where(w_s(:,j) .LE. 0) w_s(:,j) = 0.
      w_s(0 , j) = 0.
      where (w_em(:,j) .le. w_f0) w_em(:,j) = w_f0

  

    ENDDO
  
!     time loop
!    print* 'done subloop.....................................', j
  ENDDO
!   spatial coordinate loop
! print*, rank, 'done space loop'
  DEALLOCATE( gamma_s_kd, alpha_is_d, alpha_em, gamma_d_f)

END SUBROUTINE soundwaves_with_em


SUBROUTINE harmonic(w_l, w_em_new, omega, density)

  USE CONSTANTS

  USE PARAMETERS

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(-n_kv:, :), INTENT(INOUT) :: W_L, w_em_new
  REAl(KIND=8), DIMENSION(:),INTENT(IN) :: omega, density
  REAL(KIND=8), DIMENSION(-n_kv:n_kv):: a_em, n_l, w_h0, omega_emh_om, n_em, w_em_tmp,omega_emh2,W_h0_om
  INTEGER(KIND=4)::i, j
  REAL(KIND=8):: alpha_em, gamma_d, omega_sq, cosine, k_1, k_2, a_em_tmp


! define space unvarying parts of various quantities.....................................................

  omega_emH_om = sqrt(1.+v_c**2/v_t**2 *k_h**2)
  W_h0_om = kbTe_1d *k_h*k_h/(v_t*v_t)
  
  cosine = sqrt(2.)/2.
  


  DO j = 3, n_xp+2
!     timesaveing variable
    omega_sq=omega(j)*omega(j)

    omega_emH2=omega_emH_om*omega(j)
    alpha_em = omega(j)/density(j)*alpha_em_h  
    gamma_d=omega_sq*gamma_d_em
  !     coefficients

    w_h0=w_h0_om*omega_sq
    ! thermal level

    n_em = w_em_new(:,j)/omega_emH2
  !  Defining n_em but it appears to be too low to make a difference to the result

    n_l = w_l(:,j)/(omega(j) + 3.*omega(j)*k_x**2/2.)
  ! number density of Langmuir waves

    a_em = 0.


    k_1=k_x(ind_h1(-n_kv))
    k_2=k_x(ind_h2(-n_kv))
    
    a_em(-n_kv) = (omega(j)/v_t)**2*(k_1**2 - k_2**2)**2/k_2**2/abs(3.*v_t*(2.*k_1-abs(k_h(-n_kv))*cosine))* &
      (k_h(-n_kv)**2/k_2**2*16.*(n_l(ind_h1(-n_kv))*n_l(-ind_h2(-n_kv))) -n_em(-n_kv)*(n_l(ind_h1(-n_kv))+ n_l(-ind_h2(-n_kv)) ))

    k_1=k_x(ind_h1(-n_kv+1))
    k_2=k_x(ind_h2(-n_kv+1))
    
    a_em_tmp = (omega(j)/v_t)**2*(k_1**2 - k_2**2)**2/k_2**2/abs(3.*v_t*(2.*k_1-abs(k_h(-n_kv+1))*cosine))* &
      (k_h(-n_kv+1)**2/k_2**2*16.*(n_l(ind_h1(-n_kv+1))*n_l(-ind_h2(-n_kv+1))) &
  -n_em(-n_kv+1)*(n_l(ind_h1(-n_kv+1))+ n_l(-ind_h2(-n_kv+1)) ))

    a_em(-n_kv) = a_em(-n_kv) + a_em_tmp
    

    DO i=-n_kv+1, n_kv-1
     
      k_1=k_x(ind_h1(i+1))
      k_2=k_x(ind_h2(i+1))

      a_em(i) = a_em_tmp

      a_em_tmp = (omega(j)/v_t)**2*(k_1**2- &
      k_2**2)**2/k_2**2/abs(3.*v_t*(2.*k_1-abs(k_h(i))*cosine))* &
      (k_h(i)**2/k_2**2*16.*(n_l(ind_h1(i+1))*n_l(-ind_h2(i+1))) -n_em(i+1)*(n_l(ind_h1(i+1))+ n_l(-ind_h2(i+1)) ))
  
      a_em(i)=(a_em(i) + a_em_tmp)*0.5

    ENDDO

    k_1=k_x(ind_h1(n_kv))
    k_2=k_x(ind_h2(n_kv))
    
    a_em(n_kv) = (omega(j)/v_t)**2*(k_1**2 - k_2**2)**2/k_2**2/abs(3.*v_t*(2.*k_1-abs(k_h(n_kv))*cosine))* &
      (k_h(n_kv)**2/k_2**2*16.*(n_l(ind_h1(n_kv))*n_l(-ind_h2(n_kv))) -n_em(n_kv)*(n_l(ind_h1(n_kv))+ n_l(-ind_h2(n_kv)) ))



    a_em(0) = 0.

    w_em_new(:,j) = w_em_new(:,j)  + alpha_em*dt*a_em*k2*omega_emH2   - (w_em_new(:,j)- w_h0)*gamma_d/omega_emH_om**2 *dt 

    where (w_em_new(:,j) .le. w_h0) w_em_new(:,j) = w_h0

  ENDDO


END SUBROUTINE harmonic
 

END MODULE NONLIN
