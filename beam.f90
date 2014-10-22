PROGRAM QBEAM_PLASMA

  ! This program solves the one-dimensional kinetic equations
  ! of quasilinear relaxation using implicit
  ! difference scheme (second order approximation over coordinate
  ! and velocity, and first order approximation over time)

  ! MPI Version created by :
  ! H. Ratcliffe/C.S.Brady, 2014
  ! Based on a previous code by Reid/Kontar, University of Glasgow
 
  USE CONSTANTS
  USE PARAMETERS
  USE I_O
  USE INITIALISE
  USE DISTRIBSINIT
  USE CALC_INDICES
  USE QUASILINEAR_TERMS
  USE ELECTRON_SOURCE
  USE NONLIN
  USE EMPROP
  USE APPLYBOUNDS
  USE MPI

  IMPLICIT NONE 

  INTEGER:: WriteStat, i, j, t1, counter, r_zero
  ! loop and sevice variables
  
  INTEGER, DIMENSION(1):: min_tmp
! temporary array for finding minlocs

  REAL(KIND=8):: t_ql, time_begin, time_prev, t_new
!   time variables
 
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE:: R_x, delta_x
  !whole domain grid and dx
	
  REAL(KIND=8), DIMENSION(-n_v:n_v) :: log_v_v_t = 0.0
  ! Velocity arrays
  REAL(KIND=8), DIMENSION(-n_kv:n_kv) ::  coeff_h, coeff_f

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE:: omega_h_om, omega_f_om

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE:: r_x_p, delta_x_p, omega_p, density_p, density_smooth_p, dn_dx_p
!   per processor space grid, space grid dx, frequency, density, smooth density, dn/dx

  REAL(KIND=8), DIMENSION(:), ALLOCATABLE:: a2_p, l_p, gamma_cp_p, gamma_cw_p, b2_p, wind_p
! per processor coefficients

  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE::F_p, W_p,W_t_p,f_new_p,w_new_p,f_p_ql,w_p_ql, landau_arr_p
!   per processor electron distributions, spectral energy densities, coefficient
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE::w_s_p, w_s_old_p,w_em_p,w_em_f_p

!   mpi service variables
  
  REAL(KIND=8):: max_w, max_w_p, w_max_all, clip
!   maximum Langmuir wave energy density

!------------  MPI initialising......................................................................

CALL MPI_init(ecode)

CALL MPI_comm_rank(MPI_COMM_WORLD, rank, ecode)

CALL MPI_comm_size(MPI_COMM_WORLD, n_p, ecode)

IF (rank==0) THEN
  time_begin = MPI_WTIME()
  time_prev=time_begin
ENDIF
! log starting time of main calculations within MPI

IF (rank .NE. 0) left=rank-1
IF (rank .NE. n_p-1) right=rank+1
! define number of processor to left and right. NB initialised to null (do nothing)

tag=30
! generic MPI tag number

IF (rank == 0) CALL ReadInit
!   reading initial parameters on root

! share initial parameters between processors..............................................................
CALL MPI_BCAST(T_max, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ecode)
CALL MPI_BCAST(d, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ecode)
CALL MPI_BCAST(x_0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ecode)
CALL MPI_BCAST(x_min, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ecode)
CALL MPI_BCAST(time_print, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ecode)
CALL MPI_BCAST(time_save, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ecode)
CALL MPI_BCAST(t_start, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ecode)
CALL MPI_BCAST(restart_t, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ecode)
CALL MPI_BCAST(dyn_spec_t, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ecode)
CALL MPI_BCAST(n_beam, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ecode)
CALL MPI_BCAST(v_beam, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ecode)
CALL MPI_BCAST(tau, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ecode)
CALL MPI_BCAST(tau2, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ecode)
CALL MPI_BCAST(t_e, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ecode)
CALL MPI_BCAST(t_i, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ecode)

CALL MPI_BCAST(n_x, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ecode)
CALL MPI_BCAST(N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ecode)
CALL MPI_BCAST(dumpArr, 6, MPI_INTEGER, 0, MPI_COMM_WORLD, ecode)
CALL MPI_BCAST(dumpSingle, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ecode)
CALL MPI_BCAST(rollRestart, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ecode)

CALL MPI_BCAST(path, 50, MPI_CHARACTER, 0, MPI_COMM_WORLD, ecode)
! ---------------------------------------------------------------------------------------------------------

n_xp = ((n_x-3)/n_p)
! number of x points each processor has
IF(rank == 0) CALL num_points_check
CALL MPI_BCAST(n_x, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ecode)

!define position and delta x over whole domain, subdivide the deallocate long arrays...............................
ALLOCATE (r_x(n_x))
ALLOCATE (delta_x(n_x))

CALL space_grid_define(x_min, r_x, delta_x)

min_dx=minval(abs(delta_x))

ALLOCATE (r_x_p(n_xp+3))
ALLOCATE (delta_x_p(n_xp+3))

r_x_p = r_x(rank*n_xp+1:(rank+1)*n_xp+3)
delta_x_p = delta_x(rank*n_xp+1:(rank+1)*n_xp+3)
  
IF(minval(abs(r_x)) .GE. min_dx/2.) THEN
 print*, "Duude, negative space!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
 stop
 !abort if we don't resolve the injection region...
ENDIF

min_tmp=minloc(abs(r_x_p))
r_zero=0
IF(minval(abs(r_x_p)) .LE. min_dx/2.) r_zero = rank
CALL MPI_ALLREDUCE(MPI_IN_PLACE, r_zero, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ecode)
! location of r_x=~0, including rank of processor

DEALLOCATE(r_x,delta_x)
! deallocate whole domain space gridding
!=--------------------------------------------------------------------------------------------------------

CALL velocity_define
! defining velocity and wavenumber (note global)

CALL em_wavevec_define
!defines EM wavenumbers k_f (fundamental) and k_h (harmonic) from k_x

! allocate main arrays ...................................................................................
ALLOCATE (f_p(-n_v:n_v,n_xp+3))
ALLOCATE (f_new_p(-n_v:n_v,n_xp+3))
ALLOCATE (f_p_ql(-n_v:n_v,n_xp+3))
ALLOCATE (w_p(-n_kv:n_kv,n_xp+3))
ALLOCATE (w_new_p(-n_kv:n_kv,n_xp+3))
ALLOCATE (w_t_p(-n_kv:n_kv,n_xp+3))
ALLOCATE (w_p_ql(-n_kv:n_kv,n_xp+3))
ALLOCATE (w_s_p(-n_kv:n_kv,n_xp+3))
ALLOCATE (w_em_p(-n_kv:n_kv,n_xp+3))
ALLOCATE (w_em_f_p(-n_kv:n_kv,n_xp+3))
! electron and wave distributions

ALLOCATE (omega_p(n_xp+3))
ALLOCATE (density_p(n_xp+3))
ALLOCATE (density_smooth_p(n_xp+3))
ALLOCATE (dn_dx_p(n_xp+3))
! frequency and density

ALLOCATE (wind_p(n_xp+3))
ALLOCATE (a2_p(n_xp+3))
ALLOCATE (b2_p(n_xp+3))
ALLOCATE (l_p(n_xp+3))
ALLOCATE (gamma_cw_p(n_xp+3))
ALLOCATE (gamma_cp_p(n_xp+3))
ALLOCATE (landau_arr_p(-n_v:n_v,n_xp+3))
! coefficients

! -------------------------------------------------------------------------------------------------
  CALL user_mpi_init

  CALL density_define(density_p, density_smooth_p, omega_p, wind_p, R_x_p, dn_dx_p, delta_x_p, l_p)
! define coronal density, smooth and fluctuating, plasma frequency, solar windspeed

  IF (t_start .GT. 0.) THEN
      CALL restart_density(density_p, density_smooth_p, omega_p, wind_p, R_x_p, dn_dx_p, delta_x_p, l_p)
  ENDIF
!   if restarting, read old density profile and setup SAME fluctuaitons
 
  IF(rank == r_zero) omega_0 =  omega_p(min_tmp(1))

  CALL MPI_BCAST(omega_0, 1, MPI_DOUBLE_PRECISION, r_zero, MPI_COMM_WORLD, ecode)
! define reference plasma frequency at injection point and bcast to all processors

  CALL write_density(r_x_p,density_p, density_smooth_p, wind_p, dn_dx_p, omega_p, l_p)
!   write density profile, increment etc to file
  
  CALL coeff(omega_p,landau_arr_p,gamma_cw_p,gamma_cp_p,density_p,log_v_v_t,a2_p,b2_p)
 ! calculates coefficients for quasilinear processes
 
 
 
  CALL coeff_sound
  CALL coeff_emh
  CALL coeff_emf
! calculates coefficients and wavenumber matching conditions for wave-wave processes
! if running only soundwaves (not soundwaves_with_em) can comment out coeff_emf

  CALL write_coeffs

  CALL distrib_init(f_p, w_p, w_s_p, w_t_p, w_em_p, w_em_f_p, omega_p, r_x_p)
  ! initialises the distribution functions

  max_l=maxval(abs(l_p))
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, max_l, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ecode)
  CALL em_prop_dt
  CALL MPI_ALLREDUCE(MPI_IN_PLACE, dt, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ecode)
! these determine main timestep so em propagation works ok... otherwise uses default 2e-3 from coeff routine

  dt_ql = dt
  t_ql = 0.d0
  t = 0.d0
  clip = 0.d0

  IF (time_save .LE. dt) THEN
    time_save=dt+tiny(0.d0)
    IF(rank==0) write(*,'(A, E10.3, A)'), 'Save time too short, changed to', time_save, ' s'
  ENDIF
!     because we should't save on time shorter than main timestep dt
  
  IF (rank == 0) CALL write_velocity_params
! write velocity, wavenumbers, sizes and other parameters to file

  IF (t_start .GT. 0.) THEN
 
    t = t_start
    t_ql = t_start
    print*, FLOOR(t_start/restart_t)
    CALL restart(f_p, w_p, w_s_p, w_em_p, w_em_f_p, FLOOR(t_start/restart_t))

    IF(rank == 0) print*, 'READ OLD DATA, RESTARTING...................................................'

  ENDIF
 
  f_new_p = f_p
  w_new_p = w_p


  IF (dumpSingle .LE. 0) THEN
    CALL write_profiles (f_p, w_p, w_t_p, w_s_p, w_em_p, w_em_f_p, FLOOR(t/time_save), omega_p, delta_x_p)
  ELSE
    CALL write_profiles_single (REAL(f_p, 4), REAL(w_p, 4), REAL(w_t_p, 4), REAL(w_s_p, 4) &
  , REAL(w_em_p,4), REAL(w_em_f_p,4), FLOOR(t/time_save), REAL(omega_p,4), REAL(delta_x_p,4))
 ! writing inital profiles to file
  ENDIF


 IF (rank == 0) CALL write_parameters 
 ! Writing some key variables to the screen
  


! MAIN LOOP starts here ***********************************************

main :    DO WHILE (t .le. t_max)

!   print*, 'main_loop'
  
  IF (time_print - MOD(t, time_print) .LE. dt) THEN 
    IF(rank == 0) THEN 
      CALL write_progress(time_begin, time_prev)
      IF(clip .NE. 0) CALL write_warnings(1)
    ENDIF
    
  ENDIF

  IF (time_save - MOD(t, time_save) .LE. dt) THEN 
  
   IF (dumpSingle .LE. 0) THEN
       CALL write_profiles (f_p, w_p, w_t_p, w_s_p, w_em_p, w_em_f_p, FLOOR(t/time_save), omega_p, delta_x_p)
    ELSE
      CALL write_profiles_single (REAL(f_p, 4), REAL(w_p, 4), REAL(w_t_p, 4), REAL(w_s_p, 4) & 
  ,REAL(w_em_p,4), REAL(w_em_f_p,4), FLOOR(t/time_save), REAL(omega_p,4), REAL(delta_x_p,4))
  ! writing profiles from each processor to a combined data file
  ! select former for double precision, latter for single AND ALSO ABOVE
    ENDIF
  ENDIF

  IF (restart_t - MOD(t, restart_t) .LE. dt) THEN 
    IF(rank == 0) print*, 'Writing RESTARTABLE', FLOOR(t/restart_t)
    CALL write_profiles_restart (f_p, w_p, w_t_p, w_s_p, w_em_p, w_em_f_p, FLOOR(t/restart_t), &
omega_p, delta_x_p)
!     writing full precision, all data for restarting with...

  ENDIF
  IF (dyn_spec_t - MOD(t, dyn_spec_t) .LE. dt) THEN 
    CALL write_dynSpec(w_em_p, w_em_f_p,FLOOR(t/dyn_spec_t),omega_p, delta_x_p)
  ENDIF

  t = t + dt

!  f_p_ql = f_p
!  w_p_ql = w_p
  
  !f_new_p=f_p
  
  dt_ql_av = 0.
  counter = 0.
  clip = 0.d0

! quasilinear subcycle..........................................................................................

DO while (abs(t - t_ql) .GE. tiny(1.e0)) 
  max_w_p = maxval(w_p)
!   print*, max_w_p
  dt_ql = min(dt_ql_const/max_w_p, 2d-4)
  clip=dt_ql
  dt_ql = max(1d-8, dt_ql)
  clip=clip-dt_ql

!   clip timestep if it gets stupidly low.................................................

  CALL MPI_ALLREDUCE(MPI_IN_PLACE, dt_ql, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ecode)

  dt_ql_av = dt_ql_av + dt_ql
  counter = counter + 1

  t_new = min(t_ql + dt_ql, t)
  dt_ql = t_new - t_ql
  t_ql = t_new


!   routines for self-consitent evolution of electrons and Langmuir waves ....................

!  CALL spontaneous(f_p_ql , w_new_p , log_v_v_t , b2_p)
!print*, dt_ql 
!print*, gamma_cp_p(20)
 CALL collision(f_p , f_new_p , w_new_p , gamma_cp_p , gamma_cw_p , w_p , w_t_p)

 CALL inhomogeneity(w_p_ql , w_new_p , l_p)

 CALL quasilinear(f_p , w_p , f_new_p , w_new_p , a2_p)
 
 CALL landau_background(w_p , w_new_p , landau_arr_p)

 CALL electron_bounds(f_new_p)
 CALL langmuir_bounds(w_new_p, w_t_p)


!  --------------------------------------------------------------------------------------------------

!routines for processes involving ion-sound waves ..............................................
! this emcompasses all of the IS wave processes. Because these are the time-step determing process now
!   we want to absolutely minimise the calculation in there. So we subcycle WITHIN this routine as necessary... number of iterations is t_is which is internall determined


 ! CALL soundwaves(w_s_p, w_new_p, omega_p, density_p)
!  includes langmuir wave decay to L +S, S wave damping 

 CALL soundwaves_with_em(w_s_p, w_new_p, w_em_f_p, omega_p, density_p)
! sound wave routines including fundamental emission. use THIS OR soundwaves
CALL langmuir_bounds(w_new_p, w_t_p)

 CALL MPI_BARRIER(MPI_COMM_WORLD, ecode)
! within soundwaves we have differing timestep over space. protect against one processor getting here earlier
! ----------------------------------------------------------------------------------------------------

!routine for Langmuir wave scattering off ions______________________________________________



!CALL ion_scatter(w_p_ql, w_new_p, omega_p, density_p)




!-----------------------------------------------------------------------------------





  where (w_new_p(-n_kv:n_kv,:) .le. w_t_p(-n_kv:n_kv,:)) w_new_p(-n_kv:n_kv,:) = w_t_p(-n_kv:n_kv,:)

  where (f_new_p(-n_v+2:n_v-2,:) .le. 1d-30) f_new_p(-n_v+2:n_v-2,:)=0. !1d-60
  
  !bnd conds may introduce -ves, so make sure they go after these zero protects
  !and keep the edge exclusions 
 
  CALL electron_bounds(f_new_p)
  CALL langmuir_bounds(w_new_p, w_t_p)
  CALL fundamental_bounds(w_em_f_p)
 
  
  w_p=w_new_p	
  f_p=f_new_p

!   where (w_p_ql(-n_kv:n_kv,:) .le. w_t_p(-n_kv:n_kv,:)) w_p_ql(-n_kv:n_kv,:) = w_t_p(-n_kv:n_kv,:)
!   w_new_p=w_p_ql  

! end of quasilinear subcycle...
ENDDO
dt_ql_av = dt_ql_av/counter


  CALL harmonic(w_new_p, w_em_p, omega_p, density_p)

  CALL upwind_harm(w_em_p, omega_p, l_p, delta_x_p)

  CALL upwind_fund(w_em_f_p, omega_p, l_p, delta_x_p)

! propagation of waves and electrons.............................................................................
CALL vanleer(f_p , f_new_p , delta_x_p)

CALL radial(f_p , f_new_p , r_x_p)

  CALL upwind(w_p , w_new_p , delta_x_p)

  IF (t <= 4*tau+4*tau2) CALL fsource(f_new_p, r_x_p)
  
  CALL electron_bounds(f_new_p)
  CALL langmuir_bounds(w_new_p, w_t_p)
  CALL fundamental_bounds(w_em_f_p)
  CALL harmonic_bounds(w_em_p)
  
  
  
  where (w_new_p(-n_kv:n_kv,:) .le. w_t_p(-n_kv:n_kv,:)) w_new_p(-n_kv:n_kv,:) = w_t_p(-n_kv:n_kv,:)


! Exchange end points of space domains between processors
!  
  CALL MPI_SENDRECV(f_new_p(:,n_xp+1:n_xp+2),(2*n_v+1)*2,MPI_DOUBLE_PRECISION, right, tag, &
            f_new_p(:,1:2),(2*n_v+1)*2,MPI_DOUBLE_PRECISION, left, tag,MPI_COMM_WORLD, status, ecode)

  CALL MPI_SENDRECV(f_new_p(:,3),(2*n_v+1),MPI_DOUBLE_PRECISION, left, tag, &
            f_new_p(:,n_xp+3),(2*n_v+1),MPI_DOUBLE_PRECISION, right, tag,MPI_COMM_WORLD, status, ecode)


  CALL MPI_SENDRECV(w_new_p(:,n_xp+1:n_xp+2),(2*n_kv+1)*2,MPI_DOUBLE_PRECISION, right, tag, &
            w_new_p(:,1:2),(2*n_kv+1)*2,MPI_DOUBLE_PRECISION, left, tag,MPI_COMM_WORLD, status, ecode)

  CALL MPI_SENDRECV(w_new_p(:,3),(2*n_kv+1),MPI_DOUBLE_PRECISION, left, tag, &
            w_new_p(:,n_xp+3),(2*n_kv+1),MPI_DOUBLE_PRECISION, right, tag,MPI_COMM_WORLD, status, ecode)


  CALL MPI_SENDRECV(w_em_p(:,n_xp+1:n_xp+2),(2*n_kv+1)*2,MPI_DOUBLE_PRECISION, right, tag, &
            w_em_p(:,1:2),(2*n_kv+1)*2,MPI_DOUBLE_PRECISION, left, tag,MPI_COMM_WORLD, status, ecode)

  CALL MPI_SENDRECV(w_em_p(:,3),(2*n_kv+1),MPI_DOUBLE_PRECISION, left, tag, &
            w_em_p(:,n_xp+3),(2*n_kv+1),MPI_DOUBLE_PRECISION, right, tag,MPI_COMM_WORLD, status, ecode)


  CALL MPI_SENDRECV(w_em_f_p(:,n_xp+1:n_xp+2),(2*n_kv+1)*2,MPI_DOUBLE_PRECISION, right, tag, &
            w_em_f_p(:,1:2),(2*n_kv+1)*2,MPI_DOUBLE_PRECISION, left, tag,MPI_COMM_WORLD, status, ecode)

  CALL MPI_SENDRECV(w_em_f_p(:,3),(2*n_kv+1),MPI_DOUBLE_PRECISION, left, tag, &
            w_em_f_p(:,n_xp+3),(2*n_kv+1),MPI_DOUBLE_PRECISION, right, tag,MPI_COMM_WORLD, status, ecode)

	f_p = f_new_p
	w_p = w_new_p


END DO main
! End of main loop

!  System clock time at end
IF (rank == 0) THEN
  WRITE(*,*) "CALCULATIONS COMPLETED     -    OK !"
  CALL write_total_time(time_begin)
ENDIF

CALL user_mpi_fin

CALL MPI_Finalize(ecode)

DEALLOCATE(r_x_p, delta_x_p, omega_p, density_p,  density_smooth_p, dn_dx_p)
! deallocate per processor grid, frequency, density, smooth density, dn/dx

DEALLOCATE(f_p, f_p_ql, f_new_p, w_p, w_p_ql, w_new_p, w_t_p, w_s_p, w_em_p, w_em_f_p)
! deallocate electron and wave distributions

DEALLOCATE(a2_p, b2_p, landau_arr_p, l_p, Wind_p)
! deallocate coefficients


END PROGRAM
