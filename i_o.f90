MODULE I_O
! Subroutines to write data to files, read input file etc

IMPLICIT NONE

CONTAINS

SUBROUTINE user_mpi_init
!   create the mpi types required for data read write, both double and single precision versions

  USE CONSTANTS
  USE PARAMETERS 
  USE MPI

  IMPLICIT NONE


 sizes=(/2*n_v+1,n_x-3/)
 subsizes=(/2*n_v+1,n_xp/)
 starts=(/0, n_xp*rank/)
 CALL MPI_TYPE_CREATE_SUBARRAY(2,sizes,subsizes,starts,MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mpitype, ecode)
 CALL MPI_TYPE_COMMIT(mpitype, ecode) 

 sizes=(/2*n_kv+1,n_x-3/)
 subsizes=(/2*n_kv+1,n_xp/)
 starts=(/0, n_xp*rank/)
 CALL MPI_TYPE_CREATE_SUBARRAY(2,sizes,subsizes,starts,MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mpitypeL, ecode)
 CALL MPI_TYPE_COMMIT(mpitypeL, ecode)

 sizes1=(/n_x-3/)
 subsizes1=(/n_xp/)
 starts1=(/n_xp*rank/)
 CALL MPI_TYPE_CREATE_SUBARRAY(1,sizes1,subsizes1,starts1,MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mpitype1, ecode)
 CALL MPI_TYPE_COMMIT(mpitype1, ecode)

 sizes=(/2*n_v+1,n_x-3/)
 subsizes=(/2*n_v+1,n_xp/)
 starts=(/0, n_xp*rank/)
 CALL MPI_TYPE_CREATE_SUBARRAY(2,sizes,subsizes,starts,MPI_ORDER_FORTRAN, MPI_FLOAT, mpitypeF, ecode)
 CALL MPI_TYPE_COMMIT(mpitypeF, ecode) 
 
 sizes=(/2*n_kv+1,n_x-3/)
 subsizes=(/2*n_kv+1,n_xp/)
 starts=(/0, n_xp*rank/)
 CALL MPI_TYPE_CREATE_SUBARRAY(2,sizes,subsizes,starts,MPI_ORDER_FORTRAN, MPI_FLOAT, mpitypeLF, ecode)
 CALL MPI_TYPE_COMMIT(mpitypeLF, ecode)
 
 sizes1=(/n_x-3/)
 subsizes1=(/n_xp/)
 starts1=(/n_xp*rank/)
 CALL MPI_TYPE_CREATE_SUBARRAY(1,sizes1,subsizes1,starts1,MPI_ORDER_FORTRAN, MPI_FLOAT, mpitype1f, ecode)
 CALL MPI_TYPE_COMMIT(mpitype1F, ecode)


END SUBROUTINE user_mpi_init


SUBROUTINE user_mpi_fin
!   free the mpi types after completion

  USE PARAMETERS 
  USE MPI

  IMPLICIT NONE


 CALL MPI_TYPE_FREE(mpitype, ecode)
 CALL MPI_TYPE_FREE(mpitypeL, ecode)
 CALL MPI_TYPE_FREE(mpitype1, ecode)
 CALL MPI_TYPE_FREE(mpitypeF, ecode)
 CALL MPI_TYPE_FREE(mpitypeLF, ecode)
 CALL MPI_TYPE_FREE(mpitype1F, ecode)

END SUBROUTINE user_mpi_fin

SUBROUTINE  write_profiles (fxyz, wxyz, w_t, w_s, w_em, w_em_f, t_step, omega_p,delta_x_p)
! The procedure is intended to write profiles of electron distribution function and spectral energy density of waves to the respective disk files in DOUBLE PRECISION 

  USE CONSTANTS
  USE PARAMETERS 
  USE MPI

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION (-n_v:, :), intent (in) :: fxyz
  REAL(KIND=8), DIMENSION (-n_kv:,:), intent (in) :: wxyz ,w_t,w_s,w_em, w_em_f
  REAL(KIND=8), DIMENSION (:), intent (in) :: omega_p,delta_x_p
  REAL(KIND=8), DIMENSION (n_xp+3) :: density, e_beam, e_wave
  integer ::  i, j
  integer, intent(in):: t_step
  REAL(KIND=8) :: norm_x, norm_dens, norm_eb, norm_ew
  REAL(KIND=8) ::  dens_tot, eb_tot, ew_tot
  character (LEN= 32):: vxfwd_flush, fxv_flush, wxv_flush, ts_flush, w_s_flush, w_em_flush,w_em_f_flush
  character (LEN= 5) :: Findex
  ! Findex - the output file numbers
   

density 	= 0.
e_beam 		= 0.
e_wave 		= 0.
dens_tot 	= 0.
eb_tot 		= 0.
ew_tot 		= 0.
!  Initialising some of the arrays

  write (Findex,'(i5.5)')  t_step
  Findex=trim(Findex)
  vxfwd_flush= trim(vxfwd_sequence//Findex//'.dat')  
  fxv_flush= trim(fxv_matrix//Findex//'.dat')  
  wxv_flush= trim(wxv_matrix//Findex//'.dat')  
  w_em_flush= trim(w_em_matrix//Findex//'.dat')  
  w_em_f_flush= trim(w_em_f_matrix//Findex//'.dat')  
  w_s_flush= trim(w_s_matrix//Findex//'.dat')  
 
! write electron distribution to files size n_v using MPI IO
  offset=0

IF (dumpArr(1) .GE. 1) THEN 
 CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//fxv_flush,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL, handle, ecode)
 CALL MPI_FILE_SET_VIEW(handle, offset, MPI_DOUBLE_PRECISION, mpitype, 'native', MPI_INFO_NULL, ecode)
 CALL MPI_FILE_WRITE(handle, Fxyz(:,3:n_xp+2), (2*n_v+1)*(n_xp), MPI_DOUBLE_PRECISION, status, ecode)
 CALL MPI_FILE_CLOSE(handle, ecode)
! 
ENDIF
 offset=0
! write wave distributions. sizes n_kv
IF (dumpArr(2) .GE. 1) THEN 
 CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//wxv_flush,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL, handle, ecode)
 CALL MPI_FILE_SET_VIEW(handle, offset, MPI_DOUBLE_PRECISION, mpitypeL, 'native', MPI_INFO_NULL, ecode)
 CALL MPI_FILE_WRITE(handle, wxyz(:,3:n_xp+2), (2*n_kv+1)*(n_xp), MPI_DOUBLE_PRECISION, status, ecode)
 CALL MPI_FILE_CLOSE(handle, ecode)
ENDIF

IF (dumpArr(3) .GE. 1) THEN 
 CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//w_s_flush,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL, handle, ecode)
 CALL MPI_FILE_SET_VIEW(handle, offset, MPI_DOUBLE_PRECISION, mpitypeL, 'native', MPI_INFO_NULL, ecode)
 CALL MPI_FILE_WRITE(handle, w_s(:,3:n_xp+2), (2*n_kv+1)*(n_xp), MPI_DOUBLE_PRECISION, status, ecode)
 CALL MPI_FILE_CLOSE(handle, ecode)
ENDIF

IF (dumpArr(4) .GE. 1) THEN 
 CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//w_em_flush,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL, handle, ecode)
 CALL MPI_FILE_SET_VIEW(handle, offset, MPI_DOUBLE_PRECISION, mpitypeL, 'native', MPI_INFO_NULL, ecode)
 CALL MPI_FILE_WRITE(handle, w_em(:,3:n_xp+2), (2*n_kv+1)*(n_xp), MPI_DOUBLE_PRECISION, status, ecode)
 CALL MPI_FILE_CLOSE(handle, ecode)
ENDIF

IF (dumpArr(5) .GE. 1) THEN 
 CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//w_em_f_flush,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL, handle, ecode)
 CALL MPI_FILE_SET_VIEW(handle, offset, MPI_DOUBLE_PRECISION, mpitypeL, 'native', MPI_INFO_NULL, ecode)
 CALL MPI_FILE_WRITE(handle, w_em_f(:,3:n_xp+2), (2*n_kv+1)*(n_xp), MPI_DOUBLE_PRECISION, status, ecode)
 CALL MPI_FILE_CLOSE(handle, ecode)
ENDIF


! calculate total beam and wave energies at each point in space etc------------------------------------------------------------------

DO j= 3, n_xp+2
  density(j) = sum(fxyz(:,j))*dv*k1
  DO i = -n_v, n_v
    IF (i .ne. 0) e_wave(j) = omega_p(j) * dv/(velocity(i)*velocity(i)) * wxyz(i,j) *k2 + e_wave(j)
    e_beam(j) = e_beam(j) + fxyz(i,j)*velocity(i)*velocity(i)
  ENDDO

  DO i = n_v+1, n_kv
    e_wave(j) = omega_p(j) / v_t *k2 *dk * w_t(i,j) + e_wave(j)
  ENDDO

  DO i = -n_kv, -n_v-1
    e_wave(j) = omega_p(j) / v_t *k2 *dk * w_t(i,j) + e_wave(j)
  ENDDO

ENDDO  
e_beam = e_beam*dv*k1*m_e/2.0

! write total energies to file---------------------------------------------------------------------------------------
offset = 0

  CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//vxfwd_flush,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL, handle, ecode)
  CALL MPI_FILE_SET_VIEW(handle, offset, MPI_DOUBLE_PRECISION, mpitype1, 'native', MPI_INFO_NULL, ecode)
  CALL MPI_FILE_WRITE(handle, omega_p(3:n_xp+2), n_xp, MPI_DOUBLE_PRECISION, status, ecode)
  CALL MPI_FILE_WRITE(handle, density(3:n_xp+2)/n_beam, n_xp, MPI_DOUBLE_PRECISION, status, ecode)
  CALL MPI_FILE_WRITE(handle, e_beam(3:n_xp+2), n_xp, MPI_DOUBLE_PRECISION, status, ecode)
  CALL MPI_FILE_WRITE(handle, e_wave(3:n_xp+2), n_xp, MPI_DOUBLE_PRECISION, status, ecode)

  CALL MPI_FILE_CLOSE(handle, ecode)


 density = density*delta_x_p
 e_beam  = e_beam*delta_x_p


!total of energies over all processors
 CALL MPI_ALLREDUCE(sum(density), dens_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ecode)
 CALL MPI_ALLREDUCE(sum(e_beam), eb_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ecode)
 CALL MPI_ALLREDUCE(sum(e_wave), ew_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ecode)

 IF(rank==0) WRITE(*,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') 'Density ',dens_tot,'   Beam E ',eb_tot,&
  '   Wave E ',ew_tot/ew_0,'   Tot E ',eb_tot+ew_tot
!  writing to screen


END SUBROUTINE  write_profiles

SUBROUTINE  write_profiles_single (fxyz, wxyz, w_t, w_s, w_em, w_em_f, t_step, omega_p,delta_x_p)
! The procedure is intended to write profiles of electron distribution function and spectral energy density of waves to the respective disk files in SINGLE PRECISION

  USE CONSTANTS
  USE PARAMETERS 
  USE MPI

  IMPLICIT NONE

  REAL(KIND=4), DIMENSION (-n_v:, :), intent (in) :: fxyz
  REAL(KIND=4), DIMENSION (-n_kv:,:), intent (in) :: wxyz ,w_t,w_s,w_em, w_em_f
  REAL(KIND=4), DIMENSION (:), intent (in) :: omega_p,delta_x_p
  REAL(KIND=4), DIMENSION (n_xp+3) :: density, e_beam, e_wave
  integer ::  i, j
  integer, intent(in):: t_step
  REAL(KIND=4) :: norm_x, norm_dens, norm_eb, norm_ew
  REAL(KIND=4) ::  dens_tot, eb_tot, ew_tot
  character (LEN= 32):: vxfwd_flush, fxv_flush, wxv_flush, ts_flush, w_s_flush, w_em_flush,w_em_f_flush
  character (LEN= 5) :: Findex
   

density 	= 0.
e_beam 		= 0.
e_wave 		= 0.
dens_tot 	= 0.
eb_tot 		= 0.
ew_tot 		= 0.
!  Initialising some of the arrays

  write (Findex,'(i5.5)')  t_step
  Findex=trim(Findex)
  vxfwd_flush= trim(vxfwd_sequence//Findex//'.dat')  
  fxv_flush= trim(fxv_matrix//Findex//'.dat')  
  wxv_flush= trim(wxv_matrix//Findex//'.dat')  
  w_em_flush= trim(w_em_matrix//Findex//'.dat')  
  w_em_f_flush= trim(w_em_f_matrix//Findex//'.dat')  
  w_s_flush= trim(w_s_matrix//Findex//'.dat')  

! write electron distribution to files size n_v using MPI IO
offset=0

IF (dumpArr(1) .GE. 1) THEN 
 CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//fxv_flush,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL, handle, ecode)
 CALL MPI_FILE_SET_VIEW(handle, offset, MPI_FLOAT, mpitypeF, 'native', MPI_INFO_NULL, ecode)
 CALL MPI_FILE_WRITE(handle, Fxyz(:,3:n_xp+2), (2*n_v+1)*(n_xp), MPI_FLOAT, status, ecode)
 CALL MPI_FILE_CLOSE(handle, ecode)
ENDIF
! write wave distributions. sizes n_kv

IF (dumpArr(2) .GE. 1) THEN  
 CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//wxv_flush,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL, handle, ecode)
 CALL MPI_FILE_SET_VIEW(handle, offset, MPI_FLOAT, mpitypeLF, 'native', MPI_INFO_NULL, ecode)
 CALL MPI_FILE_WRITE(handle, wxyz(:,3:n_xp+2), (2*n_kv+1)*(n_xp), MPI_FLOAT, status, ecode)
 CALL MPI_FILE_CLOSE(handle, ecode)
ENDIF

IF (dumpArr(3) .GE. 1) THEN 
 CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//w_s_flush,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL, handle, ecode)
 CALL MPI_FILE_SET_VIEW(handle, offset, MPI_FLOAT, mpitypeLF, 'native', MPI_INFO_NULL, ecode)
 CALL MPI_FILE_WRITE(handle, w_s(:,3:n_xp+2), (2*n_kv+1)*(n_xp), MPI_FLOAT, status, ecode)
 CALL MPI_FILE_CLOSE(handle, ecode)
ENDIF

IF (dumpArr(4) .GE. 1) THEN 
 CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//w_em_flush,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL, handle, ecode)
 CALL MPI_FILE_SET_VIEW(handle, offset, MPI_FLOAT, mpitypeLF, 'native', MPI_INFO_NULL, ecode)
 CALL MPI_FILE_WRITE(handle, w_em(:,3:n_xp+2), (2*n_kv+1)*(n_xp), MPI_FLOAT, status, ecode)
 CALL MPI_FILE_CLOSE(handle, ecode)
ENDIF

IF (dumpArr(5) .GE. 1) THEN 
 CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//w_em_f_flush,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL, handle, ecode)
 CALL MPI_FILE_SET_VIEW(handle, offset, MPI_FLOAT, mpitypeLF, 'native', MPI_INFO_NULL, ecode)
 CALL MPI_FILE_WRITE(handle, w_em_f(:,3:n_xp+2), (2*n_kv+1)*(n_xp), MPI_FLOAT, status, ecode)
 CALL MPI_FILE_CLOSE(handle, ecode)
ENDIF


! calculate total beam and wave energies at each point in space etc------------------------------------------------------------------

DO j= 3, n_xp+2
  density(j) = sum(fxyz(:,j))*dv*k1
  DO i = -n_v, n_v
    IF (i .ne. 0) e_wave(j) = omega_p(j) * dv/(velocity(i)*velocity(i)) * wxyz(i,j) *k2 + e_wave(j)
    e_beam(j) = e_beam(j) + fxyz(i,j)*velocity(i)*velocity(i)
  ENDDO

  DO i = n_v+1, n_kv
    e_wave(j) = omega_p(j) / v_t *k2 *dk * w_t(i,j) + e_wave(j)
  ENDDO

  DO i = -n_kv, -n_v-1
    e_wave(j) = omega_p(j) / v_t *k2 *dk * w_t(i,j) + e_wave(j)
  ENDDO

ENDDO  
e_beam = e_beam*dv*k1*m_e/2.0

! write total energies to file---------------------------------------------------------------------------------------

  CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//vxfwd_flush,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL, handle, ecode)
  CALL MPI_FILE_SET_VIEW(handle, offset, MPI_FLOAT, mpitype1F, 'native', MPI_INFO_NULL, ecode)
  CALL MPI_FILE_WRITE(handle, omega_p(3:n_xp+2), n_xp, MPI_FLOAT, status, ecode)
  CALL MPI_FILE_WRITE(handle, density(3:n_xp+2)/n_beam, n_xp, MPI_FLOAT, status, ecode)
  CALL MPI_FILE_WRITE(handle, e_beam(3:n_xp+2), n_xp, MPI_FLOAT, status, ecode)
  CALL MPI_FILE_WRITE(handle, e_wave(3:n_xp+2), n_xp, MPI_FLOAT, status, ecode)

  CALL MPI_FILE_CLOSE(handle, ecode)

 density = density*delta_x_p
 e_beam  = e_beam*delta_x_p

!total of energies over all processors
 CALL MPI_ALLREDUCE(sum(density), dens_tot, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, ecode)
 CALL MPI_ALLREDUCE(sum(e_beam), eb_tot, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, ecode)
 CALL MPI_ALLREDUCE(sum(e_wave), ew_tot, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, ecode)

 IF(rank==0) WRITE(*,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') 'Density ',dens_tot,'   Beam E ',eb_tot,&
  '   Wave E ',ew_tot/ew_0,'   Tot E ',eb_tot+ew_tot
!  writing to screen


END SUBROUTINE  write_profiles_single


SUBROUTINE  write_profiles_restart (fxyz, wxyz, w_t, w_s, w_em, w_em_f, t_step, omega_p,delta_x_p)
! Dump full precision data and everything to file every restart_t

  USE CONSTANTS
  USE PARAMETERS 
  USE MPI

  IMPLICIT NONE

  REAL(KIND=8), DIMENSION (-n_v:, :), intent (in) :: fxyz
  REAL(KIND=8), DIMENSION (-n_kv:,:), intent (in) :: wxyz ,w_t,w_s,w_em, w_em_f
  REAL(KIND=8), DIMENSION (:), intent (in) :: omega_p,delta_x_p
  integer, intent(in):: t_step
  character (LEN= 32):: vxfwd_flush, fxv_flush, wxv_flush, ts_flush, w_s_flush, w_em_flush,w_em_f_flush
  character (LEN= 5) :: Findex
  ! Findex - the output file numbers

  write (Findex,'(i5.5)')  t_step
  Findex=trim(Findex)
  vxfwd_flush= trim(vxfwd_sequence//Findex//'R.dat')  
  fxv_flush= trim(fxv_matrix//Findex//'R.dat')  
  wxv_flush= trim(wxv_matrix//Findex//'R.dat')  
  w_em_flush= trim(w_em_matrix//Findex//'R.dat')  
  w_em_f_flush= trim(w_em_f_matrix//Findex//'R.dat')  
  w_s_flush= trim(w_s_matrix//Findex//'R.dat')  
 


! write electron distribution to files size n_v using MPI IO
  offset = 0
 CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//fxv_flush,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL, handle, ecode)
 CALL MPI_FILE_SET_VIEW(handle, offset, MPI_DOUBLE_PRECISION, mpitype, 'native', MPI_INFO_NULL, ecode)
 CALL MPI_FILE_WRITE(handle, Fxyz(:,3:n_xp+2), (2*n_v+1)*(n_xp), MPI_DOUBLE_PRECISION, status, ecode)
 CALL MPI_FILE_CLOSE(handle, ecode)
! 

! write wave distributions. sizes n_kv

 CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//wxv_flush,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL, handle, ecode)
 CALL MPI_FILE_SET_VIEW(handle, offset, MPI_DOUBLE_PRECISION, mpitypeL, 'native', MPI_INFO_NULL, ecode)
 CALL MPI_FILE_WRITE(handle, wxyz(:,3:n_xp+2), (2*n_kv+1)*(n_xp), MPI_DOUBLE_PRECISION, status, ecode)
 CALL MPI_FILE_CLOSE(handle, ecode)
! 
 CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//w_s_flush,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL, handle, ecode)
 CALL MPI_FILE_SET_VIEW(handle, offset, MPI_DOUBLE_PRECISION, mpitypeL, 'native', MPI_INFO_NULL, ecode)
 CALL MPI_FILE_WRITE(handle, w_s(:,3:n_xp+2), (2*n_kv+1)*(n_xp), MPI_DOUBLE_PRECISION, status, ecode)
 CALL MPI_FILE_CLOSE(handle, ecode)

 CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//w_em_flush,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL, handle, ecode)
 CALL MPI_FILE_SET_VIEW(handle, offset, MPI_DOUBLE_PRECISION, mpitypeL, 'native', MPI_INFO_NULL, ecode)
 CALL MPI_FILE_WRITE(handle, w_em(:,3:n_xp+2), (2*n_kv+1)*(n_xp), MPI_DOUBLE_PRECISION, status, ecode)
 CALL MPI_FILE_CLOSE(handle, ecode)

 CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//w_em_f_flush,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL, handle, ecode)
 CALL MPI_FILE_SET_VIEW(handle, offset, MPI_DOUBLE_PRECISION, mpitypeL, 'native', MPI_INFO_NULL, ecode)
 CALL MPI_FILE_WRITE(handle, w_em_f(:,3:n_xp+2), (2*n_kv+1)*(n_xp), MPI_DOUBLE_PRECISION, status, ecode)
 CALL MPI_FILE_CLOSE(handle, ecode)


  IF(rollRestart .GE. 1) THEN
  	IF(rank == 0) THEN
  !if we're rolling restarts, then we delete the second to last one.
  !We ALWAYS keep the current one, and the previous, just in case we manage to terminate mid-dump...
    
    
    write (Findex,'(i5.5)')  t_step-2
    IF(t_step .LT. 2) write (Findex,'(i5.5)')  0
    Findex=trim(Findex)
	print*, 'Deleting Restart no.', findex
    vxfwd_flush= trim(vxfwd_sequence//Findex//'R.dat')  
    fxv_flush= trim(fxv_matrix//Findex//'R.dat')  
    wxv_flush= trim(wxv_matrix//Findex//'R.dat')  
    w_em_flush= trim(w_em_matrix//Findex//'R.dat')  
    w_em_f_flush= trim(w_em_f_matrix//Findex//'R.dat')  
    w_s_flush= trim(w_s_matrix//Findex//'R.dat')  

 
    CALL MPI_FILE_DELETE(TRIM(path)//fxv_flush, MPI_INFO_NULL, ecode)
	CALL MPI_FILE_DELETE(TRIM(path)//wxv_flush, MPI_INFO_NULL, ecode)
  	CALL MPI_FILE_DELETE(TRIM(path)//w_s_flush, MPI_INFO_NULL, ecode)
	CALL MPI_FILE_DELETE(TRIM(path)//w_em_flush, MPI_INFO_NULL, ecode)
	CALL MPI_FILE_DELETE(TRIM(path)//w_em_f_flush, MPI_INFO_NULL, ecode)
 
	ENDIF
  ENDIF


END SUBROUTINE  write_profiles_restart


SUBROUTINE write_density(r_xi,dens_xi,dens_smooth_x, wind_x,dn_dx, omega_x, l_x)

USE CONSTANTS
USE PARAMETERS
USE MPI

IMPLICIT NONE

 REAL(KIND=8), DIMENSION(:), INTENT(IN)::r_xi, dens_xi, wind_x, dn_dx, dens_smooth_x, omega_x, l_x
   
  offset = 0
 
  CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(path)//corona_file,MPI_MODE_CREATE+MPI_MODE_WRONLY,MPI_INFO_NULL, handle, ecode)
  CALL MPI_FILE_SET_VIEW(handle, offset, MPI_DOUBLE_PRECISION, mpitype1, 'native', MPI_INFO_NULL, ecode)
  CALL MPI_FILE_WRITE(handle, (r_xi(3:n_xp+2) + r_s + x_0)/r_s, n_xp, MPI_DOUBLE_PRECISION, status, ecode)
  CALL MPI_FILE_WRITE(handle, dens_xi(3:n_xp+2), n_xp, MPI_DOUBLE_PRECISION, status, ecode)
  CALL MPI_FILE_WRITE(handle, omega_x(3:n_xp+2), n_xp, MPI_DOUBLE_PRECISION, status, ecode)
  CALL MPI_FILE_WRITE(handle, wind_x(3:n_xp+2), n_xp, MPI_DOUBLE_PRECISION, status, ecode)
  CALL MPI_FILE_WRITE(handle, dens_smooth_x(3:n_xp+2), n_xp, MPI_DOUBLE_PRECISION, status, ecode)
  CALL MPI_FILE_WRITE(handle, dn_dx(3:n_xp+2), n_xp, MPI_DOUBLE_PRECISION, status, ecode)
  CALL MPI_FILE_WRITE(handle, l_x(3:n_xp+2), n_xp, MPI_DOUBLE_PRECISION, status, ecode)

  CALL MPI_FILE_CLOSE(handle, ecode)

! omega_x(3:n_xp+2)/(2.*pi*1.e6)
! wind_x(3:n_xp+2)/1.e5
 
END SUBROUTINE write_density

SUBROUTINE write_dynSpec(w_em, w_em_f,t_step,omega, delta_x)

  USE CONSTANTS
  USE PARAMETERS

  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(-n_kv:, :):: w_em, w_em_f
  REAL(KIND=8), DIMENSION(:):: omega, delta_x
  INTEGER(KIND=4):: t_step
  CHARACTER (LEN= 32):: filename
  CHARACTER (LEN= 5) :: Findex

  write (Findex,'(i5.5)')  t_step
  Findex = trim(Findex)
  filename = trim(dynspecf_file//Findex//'.dat')

  CALL MPI_FILE_OPEN(MPI_COMM_WORLD, trim(path)//filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, handle, ecode)
  CALL MPI_FILE_SET_VIEW(handle, offset, MPI_DOUBLE_PRECISION, mpitype1, 'native', MPI_INFO_NULL, ecode)
  CALL MPI_FILE_WRITE(handle, w_em_f(-n_kv+2,:), n_xp, MPI_DOUBLE_PRECISION, status, ecode)
  CALL MPI_FILE_WRITE(handle, w_em_f(n_kv-2,:), n_xp, MPI_DOUBLE_PRECISION, status, ecode)
  CALL MPI_FILE_CLOSE(handle, ecode)

  filename = trim(dynspech_file//Findex//'.dat')

  CALL MPI_FILE_OPEN(MPI_COMM_WORLD, trim(path)//filename, MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, handle, ecode)
  CALL MPI_FILE_SET_VIEW(handle, offset, MPI_DOUBLE_PRECISION, mpitype1, 'native', MPI_INFO_NULL, ecode)
  CALL MPI_FILE_WRITE(handle, w_em(-n_kv+2,:), n_xp, MPI_DOUBLE_PRECISION, status, ecode)
  CALL MPI_FILE_WRITE(handle, w_em(n_kv-2,:), n_xp, MPI_DOUBLE_PRECISION, status, ecode)
  CALL MPI_FILE_CLOSE(handle, ecode)

END SUBROUTINE write_dynSpec

SUBROUTINE write_parameters   

  USE CONSTANTS
  USE PARAMETERS

  IMPLICIT NONE



  write(*,*) '------------------------------------------'
  write(*,*) 'Finite diference scheme parameters :' 
  write(*,'(A,F14.10,A)')'Time step dt		=', dt,  ' sec'
!   write(*,'(A,F6.4,A)')'X-coodinate step dx	=', dx/d,'*d cm'
  write(*,'(A,F6.2,A)')'Velocity step dv	=', dv/v_beam,'*Vbeam cm/sec'
  write(*,*) '------------------------------------------'
  Write(*,*)'ELECTRON BEAM & STARTING PLASMA PARAMETERS:'
  write(*,'(A,F9.1,A)')'Plasma frequency	=', Omega_0/(2*PI*1000000), ' MHz'
  write(*,'(A,F9.1,A)')'Thermal velocity	=', v_t/1e+9, '*10^9 cm/s'
  write(*,'(A,ES9.2,A)')'Beam density		=', n_beam, ' cm^{-3}'
!   write(*,'(A,F9.5,A)')'Maximum distance	=', r_x/1e+11,'*10^11 cm'
!   write(*,'(A,F9.5,A)')'Time Constant		= ',time_const,tqv(25),(n_v-1)**2
  write(*,*) '------------------------------------------'
  write(*,*) "Calculations started		OK";

END SUBROUTINE write_parameters

SUBROUTINE write_velocity_params

  USE CONSTANTS
  USE PARAMETERS

  IMPLICIT NONE

  INTEGER:: WriteStat

  OPEN (UNIT=20, FILE = trim(path)//params_file, STATUS = 'REPLACE', ACTION='WRITE' ) 

  WRITE (20,*, IOSTAT=WRITESTAT) v_min/v_beam
  WRITE (20,*, IOSTAT=WRITESTAT) v_max/v_beam
  WRITE (20,*, IOSTAT=WRITESTAT) n_kv*2 + 1
  WRITE (20,*, IOSTAT=WRITESTAT) 2*n_v + 1
  WRITE (20,*, IOSTAT=WRITESTAT) n_x
  WRITE (20,*, IOSTAT=WRITESTAT) time_save
  WRITE (20,*, IOSTAT=WRITESTAT) n_beam
  WRITE (20,*, IOSTAT=WRITESTAT) d
  WRITE (20,*, IOSTAT=WRITESTAT) x_0
  WRITE (20,*, IOSTAT=WRITESTAT) x_min

  CLOSE(20)

  OPEN (UNIT=20, FILE = trim(path)//velocity_file, STATUS = 'REPLACE', ACTION='WRITE' ) 

  WRITE(20,*, IOSTAT=WRITESTAT) velocity
  WRITE(20,*, IOSTAT=WRITESTAT) k_x
  WRITE(20,*, IOSTAT=WRITESTAT) k_f
  WRITE(20,*, IOSTAT=WRITESTAT) k_h

  CLOSE(20)

END SUBROUTINE write_velocity_params

SUBROUTINE write_progress(time_beg, time_prev)

!   writes progress file at each data save

 USE PARAMETERS
 USE CONSTANTS
 IMPLICIT NONE

  INTEGER::WStat
  REAL(KIND=8), INTENT(IN):: time_beg
  REAL(KIND=8), INTENT(INOUT):: time_prev
  REAL(KIND=8):: time_now

  time_now=MPI_WTIME()

  WRITE (*,'(A, F6.2,A,F6.2,A)') 'Code t =', t, 's. Time elapsed =', time_now-time_prev, 's'
  
  OPEN(UNIT=33, FILE= 'progress.dat', STATUS='OLD', POSITION='APPEND', ACTION='WRITE')
    WRITE(33, '(A, F7.2,A, F10.0,A, F7.2,A)', IOSTAT=WStat) 'Code time = ', &
  t, 's. Elapsed =', time_now - time_beg, 's. Last it =', time_now-time_prev, 's'
    
  CLOSE(33)

  time_prev=time_now

END SUBROUTINE write_progress

SUBROUTINE write_warnings(clip_w)

!   writes progress file at each data save

 USE PARAMETERS
 USE CONSTANTS
 IMPLICIT NONE

  INTEGER::WStat
  INTEGER(KIND=4), INTENT(IN):: clip_w
  
  WRITE(*,*) 'WARNING, WARNING, clipping timestep'

  
  OPEN(UNIT=33, FILE= 'warnings.dat', STATUS='OLD', POSITION='APPEND', ACTION='WRITE')
    WRITE(33, '(A, F7.2,A)', IOSTAT=WStat) 'Code time = ', t, 'Clipped timesteps'
    
  CLOSE(33)

 
END SUBROUTINE write_warnings


SUBROUTINE write_total_time(time_beg)
!   writes total wallclock time file of main code

  USE PARAMETERS
  USE CONSTANTS
  IMPLICIT NONE

  INTEGER::WStat
  REAL(KIND=8), INTENT(IN):: time_beg
  REAL(KIND=8):: time_end

  time_end=MPI_WTIME()
  OPEN(UNIT=33, FILE= 'wallclock.dat', STATUS='REPLACE', ACTION='WRITE')

  WRITE(33, *, IOSTAT=WStat) 'Total time elapsed = ', time_end - time_beg, 's'
  CLOSE(33)


END SUBROUTINE write_total_time


SUBROUTINE readInit
 
 USE PARAMETERS
 USE CONSTANTS
 IMPLICIT NONE

 INTEGER::ReadStatus, RSall
 ! service variable to control read status

  NAMELIST/config/T_max,  d,  x_0, x_min, time_print, time_save, restart_t,dyn_spec_t,t_start,  n_x
  NAMELIST/beam/n_beam, v_beam, tau, tau2
  NAMELIST/plasma/N, t_e, t_i
  NAMELIST/filenames/path
  NAMELIST/datadumps/dumpArr, dumpSingle, rollRestart

  RSall=0
  
  OPEN (UNIT=10, FILE= init_params, STATUS='OLD', ACTION='READ' ) 
   ! open file for reading
   
  READ (10,NML=config, IOSTAT= ReadStatus ) 
  RSall=RSall + ReadStatus
  IF(ReadStatus .NE. 0)  print*, 'Error reading input file in block "config"'
  READ (10,NML=beam, IOSTAT= ReadStatus ) 
  RSall=RSall + ReadStatus
  IF(ReadStatus .NE. 0)  print*, 'Error reading input file in block "beam"'
  READ (10,NML=plasma, IOSTAT= ReadStatus )
  RSall=RSall + ReadStatus
  IF(ReadStatus .NE. 0)  print*, 'Error reading input file in block "plasma"'
  READ (10,NML=filenames, IOSTAT= ReadStatus ) 
  RSall=RSall + ReadStatus
  IF(ReadStatus .NE. 0)  print*, 'Error reading input file in block "filenames"'
  READ (10,NML=datadumps, IOSTAT= ReadStatus ) 
  RSall=RSall + ReadStatus
  IF(ReadStatus .NE. 0)  print*, 'Error reading input file in block "datadumps"'
  
  CLOSE(10)

  IF(RSall .NE. 0) THEN
    print*, 'Error reading input file, aborting'
    stop
  ENDIF

  OPEN(UNIT=33, FILE= 'progress.dat', STATUS='REPLACE', ACTION='WRITE')
    WRITE(33, *) 'Beginning calculations'
  CLOSE(33)

  OPEN(UNIT=33, FILE= 'warnings.dat', STATUS='REPLACE', ACTION='WRITE')
    WRITE(33, *) 'Beginning calculations'
  CLOSE(33)

! IF YOU ADD ANYTHING HERE, MAKE SURE TO ADD A BROADCAST TOO......  

END SUBROUTINE ReadInit


SUBROUTINE restart(fxyz, wxyz, w_s, w_em, w_em_f, time)

 USE PARAMETERS
 USE CONSTANTS
 USE MPI
 IMPLICIT NONE



  REAL(KIND=8), DIMENSION (-n_v:, :), intent (out) :: fxyz
  REAL(KIND=8), DIMENSION (-n_kv:,:), intent (out) :: wxyz ,w_s,w_em, w_em_f
  INTEGER(KIND=4):: time
  character (LEN= 5) :: Findex
  character (LEN= 32):: vxfwd_flush, fxv_flush, wxv_flush, ts_flush, w_s_flush, w_em_flush,w_em_f_flush
  ! Findex - the output file numbers
   

  IF(rank == 0) print*, 'READING OLD DATA FILES:'

  write (Findex,'(i5.5)')  time
  Findex=trim(Findex)
  vxfwd_flush= trim(vxfwd_sequence//Findex//'R.dat')  
  fxv_flush= trim(fxv_matrix//Findex//'R.dat')  
  wxv_flush= trim(wxv_matrix//Findex//'R.dat')  
  w_em_flush= trim(w_em_matrix//Findex//'R.dat')  
  w_em_f_flush= trim(w_em_f_matrix//Findex//'R.dat')  
  w_s_flush= trim(w_s_matrix//Findex//'R.dat')  
 
! read electron distribution to files size n_v using MPI IO

  offset = 0
 CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//fxv_flush,MPI_MODE_RDONLY,MPI_INFO_NULL, handle, ecode)
 CALL MPI_FILE_SET_VIEW(handle, offset, MPI_DOUBLE_PRECISION, mpitype, 'native', MPI_INFO_NULL, ecode)
 CALL MPI_FILE_READ(handle, Fxyz(:,3:n_xp+2), (2*n_v+1)*(n_xp), MPI_DOUBLE_PRECISION, status, ecode)
 CALL MPI_FILE_CLOSE(handle, ecode)
  print*, rank, 'f', status, ecode

! read wave distributions. sizes n_kv

 
 CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//wxv_flush,MPI_MODE_RDONLY,MPI_INFO_NULL, handle, ecode)
 CALL MPI_FILE_SET_VIEW(handle, offset, MPI_DOUBLE_PRECISION, mpitypeL, 'native', MPI_INFO_NULL, ecode)
 CALL MPI_FILE_READ(handle, wxyz(:,3:n_xp+2), (2*n_kv+1)*(n_xp), MPI_DOUBLE_PRECISION, status, ecode)
 CALL MPI_FILE_CLOSE(handle, ecode)
  print*, rank, 'w', status, ecode

! 
 CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//w_s_flush,MPI_MODE_RDONLY,MPI_INFO_NULL, handle, ecode)
 CALL MPI_FILE_SET_VIEW(handle, offset, MPI_DOUBLE_PRECISION, mpitypeL, 'native', MPI_INFO_NULL, ecode)
 CALL MPI_FILE_READ(handle, w_s(:,3:n_xp+2), (2*n_kv+1)*(n_xp), MPI_DOUBLE_PRECISION, status, ecode)
 CALL MPI_FILE_CLOSE(handle, ecode)
  print*, rank, 'w_s',status, ecode


 CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//w_em_flush,MPI_MODE_RDONLY,MPI_INFO_NULL, handle, ecode)
 CALL MPI_FILE_SET_VIEW(handle, offset, MPI_DOUBLE_PRECISION, mpitypeL, 'native', MPI_INFO_NULL, ecode)
 CALL MPI_FILE_READ(handle, w_em(:,3:n_xp+2), (2*n_kv+1)*(n_xp), MPI_DOUBLE_PRECISION, status, ecode)
 CALL MPI_FILE_CLOSE(handle, ecode)
  print*, rank, 'w_em',status, ecode


 CALL MPI_FILE_OPEN(MPI_COMM_WORLD,TRIM(path)//w_em_f_flush,MPI_MODE_RDONLY,MPI_INFO_NULL, handle, ecode)
 CALL MPI_FILE_SET_VIEW(handle, offset, MPI_DOUBLE_PRECISION, mpitypeL, 'native', MPI_INFO_NULL, ecode)
 CALL MPI_FILE_READ(handle, w_em_f(:,3:n_xp+2), (2*n_kv+1)*(n_xp), MPI_DOUBLE_PRECISION, status, ecode)
 CALL MPI_FILE_CLOSE(handle, ecode)
  print*, rank, 'w_em_f',status, ecode

    

    CALL MPI_SENDRECV(fxyz(:,n_xp+1:n_xp+2),(2*n_v+1)*2,MPI_DOUBLE_PRECISION, right, tag, &
            fxyz(:,1:2),(2*n_v+1)*2,MPI_DOUBLE_PRECISION, left, tag,MPI_COMM_WORLD, status, ecode)

    CALL MPI_SENDRECV(fxyz(:,3),(2*n_v+1),MPI_DOUBLE_PRECISION, left, tag, &
            fxyz(:,n_xp+3),(2*n_v+1),MPI_DOUBLE_PRECISION, right, tag,MPI_COMM_WORLD, status, ecode)


    CALL MPI_SENDRECV(wxyz(:,n_xp+1:n_xp+2),(2*n_kv+1)*2,MPI_DOUBLE_PRECISION, right, tag, &
            wxyz(:,1:2),(2*n_kv+1)*2,MPI_DOUBLE_PRECISION, left, tag,MPI_COMM_WORLD, status, ecode)

    CALL MPI_SENDRECV(wxyz(:,3),(2*n_kv+1),MPI_DOUBLE_PRECISION, left, tag, &
            wxyz(:,n_xp+3),(2*n_kv+1),MPI_DOUBLE_PRECISION, right, tag,MPI_COMM_WORLD, status, ecode)


    CALL MPI_SENDRECV(w_em(:,n_xp+1:n_xp+2),(2*n_kv+1)*2,MPI_DOUBLE_PRECISION, right, tag, &
            w_em(:,1:2),(2*n_kv+1)*2,MPI_DOUBLE_PRECISION, left, tag,MPI_COMM_WORLD, status, ecode)

    CALL MPI_SENDRECV(w_em(:,3),(2*n_kv+1),MPI_DOUBLE_PRECISION, left, tag, &
            w_em(:,n_xp+3),(2*n_kv+1),MPI_DOUBLE_PRECISION, right, tag,MPI_COMM_WORLD, status, ecode)


    CALL MPI_SENDRECV(w_em_f(:,n_xp+1:n_xp+2),(2*n_kv+1)*2,MPI_DOUBLE_PRECISION, right, tag, &
            w_em_f(:,1:2),(2*n_kv+1)*2,MPI_DOUBLE_PRECISION, left, tag,MPI_COMM_WORLD, status, ecode)

    CALL MPI_SENDRECV(w_em_f(:,3),(2*n_kv+1),MPI_DOUBLE_PRECISION, left, tag, &
            w_em_f(:,n_xp+3),(2*n_kv+1),MPI_DOUBLE_PRECISION, right, tag,MPI_COMM_WORLD, status, ecode)


END SUBROUTINE restart

SUBROUTINE restart_density(dens_xi, dens_smooth_x, omega_x, wind_x, R_xi, dn_dx_in, delta_x_in, l_in)


 USE PARAMETERS
 USE CONSTANTS
 USE MPI
 IMPLICIT NONE

  REAL(KIND=8), DIMENSION(:), INTENT(OUT):: dens_xi, wind_x, dn_dx_in, dens_smooth_x, omega_x, l_in
  REAL(KIND=8), DIMENSION(:)::r_xi, delta_x_in
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE::r_xi_dummy
  
  ALLOCATE(r_xi_dummy(n_xp+3))

  offset = 0
   
  CALL MPI_FILE_OPEN(MPI_COMM_WORLD,trim(path)//corona_file,MPI_MODE_RDONLY,MPI_INFO_NULL, handle, ecode)
  CALL MPI_FILE_SET_VIEW(handle, offset, MPI_DOUBLE_PRECISION, mpitype1, 'native', MPI_INFO_NULL, ecode)
  CALL MPI_FILE_READ(handle, r_xi_dummy(3:n_xp+2), n_xp, MPI_DOUBLE_PRECISION, status, ecode)
  CALL MPI_FILE_READ(handle, dens_xi(3:n_xp+2), n_xp, MPI_DOUBLE_PRECISION, status, ecode)
  CALL MPI_FILE_READ(handle, omega_x(3:n_xp+2), n_xp, MPI_DOUBLE_PRECISION, status, ecode)
  CALL MPI_FILE_READ(handle, wind_x(3:n_xp+2), n_xp, MPI_DOUBLE_PRECISION, status, ecode)
  CALL MPI_FILE_READ(handle, dens_smooth_x(3:n_xp+2), n_xp, MPI_DOUBLE_PRECISION, status, ecode)
  CALL MPI_FILE_READ(handle, dn_dx_in(3:n_xp+2), n_xp, MPI_DOUBLE_PRECISION, status, ecode)
  CALL MPI_FILE_READ(handle, l_in(3:n_xp+2), n_xp, MPI_DOUBLE_PRECISION, status, ecode)
  
  CALL MPI_FILE_CLOSE(handle, ecode)

     CALL MPI_SENDRECV(dens_xi(n_xp+1:n_xp+2),2,MPI_DOUBLE_PRECISION, right, tag, &
            dens_xi(1:2),2,MPI_DOUBLE_PRECISION, left, tag,MPI_COMM_WORLD, status, ecode)

    CALL MPI_SENDRECV(dens_xi(3),1,MPI_DOUBLE_PRECISION, left, tag, &
            dens_xi(n_xp+3),1,MPI_DOUBLE_PRECISION, right, tag,MPI_COMM_WORLD, status, ecode)
    
    CALL MPI_SENDRECV(omega_x(n_xp+1:n_xp+2),2,MPI_DOUBLE_PRECISION, right, tag, &
            omega_x(1:2),2,MPI_DOUBLE_PRECISION, left, tag,MPI_COMM_WORLD, status, ecode)

    CALL MPI_SENDRECV(omega_x(3),1,MPI_DOUBLE_PRECISION, left, tag, &
            omega_x(n_xp+3),1,MPI_DOUBLE_PRECISION, right, tag,MPI_COMM_WORLD, status, ecode)

    CALL MPI_SENDRECV(wind_x(n_xp+1:n_xp+2),2,MPI_DOUBLE_PRECISION, right, tag, &
            wind_x(1:2),2,MPI_DOUBLE_PRECISION, left, tag,MPI_COMM_WORLD, status, ecode)

    CALL MPI_SENDRECV(wind_x(3),1,MPI_DOUBLE_PRECISION, left, tag, &
            wind_x(n_xp+3),1,MPI_DOUBLE_PRECISION, right, tag,MPI_COMM_WORLD, status, ecode)

    CALL MPI_SENDRECV(dens_smooth_x(n_xp+1:n_xp+2),2,MPI_DOUBLE_PRECISION, right, tag, &
            dens_smooth_x(1:2),2,MPI_DOUBLE_PRECISION, left, tag,MPI_COMM_WORLD, status, ecode)

    CALL MPI_SENDRECV(dens_smooth_x(3),1,MPI_DOUBLE_PRECISION, left, tag, &
            dens_smooth_x(n_xp+3),1,MPI_DOUBLE_PRECISION, right, tag,MPI_COMM_WORLD, status, ecode)

    CALL MPI_SENDRECV(dn_dx_in(n_xp+1:n_xp+2),2,MPI_DOUBLE_PRECISION, right, tag, &
            dn_dx_in(1:2),2,MPI_DOUBLE_PRECISION, left, tag,MPI_COMM_WORLD, status, ecode)

    CALL MPI_SENDRECV(dn_dx_in(3),1,MPI_DOUBLE_PRECISION, left, tag, &
            dn_dx_in(n_xp+3),1,MPI_DOUBLE_PRECISION, right, tag,MPI_COMM_WORLD, status, ecode)

    CALL MPI_SENDRECV(l_in(n_xp+1:n_xp+2),2,MPI_DOUBLE_PRECISION, right, tag, &
            l_in(1:2),2,MPI_DOUBLE_PRECISION, left, tag,MPI_COMM_WORLD, status, ecode)

    CALL MPI_SENDRECV(l_in(3),1,MPI_DOUBLE_PRECISION, left, tag, &
            l_in(n_xp+3),1,MPI_DOUBLE_PRECISION, right, tag,MPI_COMM_WORLD, status, ecode)



  DEALLOCATE(r_xi_dummy)
 
END SUBROUTINE restart_density



END MODULE I_O
