MODULE PARAMETERS

! Various parameters used in many locations or across subroutines

USE CONSTANTS
USE MPI

IMPLICIT NONE

REAL(KIND=8) :: k1, k2, a1
! normalisation coeficients
REAL(KIND=8) :: v, x, dx, t, dt,  dv, t_0, dvdv,dk
! intrinsic variables 
REAL(KIND=8) :: ew_0
! intial Langmuir wave energy
REAL(KIND=8) :: dt_ql_const,dt_ql,dt_ql_av
REAL(KIND=8) :: n_e, n_b, e_b
	! plasma density, beam density and beam energy density
REAL(KIND=8):: alpha_is,gamma_s, k_s_stretch, gamma_d_emf

REAL(KIND=8):: omega_0
! reference plasma frequency

! variables read from input file see InitL.par............................................................

  REAL(KIND=8):: T_max, d, x_0, x_min, time_save, t_start,restart_t, dyn_spec_t, time_print
  INTEGER(KIND=4):: n_x
  REAL(KIND=8):: n_beam, v_beam, tau, tau2
  INTEGER(KIND=4):: N
  REAL(KIND=8):: t_e, t_i
  CHARACTER(LEN=50):: path
  INTEGER, DIMENSION(6):: dumpArr

! ----------------------------------------------------------------------------------------------

REAL(KIND=8):: v_t,three_v_t_sq,v_ti,v_s, kbTe_1d, gamma_d_em, alpha_em_h, alpha_em_f, min_dx, dk_h, max_l
! Thermal electron beam velocity

INTEGER:: n_xp, v_beam_i, n_p, lpc, lpe
  

!   MPI service variables.........................................................................
  INTEGER(KIND=4) :: ecode,rank
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER :: left = MPI_PROC_NULL, right=MPI_PROC_NULL, tag
  INTEGER, DIMENSION(2) :: starts, sizes, subsizes
  INTEGER, DIMENSION(1) :: starts1, sizes1, subsizes1
  INTEGER :: mpitype, handle, mpitypeL, mpitype1, mpitypeF, mpitypeLF, mpitype1F, dumpSingle, rollRestart
  INTEGER(KIND=MPI_OFFSET_KIND) :: offset

! --------------------------------------------------------------------------------------------------

REAL(KIND=8), DIMENSION(-n_v:n_v):: velocity
REAL(KIND=8), DIMENSION(-n_kv:n_kv):: k_x, dk_x, k_f, k_h
!velocity and wavenumber
INTEGER(KIND=4), DIMENSION(-n_kv:n_kv):: ind_2pl,ind_pls, ind_fs, ind_fe, ind_fe2,ind_h1, ind_h2,up_2pl
INTEGER(KIND=4), DIMENSION(-n_kv:n_kv):: up_pls, up_fe, up_fe2, dn_2pl, dn_pls
REAL(KIND=8), DIMENSION(-n_kv: n_kv, 3):: frac_2pl, frac_pls
REAL(KIND=8), DIMENSION(-n_kv:n_kv):: frac_fe, frac_fe2
!     wave wave coefficients

END MODULE PARAMETERS