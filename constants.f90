 Module CONSTANTS

 ! This module contains constant of calculations required for later programs
    
	REAL(KIND=8), PARAMETER:: PI=3.141593
	! Pi

	REAL(KIND=8), PARAMETER:: t_min = 0;
	! Initial time

	REAL(KIND=8), PARAMETER:: v0 = 5.93e+9;
	! minimal velocity

	REAL(KIND=8), PARAMETER:: v_min = 1.43e+9;
	! minimal velocity

        REAL(KIND=8), PARAMETER:: v_max = 2.1e+10;
	! Maximal velocity

        REAL(KIND=8), PARAMETER:: v_c = 3e+10;
	! Maximal velocity

	REAL(KIND=8), PARAMETER:: e = 4.8e-10;
        ! Electron charge
	
	REAL(KIND=8), PARAMETER:: v_p = 4.2e+7;
	! Plasma frequency for dimensionless parameters

	REAL(KIND=8), PARAMETER::  m_e = 9.1e-28;
	! electron mass
	
	REAL(KIND=8), PARAMETER:: m_i = 1.6726485e-24;
	! proton mass
	
        REAL(KIND=8), PARAMETER:: k_b = 1.380662e-16;
	! Boltzman constant
	
	REAL(KIND=8), PARAMETER:: G = 6.6720e-8;
	! gravitational constant
	
	REAL(KIND=8), PARAMETER:: mu = 0.6;
	! mean molecular value
	! Priest E.R.(1982) Solar Magnetohydrodynamics. Reidel,Dordrecht
       
	REAL(KIND=8), PARAMETER:: m_s = 1.99e+33;
	! Solar mass
	
	REAL(KIND=8), PARAMETER:: r_s = 6.958e+10 
	! Solar radius
	
	REAL(KIND=8), PARAMETER:: AU = 1.5d+13
	! Astronomical Unit	

	INTEGER, PARAMETER:: n_v = 50;
	! Number of V cells on one side

	INTEGER, PARAMETER:: n_k = 50;
	! Number of extra K cells on one side

	INTEGER, PARAMETER:: n_kv = 100;
	! Number of total cells on one side == n_v + n_k

	INTEGER, PARAMETER:: dbound = 51

        CHARACTER(*), PARAMETER:: vxfwd_sequence = 'fwd'
        CHARACTER(*), PARAMETER:: fxv_matrix = 'fxv'
        CHARACTER(*), PARAMETER:: wxv_matrix = 'wxv'
        CHARACTER(*), PARAMETER:: w_em_matrix = 'w_em'
        CHARACTER(*), PARAMETER:: w_em_f_matrix = 'w_em_f'
	CHARACTER(*), PARAMETER:: w_s_matrix = 'w_s'
	CHARACTER(*), PARAMETER:: corona_file = 'xdf.dat'
	CHARACTER(*), PARAMETER:: dynspecf_file = 'dynSpecF'
        CHARACTER(*), PARAMETER:: dynspech_file = 'dynSpecH'
        CHARACTER(*), PARAMETER:: params_file = 'params.dat'
	CHARACTER(*), PARAMETER:: velocity_file = 'velocity.dat'

	! OUTPUT File names

        CHARACTER(*), PARAMETER:: init_params = 'initL.par'

        ! Input file name 

        CHARACTER (8):: xvfwn, st; 

End  Module CONSTANTS
