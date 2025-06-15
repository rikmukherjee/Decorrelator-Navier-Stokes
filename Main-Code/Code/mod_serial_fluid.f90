!!             mod_serial_fluid.f90
!! This is the module containing the defination and commons for the 
!! serial code. As all the serial ffts requal arrays of same type 
!! we expect this should not be changed as we go from one machine  
!! to another.
			module mod_serial_fluid
	implicit none
	save
!! ----velocites  -------------
	real*8,allocatable,dimension(:,:,:,:) :: Omegavort,E,Aorth,Atrans,Resistspfix,Omegap,Ep
	
	real*8,allocatable,dimension(:,:,:) :: Vk1,Vk2,Vk3,VWk1,VWk2,& 
                                        VWk3,Wk1,Wk2,Wk3,tm1,tm2,tm3,&
                                        fk1,fk2,fk3,g_hat,Sed_torque
                                                 
                                                                  
!----------lagrangian trajectaries-----------------------------------
	real*8,allocatable,dimension(:,:)::Xp,Yp,Zp,Uxp,Uyp,Uzp,Vxp,Vyp,Vzp,Axp,Ayp,Azp,Qp,Rp,discrip, &
																			dx_uxp,dx_uyp,dx_uzp,dy_uxp,dy_uyp,dy_uzp,dz_uxp,dz_uyp,dz_uzp, &
																			preAxp,preAyp,preAzp,dAdtxp,dAdtyp,dAdtzp,&
                                      q1,q2,q3,q4,q1dot,q2dot,q3dot,q4dot,&
                                      Omega23,Omega31,Omega12,E23,E31,E12,&
                                      q1_prev,q2_prev,q3_prev,q4_prev,quaterNorm,costhetaEuler,thetaEuler,phiEuler,psiEuler,&
                                      omega_x,omega_y,omega_z,omega_x2,omega_y2,omega_z2,&
                                      omega_x3,omega_y3,omega_z3,omega_x4,omega_y4,omega_z4,Angx1,Angy1,Angz1,Angx2,Angy2,Angz2,&
                                      Angx3,Angy3,Angz3,Angx4,Angy4,Angz4,vort_mag

	real*8,allocatable,dimension(:)::tau,one_by_tau,Gr,one_by_taupRot,GrRot
  logical::forcing,particl
  real*8,allocatable,dimension(:,:,:)::dx_ux,dx_uy,dx_uz, &
                              dy_ux,dy_uy,dy_uz,dz_ux,dz_uy,dz_uz

	real*8,allocatable,dimension(:)::uxe,uye,uze,dx_uxe,dx_uye,dx_uze,dy_uxe,dy_uye,dy_uze,dz_uxe,dz_uye,dz_uze
	real*8,allocatable,dimension(:,:):: kq11,kq21,kq31,kq41
	real*8,allocatable,dimension(:,:):: kq12,kq22,kq32,kq42
	real*8,allocatable,dimension(:,:):: kq13,kq23,kq33,kq43
	real*8,allocatable,dimension(:,:):: kq14,kq24,kq34,kq44
	real*8,allocatable,dimension(:,:):: q12,q22,q32,q42
	real*8,allocatable,dimension(:,:):: q13,q23,q33,q43
	real*8,allocatable,dimension(:,:):: q14,q24,q34,q44
        real*8,allocatable,dimension(:,:):: tempsum1,tempsum2,tempVx1,tempVy1,tempVz1,tempVx2,tempVy2,tempVz2
        real*8,allocatable,dimension(:,:):: Resistprt,Mobprt

	integer,allocatable,dimension(:,:,:):: mask	
	integer,allocatable,dimension(:):: N_modes	
	real*8:: D
!! --energy,viscosity,forced modes etc --------------
	real*8,allocatable,dimension(:,:) :: E_Omega  
	double precision,allocatable,dimension(:) :: iniEk 
	integer,allocatable,dimension(:,:) :: forced_modes1,forced_modes2,forced_modes3
	integer :: iter,count1,count2,count3
!! ---parameters of simulations -----------------
	real*8 :: edr,tlr_micro_scale,dissipation_scale,vkrms,Rlambda, &
              tau_eddy,integral_scale,dx,dy,dz,t1,t2,t,dt_by_dx
!! -----real space structure functions -----------------------
	real*8,allocatable,dimension(:,:) :: S_p
	integer :: pmax,nthreads,count_out,cnt
!! -----common parameters ----------------------------
	integer :: thousand,nn,maxiter,nrun,n1,n2,n3,n1d,n1h,n1hf, &
               nshell,navg,nalias_sqr,nprtcl,nstokes,nrunpart
	real*8::r1,r2! axes-ratios of the ellipsoid to be read in from a file
	real*8 ::constA,constB,vis,vis2,delta,pi,length,factor,umax1,umax2,umax3,umax,c_cfl,Fr
	double precision :: fixed_energy1,fixed_energy2,fixed_energy3
	double precision :: Xa,Ya,onebyXa,lk,Xc,Yc,Yh,Re_lambda,Torque_sed,tauK,Uo
	complex*16 :: zi,zone 
!! -----global arrays --------------
	real*8,allocatable,dimension(:) :: time_increment
	integer,allocatable,dimension(:) :: den_state
!! ----------------------------------------------------------
!! -----------For FFTW --------------------------------------
	integer*8:: pfor,pinv
	integer,dimension(3) :: dim
	double precision :: scale
	include "fftw_f77.i"
!! --------------------------------------------------------
			end module mod_serial_fluid 
