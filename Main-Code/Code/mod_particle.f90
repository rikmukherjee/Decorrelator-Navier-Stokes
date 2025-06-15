	module mod_particle
	use mod_serial_fluid
	contains
	
	subroutine allocate_particle
	implicit none
	integer::i1

	allocate(tau(nstokes-1),one_by_tau(nstokes-1),Gr(nstokes-1),one_by_taupRot(nstokes-1),GrRot(nstokes-1))

	open(unit = 11,file = '../../fluid_parameters.in',status = 'old')
	read(11,*) Re_lambda
	read(11,*) lk
	read(11,*) tauK
	read(11,*) Uo
	close(11)

            Xa=4.0*(1-r2*r2)**(1.5)/(3.0*(r2*sqrt(1-r2*r2)+acos(r2)-2*r2*r2*acos(r2)))
            Ya=-8.0*(1-r2*r2)**(1.5)/(3.0*(r2*sqrt(1-r2*r2)+(-3+2*r2*r2)*acos(r2)))
            onebyXa=1/Xa
            Xc=-2.0*(1-r2*r2)**(1.5)/(3.0*(r2*sqrt(1-r2*r2)-acos(r2)))
            Yc=-2.0*sqrt(1-r2*r2)*(r2*r2*r2*r2-1)/(3.0*(r2*sqrt(1-r2*r2)+acos(r2)-2*r2*r2*acos(r2)))
            Yh=-2.0*(1-r2*r2)**(2.5)/(3.0*(r2*sqrt(1-r2*r2)+acos(r2)-2*r2*r2*acos(r2)))

	open(unit = 11,file = 'stokes_no.in',status = 'old')
        print*,nstokes

	do i1 = 1, nstokes-1
		read(11,*) tau(i1)
		one_by_tau(i1) = 1.0d0/(tauK*tau(i1))
         	one_by_taupRot(i1) = 20*one_by_tau(i1)/(3*(1+r2*r2)*Xa)
		Gr(i1)=Fr*Uo*one_by_tau(i1)
		GrRot(i1)=30*pi*Gr(i1)*Gr(i1)*Xa/(vis*(1+r2*r2)*one_by_tau(i1)) !Sedimentation torque to be read from flu.in
	enddo
	close(11)
	
    print*,Gr(1),GrRot(1),one_by_tau(1)
	print*,one_by_taupRot(1), Torque_sed

	allocate(Xp(nprtcl,nstokes),Yp(nprtcl,nstokes),Zp(nprtcl,nstokes))
	Xp = 0.0d0; Yp = 0.0d0; Zp = 0.0d0
	
	allocate(Uxp(nprtcl,nstokes),Uyp(nprtcl,nstokes),Uzp(nprtcl,nstokes))
	Uxp = 0.0d0; Uyp = 0.0d0; Uzp = 0.0d0 

	allocate(Vxp(nprtcl,nstokes),Vyp(nprtcl,nstokes),Vzp(nprtcl,nstokes))
	Vxp = 0.0d0; Vyp = 0.0d0; Vzp = 0.0d0 

	allocate(Axp(nprtcl,nstokes),Ayp(nprtcl,nstokes),Azp(nprtcl,nstokes))
	Axp = 0.0d0; Ayp = 0.0d0; Azp = 0.0d0 

	allocate(preAxp(nprtcl,nstokes),preAyp(nprtcl,nstokes),preAzp(nprtcl,nstokes))
	preAxp = 0.0d0; preAyp = 0.0d0; preAzp = 0.0d0 

	allocate(dAdtxp(nprtcl,nstokes),dAdtyp(nprtcl,nstokes),dAdtzp(nprtcl,nstokes))
	dAdtxp = 0.0d0; dAdtyp = 0.0d0; dAdtzp = 0.0d0 

	allocate(dx_uxp(nprtcl,nstokes),dx_uyp(nprtcl,nstokes),dx_uzp(nprtcl,nstokes))
	allocate(dy_uxp(nprtcl,nstokes),dy_uyp(nprtcl,nstokes),dy_uzp(nprtcl,nstokes))
	allocate(dz_uxp(nprtcl,nstokes),dz_uyp(nprtcl,nstokes),dz_uzp(nprtcl,nstokes))
	dx_uxp = 0.0d0; dx_uyp = 0.0d0; dx_uzp = 0.0d0 
	dy_uxp = 0.0d0; dy_uyp = 0.0d0; dy_uzp = 0.0d0 
	dz_uxp = 0.0d0; dz_uyp = 0.0d0; dz_uzp = 0.0d0 

	allocate(Qp(nprtcl,nstokes),Rp(nprtcl,nstokes),discrip(nprtcl,nstokes))
	Qp = 0.0d0; Rp = 0.0d0; discrip = 0.0d0 

	allocate(costhetaEuler(nprtcl,nstokes),thetaEuler(nprtcl,nstokes),phiEuler(nprtcl,nstokes),psiEuler(nprtcl,nstokes))
  	costhetaEuler = 0.0d0;thetaEuler = 0.0d0; phiEuler=0.0d0;psiEuler=0.0d0
	
    allocate(q1(nprtcl,nstokes),q2(nprtcl,nstokes),q3(nprtcl,nstokes),q4(nprtcl,nstokes),quaterNorm(nprtcl,nstokes),Resistprt(3,3),Mobprt(3,3))
  	q1 = 0.0d0;q2 = 0.0d0; q3=0.0d0;q4=0.0d0;quaterNorm=0.0d0;Resistprt=0.0d0;Mobprt=0.0d0
  	
  	allocate(g_hat(nprtcl,nstokes,3),Sed_torque(nprtcl,nstokes,3))
	g_hat = 0.0d0; Sed_torque=0.0d0;
  	
    allocate(q1dot(nprtcl,nstokes),q2dot(nprtcl,nstokes),q3dot(nprtcl,nstokes),q4dot(nprtcl,nstokes))
  	q1dot = 0.0d0; q2dot=0.0d0;q3dot=0.0d0;q4dot=0.0d0

    allocate(omega_x(nprtcl,nstokes),omega_y(nprtcl,nstokes),omega_z(nprtcl,nstokes))
  	omega_x = 0.0d0; omega_y=0.0d0;omega_z=0.0d0
  	
    allocate(omega_x2(nprtcl,nstokes),omega_y2(nprtcl,nstokes),omega_z2(nprtcl,nstokes))
  	omega_x2 = 0.0d0; omega_y2=0.0d0;omega_z2=0.0d0
  	
    allocate(omega_x3(nprtcl,nstokes),omega_y3(nprtcl,nstokes),omega_z3(nprtcl,nstokes))
  	omega_x3 = 0.0d0; omega_y3=0.0d0;omega_z3=0.0d0
  	
    allocate(omega_x4(nprtcl,nstokes),omega_y4(nprtcl,nstokes),omega_z4(nprtcl,nstokes))
  	omega_x4 = 0.0d0; omega_y4=0.0d0;omega_z4=0.0d0
 
	allocate(Omega23(nprtcl,nstokes),Omega31(nprtcl,nstokes),Omega12(nprtcl,nstokes),vort_mag(nprtcl,nstokes))
  	Omega23= 0.0d0; Omega31=0.0d0;Omega12=0.0d0;vort_mag=0.0d0

	allocate(E23(nprtcl,nstokes),E31(nprtcl,nstokes),E12(nprtcl,nstokes))
  	E23= 0.0d0; E31=0.0d0;E12=0.0d0

	allocate(Omegavort(nprtcl,nstokes,3,3),E(nprtcl,nstokes,3,3),Aorth(nprtcl,nstokes,3,3),Atrans(nprtcl,nstokes,3,3),Resistspfix(nprtcl,nstokes,3,3),Ep(nprtcl,nstokes,3,3),Omegap(nprtcl,nstokes,3,3))
  	Omegavort= 0.0d0;E=0.0d0;Aorth=0.0d0;Resistspfix=0.0d0;Atrans=0.0d0;Omegap=0.0d0;Ep=0.0d0        
        
	allocate(tempsum1(nprtcl,nstokes),tempsum2(nprtcl,nstokes))
        tempsum1=0.0d0;tempsum2=0.0d0;

  	allocate(q1_prev(nprtcl,nstokes),q2_prev(nprtcl,nstokes),q3_prev(nprtcl,nstokes),q4_prev(nprtcl,nstokes))
  	q1_prev = 0.0d0; q2_prev=0.0d0;q3_prev=0.0d0;q4_prev=0.0d0
  	
  	allocate(kq11(nprtcl,nstokes),kq21(nprtcl,nstokes),kq31(nprtcl,nstokes),kq41(nprtcl,nstokes))
  	allocate(kq12(nprtcl,nstokes),kq22(nprtcl,nstokes),kq32(nprtcl,nstokes),kq42(nprtcl,nstokes))
  	allocate(kq13(nprtcl,nstokes),kq23(nprtcl,nstokes),kq33(nprtcl,nstokes),kq43(nprtcl,nstokes))
  	allocate(kq14(nprtcl,nstokes),kq24(nprtcl,nstokes),kq34(nprtcl,nstokes),kq44(nprtcl,nstokes))
        kq11 = 0.0d0; kq21=0.0d0;kq31=0.0d0;kq41=0.0d0 
        kq12 = 0.0d0; kq22=0.0d0;kq32=0.0d0;kq42=0.0d0 
        kq13 = 0.0d0; kq23=0.0d0;kq33=0.0d0;kq43=0.0d0 
        kq14 = 0.0d0; kq24=0.0d0;kq34=0.0d0;kq44=0.0d0 

  	allocate(q12(nprtcl,nstokes),q22(nprtcl,nstokes),q32(nprtcl,nstokes),q42(nprtcl,nstokes))
  	allocate(q13(nprtcl,nstokes),q23(nprtcl,nstokes),q33(nprtcl,nstokes),q43(nprtcl,nstokes))
  	allocate(q14(nprtcl,nstokes),q24(nprtcl,nstokes),q34(nprtcl,nstokes),q44(nprtcl,nstokes))
        q12 = 0.0d0; q22=0.0d0;q32=0.0d0;q42=0.0d0 
        q13 = 0.0d0; q23=0.0d0;q33=0.0d0;q43=0.0d0  
        q14 = 0.0d0; q24=0.0d0;q34=0.0d0;q44=0.0d0 

	allocate(Angx1(nprtcl,nstokes),Angy1(nprtcl,nstokes),Angz1(nprtcl,nstokes))
	Angx1 = 0.0d0; Angy1 = 0.0d0; Angz1 = 0.0d0 

	allocate(Angx2(nprtcl,nstokes),Angy2(nprtcl,nstokes),Angz2(nprtcl,nstokes))
	Angx2 = 0.0d0; Angy2 = 0.0d0; Angz2 = 0.0d0 
		
        allocate(Angx3(nprtcl,nstokes),Angy3(nprtcl,nstokes),Angz3(nprtcl,nstokes))
	Angx3 = 0.0d0; Angy3 = 0.0d0; Angz3 = 0.0d0 
	
        allocate(Angx4(nprtcl,nstokes),Angy4(nprtcl,nstokes),Angz4(nprtcl,nstokes))
	Angx4 = 0.0d0; Angy4 = 0.0d0; Angz4 = 0.0d0 
	

	end subroutine allocate_particle
!------------------------------------------------------------------------------

	subroutine initialize_particle
	implicit none
	integer::i1,intx,inty,intz,ip
	if (nrunpart == 1) then

	        call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,Vk1, 0)
		call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,Vk2, 0)
		call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,Vk3, 0)
!! -------normalise etc ------------------------------
	  do i1=n1+1,n1+2
                        Vk1(i1,:,:) = 0.0d0
			Vk2(i1,:,:) = 0.0d0
			Vk3(i1,:,:) = 0.0d0
	  enddo
		Vk1 = Vk1*scale
        Vk2 = Vk2*scale
		Vk3 = Vk3*scale


	  call random_seed()
	  call random_number(Xp(:,1)); call random_number(Yp(:,1)); call random_number(Zp(:,1))

	  do i1 = 1,nprtcl
            intx = floor(Xp(i1,1)*dfloat(n1))+1
			inty = floor(Yp(i1,1)*dfloat(n2))+1
			intz = floor(Zp(i1,1)*dfloat(n3))+1
	    Xp(i1,:) = (intx-1)*dx; Yp(i1,:) = (inty-1)*dy; Zp(i1,:) = (intz-1)*dz
                        Vxp(i1,:)=Vk1(intx,inty,intz)
                        Vyp(i1,:)=Vk2(intx,inty,intz)
                        Vzp(i1,:)=Vk3(intx,inty,intz)
		enddo

        call rfftwnd_f77_threads_one_real_to_complex(nthreads,pfor,Vk1, 0)
		call rfftwnd_f77_threads_one_real_to_complex(nthreads,pfor,Vk2, 0)
		call rfftwnd_f77_threads_one_real_to_complex(nthreads,pfor,Vk3, 0)

	  call random_seed()
	  call random_number(costhetaEuler(:,1));call random_number(phiEuler(:,1)); call random_number(psiEuler(:,1))
	  
        do i1=1,nprtcl
	costhetaEuler(i1,1)=2*costhetaEuler(i1,1)-1
	thetaEuler(i1,1)=acos(costhetaEuler(i1,1))
	phiEuler(i1,1)=phiEuler(i1,1)*2*pi
	psiEuler(i1,1)=psiEuler(i1,1)*2*pi
	enddo

    do i1 = 1,nprtcl
    q1(i1,1) = SIN(0.5*thetaEuler(i1,1))*SIN(0.5*(psiEuler(i1,1)-phiEuler(i1,1)))
    q2(i1,1) = SIN(0.5*thetaEuler(i1,1))*COS(0.5*(psiEuler(i1,1)-phiEuler(i1,1)))
    q3(i1,1) = COS(0.5*thetaEuler(i1,1))*SIN(0.5*(psiEuler(i1,1)+phiEuler(i1,1)))
    q4(i1,1) = COS(0.5*thetaEuler(i1,1))*COS(0.5*(psiEuler(i1,1)+phiEuler(i1,1)))
    enddo
    
    do i1 = 2,nstokes
    do ip=1,nprtcl
    q1(ip,i1) = SIN(0.5*thetaEuler(ip,1))*SIN(0.5*(psiEuler(ip,1)-phiEuler(ip,1)))
    q2(ip,i1) = SIN(0.5*thetaEuler(ip,1))*COS(0.5*(psiEuler(ip,1)-phiEuler(ip,1)))
    q3(ip,i1) = COS(0.5*thetaEuler(ip,1))*SIN(0.5*(psiEuler(ip,1)+phiEuler(ip,1)))
    q4(ip,i1) = COS(0.5*thetaEuler(ip,1))*COS(0.5*(psiEuler(ip,1)+phiEuler(ip,1)))
    enddo
    enddo
        
        !Calculating the Resistance tensors(6*pi has already been included in tau_p)
    Resistprt(1,1)=Ya
    Resistprt(2,2)=Ya
    Resistprt(3,3)=Xa
    
    Mobprt(1,1)=1/(6*pi*Ya)
    Mobprt(2,2)=1/(6*pi*Ya)
    Mobprt(3,3)=1/(6*pi*Xa)
	endif
	end subroutine initialize_particle
!------------------------------------------------------------------------

	subroutine evolv_particle 
  implicit none
  integer::ist

  call get_particle_rhs
  
!The first set of particles is always Lagrangian

!$OMP WORKSHARE
  Xp=Xp+Vxp*delta
  Yp=Yp+Vyp*delta
  Zp=Zp+Vzp*delta
!$OMP END WORKSHARE

!Getting fluid velocities at the new particle positions
call get_velocity_at_particle_position

  do ist=2, nstokes
!$OMP WORKSHARE
    Vxp(:,ist)=Vxp(:,ist)+Axp(:,ist)*delta
    Vyp(:,ist)=Vyp(:,ist)+Ayp(:,ist)*delta
    Vzp(:,ist)=Vzp(:,ist)+Azp(:,ist)*delta
!$OMP END WORKSHARE
  enddo

  end subroutine evolv_particle
!---------------------------------------------------------------------------

	subroutine get_particle_rhs
	implicit none
	integer::ist,ip,ko,lo,io,jo

  call dermat_at_particle
   
   do ist=1,nstokes
   do ip=1,nprtcl
   !Backing up the quaternions in q_prev
    q1_prev(ip,ist) = q1(ip,ist)
    q2_prev(ip,ist) = q2(ip,ist)
    q3_prev(ip,ist) = q3(ip,ist)
    q4_prev(ip,ist) = q4(ip,ist)

    !Calculating the orthogonal transformation tensor-Aorth from space-fixed to body-fixed
   Aorth(ip,ist,1,1)=-q1(ip,ist)*q1(ip,ist)+q2(ip,ist)*q2(ip,ist)-q3(ip,ist)*q3(ip,ist)+q4(ip,ist)*q4(ip,ist)
   Aorth(ip,ist,2,2)=q1(ip,ist)*q1(ip,ist)-q2(ip,ist)*q2(ip,ist)-q3(ip,ist)*q3(ip,ist)+q4(ip,ist)*q4(ip,ist)
   Aorth(ip,ist,3,3)=-q1(ip,ist)*q1(ip,ist)-q2(ip,ist)*q2(ip,ist)+q3(ip,ist)*q3(ip,ist)+q4(ip,ist)*q4(ip,ist)
   Aorth(ip,ist,1,2)=-2*(q1(ip,ist)*q2(ip,ist)-q3(ip,ist)*q4(ip,ist))
   Aorth(ip,ist,2,1)=-2*(q1(ip,ist)*q2(ip,ist)+q3(ip,ist)*q4(ip,ist))
   Aorth(ip,ist,1,3)=2*(q2(ip,ist)*q3(ip,ist)+q1(ip,ist)*q4(ip,ist))
   Aorth(ip,ist,3,1)=2*(q2(ip,ist)*q3(ip,ist)-q1(ip,ist)*q4(ip,ist))
   Aorth(ip,ist,2,3)=-2*(q1(ip,ist)*q3(ip,ist)-q2(ip,ist)*q4(ip,ist))
   Aorth(ip,ist,3,2)=-2*(q1(ip,ist)*q3(ip,ist)+q2(ip,ist)*q4(ip,ist))
  do ko=1,3
    do lo=1,3
      Atrans(ip,ist,ko,lo)=Aorth(ip,ist,lo,ko)
    enddo
  enddo
  
  !Resistance tensor transformed to the space-fixed system to be used in the force balance equations
  do ko=1,3
    do lo=1,3
      Resistspfix(ip,ist,ko,lo)=Atrans(ip,ist,ko,1)*Atrans(ip,ist,lo,1)*Resistprt(1,1)+Atrans(ip,ist,ko,2)*Atrans(ip,ist,lo,2)*Resistprt(2,2)+Atrans(ip,ist,ko,3)*Atrans(ip,ist,lo,3)*Resistprt(3,3)
    enddo
  enddo

  Omega12(ip,ist)=0.5*(dx_uyp(ip,ist)-dy_uxp(ip,ist))
  E12(ip,ist)=0.5*(dx_uyp(ip,ist)+dy_uxp(ip,ist))
  Omega23(ip,ist)=0.5*(dy_uzp(ip,ist)-dz_uyp(ip,ist))
  E23(ip,ist)=0.5*(dy_uzp(ip,ist)+dz_uyp(ip,ist))
  Omega31(ip,ist)=0.5*(dz_uxp(ip,ist)-dx_uzp(ip,ist))
  E31(ip,ist)=0.5*(dz_uxp(ip,ist)+dx_uzp(ip,ist))

!Calculating vorticity magnitude
  vort_mag(ip,ist)=2*sqrt(Omega12(ip,ist)*Omega12(ip,ist)+Omega23(ip,ist)*Omega23(ip,ist)+Omega31(ip,ist)*Omega31(ip,ist))  
  
!Calculating the vorticity tensor-Omega 
      Omegavort(ip,ist,1,2)=Omega12(ip,ist)
      Omegavort(ip,ist,2,1)=-Omegavort(ip,ist,1,2)
      Omegavort(ip,ist,2,3)=Omega23(ip,ist)
      Omegavort(ip,ist,3,2)=-Omegavort(ip,ist,2,3)
      Omegavort(ip,ist,3,1)=Omega31(ip,ist)
      Omegavort(ip,ist,1,3)=-Omegavort(ip,ist,3,1)
      
!Calculating the rate of strain tensor-E
      E(ip,ist,1,2)=E12(ip,ist)
      E(ip,ist,2,1)=E(ip,ist,1,2)
      E(ip,ist,2,3)=E23(ip,ist)
      E(ip,ist,3,2)=E(ip,ist,2,3)
      E(ip,ist,3,1)=E31(ip,ist)
      E(ip,ist,1,3)=E(ip,ist,3,1)
      E(ip,ist,1,1)=dx_uxp(ip,ist)
      E(ip,ist,2,2)=dy_uyp(ip,ist)
      E(ip,ist,3,3)=dz_uzp(ip,ist)

! Transforming rate-of-strain and vorticity tensor from space-fixed to body-fixed coordinates
 do io=1,3
    do jo=1,3
    tempsum1(ip,ist)=0.0d0
    tempsum2(ip,ist)=0.0d0
        do ko=1,3
            do lo=1,3
                tempsum1(ip,ist)=tempsum1(ip,ist)+Aorth(ip,ist,io,ko)*Aorth(ip,ist,jo,lo)*Omegavort(ip,ist,ko,lo)
                tempsum2(ip,ist)=tempsum2(ip,ist)+Aorth(ip,ist,io,ko)*Aorth(ip,ist,jo,lo)*E(ip,ist,ko,lo)
            enddo
        enddo
    Omegap(ip,ist,io,jo)=tempsum1(ip,ist)
    Ep(ip,ist,io,jo)=tempsum2(ip,ist)
    enddo
 enddo
enddo
enddo

!Updating angular velocity for tracer particles
   do ip=1,nprtcl
  omega_x(ip,1)=Omegap(ip,1,2,3)+(r1*r1-r2*r2)/(r2*r2+r1*r1)*Ep(ip,1,2,3)
  omega_y(ip,1)=Omegap(ip,1,3,1)+(r2*r2-1)/(r2*r2+1)*Ep(ip,1,3,1)
  omega_z(ip,1)=Omegap(ip,1,1,2)+(1-r1*r1)/(1+r1*r1)*Ep(ip,1,1,2)  
   enddo

   !Initializing angular velocity of inertial particles to be the same as that of tracer particles
if(cnt==1) then
  do ist=2,nstokes
    do ip=1,nprtcl
    omega_x(ip,ist)=Omegap(ip,ist,2,3)+(r1*r1-r2*r2)/(r2*r2+r1*r1)*Ep(ip,ist,2,3)
    omega_y(ip,ist)=Omegap(ip,ist,3,1)+(r2*r2-1)/(r2*r2+1)*Ep(ip,ist,3,1)
    omega_z(ip,ist)=Omegap(ip,ist,1,2)+(1-r1*r1)/(1+r1*r1)*Ep(ip,ist,1,2)
    enddo
  enddo
endif

!Solving the kinematic relation dq/dt
do ist=1,nstokes
   do ip=1,nprtcl
  kq11(ip,ist)=0.5*(-q3(ip,ist)*omega_x(ip,ist)-q4(ip,ist)*omega_y(ip,ist)+q2(ip,ist)*omega_z(ip,ist))
  kq21(ip,ist)=0.5*(q4(ip,ist)*omega_x(ip,ist)-q3(ip,ist)*omega_y(ip,ist)-q1(ip,ist)*omega_z(ip,ist))
  kq31(ip,ist)=0.5*(q1(ip,ist)*omega_x(ip,ist)+q2(ip,ist)*omega_y(ip,ist)+q4(ip,ist)*omega_z(ip,ist))
  kq41(ip,ist)=0.5*(-q2(ip,ist)*omega_x(ip,ist)+q1(ip,ist)*omega_y(ip,ist)-q3(ip,ist)*omega_z(ip,ist))
  
          q12(ip,ist)=q1(ip,ist)+kq11(ip,ist)*delta*0.5
          q22(ip,ist)=q2(ip,ist)+kq21(ip,ist)*delta*0.5
          q32(ip,ist)=q3(ip,ist)+kq31(ip,ist)*delta*0.5
          q42(ip,ist)=q4(ip,ist)+kq41(ip,ist)*delta*0.5
          
  kq12(ip,ist)=0.5*(-q32(ip,ist)*omega_x(ip,ist)-q42(ip,ist)*omega_y(ip,ist)+q22(ip,ist)*omega_z(ip,ist))
  kq22(ip,ist)=0.5*(q42(ip,ist)*omega_x(ip,ist)-q32(ip,ist)*omega_y(ip,ist)-q12(ip,ist)*omega_z(ip,ist))
  kq32(ip,ist)=0.5*(q12(ip,ist)*omega_x(ip,ist)+q22(ip,ist)*omega_y(ip,ist)+q42(ip,ist)*omega_z(ip,ist))
  kq42(ip,ist)=0.5*(-q22(ip,ist)*omega_x(ip,ist)+q12(ip,ist)*omega_y(ip,ist)-q32(ip,ist)*omega_z(ip,ist))
 
          q13(ip,ist)=q1(ip,ist)+kq12(ip,ist)*delta*0.5
          q23(ip,ist)=q2(ip,ist)+kq22(ip,ist)*delta*0.5
          q33(ip,ist)=q3(ip,ist)+kq32(ip,ist)*delta*0.5
          q43(ip,ist)=q4(ip,ist)+kq42(ip,ist)*delta*0.5
             
  kq13(ip,ist)=0.5*(-q33(ip,ist)*omega_x(ip,ist)-q43(ip,ist)*omega_y(ip,ist)+q23(ip,ist)*omega_z(ip,ist))
  kq23(ip,ist)=0.5*(q43(ip,ist)*omega_x(ip,ist)-q33(ip,ist)*omega_y(ip,ist)-q13(ip,ist)*omega_z(ip,ist))
  kq33(ip,ist)=0.5*(q13(ip,ist)*omega_x(ip,ist)+q23(ip,ist)*omega_y(ip,ist)+q43(ip,ist)*omega_z(ip,ist))
  kq43(ip,ist)=0.5*(-q23(ip,ist)*omega_x(ip,ist)+q13(ip,ist)*omega_y(ip,ist)-q33(ip,ist)*omega_z(ip,ist))

          q14(ip,ist)=q1(ip,ist)+kq13(ip,ist)*delta
          q24(ip,ist)=q2(ip,ist)+kq23(ip,ist)*delta
          q34(ip,ist)=q3(ip,ist)+kq33(ip,ist)*delta
          q44(ip,ist)=q4(ip,ist)+kq43(ip,ist)*delta
            
  kq14(ip,ist)=0.5*(-q34(ip,ist)*omega_x(ip,ist)-q44(ip,ist)*omega_y(ip,ist)+q24(ip,ist)*omega_z(ip,ist))
  kq24(ip,ist)=0.5*(q44(ip,ist)*omega_x(ip,ist)-q34(ip,ist)*omega_y(ip,ist)-q14(ip,ist)*omega_z(ip,ist))
  kq34(ip,ist)=0.5*(q14(ip,ist)*omega_x(ip,ist)+q24(ip,ist)*omega_y(ip,ist)+q44(ip,ist)*omega_z(ip,ist))
  kq44(ip,ist)=0.5*(-q24(ip,ist)*omega_x(ip,ist)+q14(ip,ist)*omega_y(ip,ist)-q34(ip,ist)*omega_z(ip,ist))
            
    ! Updating the values at the end of time-step 
      q1(ip,ist)=q1(ip,ist)+(kq11(ip,ist)+2*kq12(ip,ist)+2*kq13(ip,ist)+kq14(ip,ist))*delta/6.0;
      q2(ip,ist)=q2(ip,ist)+(kq21(ip,ist)+2*kq22(ip,ist)+2*kq23(ip,ist)+kq24(ip,ist))*delta/6.0;
      q3(ip,ist)=q3(ip,ist)+(kq31(ip,ist)+2*kq32(ip,ist)+2*kq33(ip,ist)+kq34(ip,ist))*delta/6.0;
      q4(ip,ist)=q4(ip,ist)+(kq41(ip,ist)+2*kq42(ip,ist)+2*kq43(ip,ist)+kq44(ip,ist))*delta/6.0;
      
      !Finding the quaternion norm at the end of the time step
      quaterNorm(ip,ist)=sqrt(q1(ip,ist)*q1(ip,ist)+q2(ip,ist)*q2(ip,ist)+q3(ip,ist)*q3(ip,ist)+q4(ip,ist)*q4(ip,ist))
      !Renormalizing
      q1(ip,ist)=q1(ip,ist)/quaterNorm(ip,ist)
      q2(ip,ist)=q2(ip,ist)/quaterNorm(ip,ist)
      q3(ip,ist)=q3(ip,ist)/quaterNorm(ip,ist)
      q4(ip,ist)=q4(ip,ist)/quaterNorm(ip,ist)
      
   enddo
enddo

!Solving the dynamic equation dw/dt for ist>=2 using the quaternions at the beginning of the time-step
do ist=2,nstokes
   do ip=1,nprtcl
     
    g_hat(ip,ist,1)=-2.0*(q2_prev(ip,ist)*q3_prev(ip,ist)+q1_prev(ip,ist)*q4_prev(ip,ist))
    g_hat(ip,ist,2)=2.0*(q1_prev(ip,ist)*q3_prev(ip,ist)-q2_prev(ip,ist)*q4_prev(ip,ist))
    g_hat(ip,ist,3)=q1_prev(ip,ist)*q1_prev(ip,ist)+q2_prev(ip,ist)*q2_prev(ip,ist)-q3_prev(ip,ist)*q3_prev(ip,ist)-q4_prev(ip,ist)*q4_prev(ip,ist)
    
    Sed_torque(ip,ist,1)=GrRot(ist-1)*Mobprt(2,2)*Mobprt(3,3)*g_hat(ip,ist,2)*g_hat(ip,ist,3)*Torque_sed
    Sed_torque(ip,ist,2)=-GrRot(ist-1)*Mobprt(1,1)*Mobprt(3,3)*g_hat(ip,ist,1)*g_hat(ip,ist,3)*Torque_sed
    
    !RK4 details
    Angx1(ip,ist)=omega_y(ip,ist)*omega_z(ip,ist)*(r2*r2-1)/(r2*r2+1)+Sed_torque(ip,ist,1)+one_by_taupRot(ist-1)*(Yc*(Omegap(ip,ist,2,3)-omega_x(ip,ist))-Yh*Ep(ip,ist,2,3))
    Angy1(ip,ist)=omega_x(ip,ist)*omega_z(ip,ist)*(1-r2*r2)/(r2*r2+1)+Sed_torque(ip,ist,2)+one_by_taupRot(ist-1)*(Yc*(Omegap(ip,ist,3,1)-omega_y(ip,ist))+Yh*Ep(ip,ist,1,3))
    Angz1(ip,ist)=0.5*one_by_taupRot(ist-1)*(1+r2*r2)*Xc*(Omegap(ip,ist,1,2)-omega_z(ip,ist))    
    
    omega_x2(ip,ist)=omega_x(ip,ist)+Angx1(ip,ist)*delta*0.5
    omega_y2(ip,ist)=omega_y(ip,ist)+Angy1(ip,ist)*delta*0.5
    omega_z2(ip,ist)=omega_z(ip,ist)+Angz1(ip,ist)*delta*0.5
    
    Angx2(ip,ist)=omega_y2(ip,ist)*omega_z2(ip,ist)*(r2*r2-1)/(r2*r2+1)+Sed_torque(ip,ist,1)+one_by_taupRot(ist-1)*(Yc*(Omegap(ip,ist,2,3)-omega_x2(ip,ist))-Yh*Ep(ip,ist,2,3))
    Angy2(ip,ist)=omega_x2(ip,ist)*omega_z2(ip,ist)*(1-r2*r2)/(r2*r2+1)+Sed_torque(ip,ist,2)+one_by_taupRot(ist-1)*(Yc*(Omegap(ip,ist,3,1)-omega_y2(ip,ist))+Yh*Ep(ip,ist,1,3))
    Angz2(ip,ist)=0.5*one_by_taupRot(ist-1)*(1+r2*r2)*Xc*(Omegap(ip,ist,1,2)-omega_z2(ip,ist))
    
    omega_x3(ip,ist)=omega_x(ip,ist)+Angx2(ip,ist)*delta*0.5
    omega_y3(ip,ist)=omega_y(ip,ist)+Angy2(ip,ist)*delta*0.5
    omega_z3(ip,ist)=omega_z(ip,ist)+Angz2(ip,ist)*delta*0.5
    
    Angx3(ip,ist)=omega_y3(ip,ist)*omega_z3(ip,ist)*(r2*r2-1)/(r2*r2+1)+Sed_torque(ip,ist,1)+one_by_taupRot(ist-1)*(Yc*(Omegap(ip,ist,2,3)-omega_x3(ip,ist))-Yh*Ep(ip,ist,2,3))
    Angy3(ip,ist)=omega_x3(ip,ist)*omega_z3(ip,ist)*(1-r2*r2)/(r2*r2+1)+Sed_torque(ip,ist,2)+one_by_taupRot(ist-1)*(Yc*(Omegap(ip,ist,3,1)-omega_y3(ip,ist))+Yh*Ep(ip,ist,1,3))
    Angz3(ip,ist)=0.5*one_by_taupRot(ist-1)*(1+r2*r2)*Xc*(Omegap(ip,ist,1,2)-omega_z3(ip,ist))
    
    omega_x4(ip,ist)=omega_x(ip,ist)+Angx3(ip,ist)*delta
    omega_y4(ip,ist)=omega_y(ip,ist)+Angy3(ip,ist)*delta
    omega_z4(ip,ist)=omega_z(ip,ist)+Angz3(ip,ist)*delta
    
    Angx4(ip,ist)=omega_y4(ip,ist)*omega_z4(ip,ist)*(r2*r2-1)/(r2*r2+1)+Sed_torque(ip,ist,1)+one_by_taupRot(ist-1)*(Yc*(Omegap(ip,ist,2,3)-omega_x4(ip,ist))-Yh*Ep(ip,ist,2,3))
    Angy4(ip,ist)=omega_x4(ip,ist)*omega_z4(ip,ist)*(1-r2*r2)/(r2*r2+1)+Sed_torque(ip,ist,2)+one_by_taupRot(ist-1)*(Yc*(Omegap(ip,ist,3,1)-omega_y4(ip,ist))+Yh*Ep(ip,ist,1,3))
    Angz4(ip,ist)=0.5*one_by_taupRot(ist-1)*(1+r2*r2)*Xc*(Omegap(ip,ist,1,2)-omega_z4(ip,ist))
    
!   Updating angular velocities at the end of the time-step
    omega_x(ip,ist)=omega_x(ip,ist)+(Angx1(ip,ist)+2*Angx2(ip,ist)+2*Angx3(ip,ist)+Angx4(ip,ist))*delta/6.0
    omega_y(ip,ist)=omega_y(ip,ist)+(Angy1(ip,ist)+2*Angy2(ip,ist)+2*Angy3(ip,ist)+Angy4(ip,ist))*delta/6.0
    omega_z(ip,ist)=omega_z(ip,ist)+(Angz1(ip,ist)+2*Angz2(ip,ist)+2*Angz3(ip,ist)+Angz4(ip,ist))*delta/6.0
   enddo
 enddo
 
  !Calculating qdot at the end of the time-step
do ist=1,nstokes
   do ip=1,nprtcl
  q1dot(ip,ist)=0.5*(-q3(ip,ist)*omega_x(ip,ist)-q4(ip,ist)*omega_y(ip,ist)+q2(ip,ist)*omega_z(ip,ist))
  q2dot(ip,ist)=0.5*(q4(ip,ist)*omega_x(ip,ist)-q3(ip,ist)*omega_y(ip,ist)-q1(ip,ist)*omega_z(ip,ist))
  q3dot(ip,ist)=0.5*(q1(ip,ist)*omega_x(ip,ist)+q2(ip,ist)*omega_y(ip,ist)+q4(ip,ist)*omega_z(ip,ist))
  q4dot(ip,ist)=0.5*(-q2(ip,ist)*omega_x(ip,ist)+q1(ip,ist)*omega_y(ip,ist)-q3(ip,ist)*omega_z(ip,ist))
   enddo
enddo

	call get_velocity_at_particle_position

!$OMP WORKSHARE
  Vxp(:,1) = Uxp(:,1)
  Vyp(:,1) = Uyp(:,1)
  Vzp(:,1) = Uzp(:,1)
!$OMP END WORKSHARE

  do ist=2,nstokes
  do ip=1,nprtcl

   Axp(ip,ist) = one_by_tau(ist-1)*onebyXa*(Resistspfix(ip,ist,1,1)*(Uxp(ip,ist)-Vxp(ip,ist))+Resistspfix(ip,ist,1,2)*(Uyp(ip,ist)-Vyp(ip,ist))+Resistspfix(ip,ist,1,3)*(Uzp(ip,ist)-Vzp(ip,ist)))
   
   Ayp(ip,ist) = one_by_tau(ist-1)*onebyXa*(Resistspfix(ip,ist,2,1)*(Uxp(ip,ist)-Vxp(ip,ist))+Resistspfix(ip,ist,2,2)*(Uyp(ip,ist)-Vyp(ip,ist))+Resistspfix(ip,ist,2,3)*(Uzp(ip,ist)-Vzp(ip,ist)))
   
   Azp(ip,ist) = -Gr(ist-1)+one_by_tau(ist-1)*onebyXa*(Resistspfix(ip,ist,3,1)*(Uxp(ip,ist)-Vxp(ip,ist))+Resistspfix(ip,ist,3,2)*(Uyp(ip,ist)-Vyp(ip,ist))+Resistspfix(ip,ist,3,3)*(Uzp(ip,ist)-Vzp(ip,ist)))

  enddo
  enddo

	end subroutine get_particle_rhs
!--------------------------------------------------------------------------
  subroutine get_dAdt
  implicit none
  integer::ist

  do ist=2,nstokes
!$OMP WORKSHARE
    dAdtxp(:,ist) = (Axp(:,ist)-preAxp(:,ist))/delta
    dAdtyp(:,ist) = (Ayp(:,ist)-preAyp(:,ist))/delta
    dAdtzp(:,ist) = (Azp(:,ist)-preAzp(:,ist))/delta
!$OMP END WORKSHARE
  enddo

  do ist=2,nstokes
!$OMP WORKSHARE
		preAxp(:,ist) = Axp(:,ist)
		preAyp(:,ist) = Ayp(:,ist)
		preAzp(:,ist) = Azp(:,ist)
!$OMP END WORKSHARE
  enddo

  end subroutine get_dAdt 
!--------------------------------------------------------------
	subroutine get_velocity_at_particle_position
	implicit none
	integer::ist,i1
	real*8::xxp,yyp,zzp,uuxp,uuyp,uuzp

  do ist = 1,nstokes
!$OMP PARALLEL DO PRIVATE(i1,xxp,yyp,zzp,uuxp,uuyp,uuzp) SHARED(Xp,Yp,Zp,Vk1,Vk2,Vk3,Uxp,Uyp,Uzp,n1,n2,n3,ist)
    do i1 = 1,nprtcl
      xxp = Xp(i1,ist); yyp = Yp(i1,ist); zzp = Zp(i1,ist)

      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,Vk1,uuxp)
      Uxp(i1,ist) = uuxp
      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,Vk2,uuyp)
      Uyp(i1,ist) = uuyp
      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,Vk3,uuzp)
      Uzp(i1,ist) = uuzp

    enddo
!$OMP END PARALLEL DO
  enddo

	end subroutine get_velocity_at_particle_position
!--------------------------------------------------------------------------

	subroutine linear_interp_to_offgrid(xpar,ypar,zpar,nn1,nn2,nn3,psi,psi_interpolated)
	implicit none
	real*8::xpar,ypar,zpar
	integer::nn1,nn2,nn3,xgrid,ygrid,zgrid,xgrid_nxt,ygrid_nxt,zgrid_nxt
	real*8,dimension(nn1+2,nn2,nn3)::psi
	real*8::psi_interpolated,a1,a2,b1,b2,c1,c2,xd,yd,zd

  call get_grid_loc(xpar,xgrid,xgrid_nxt,xd)
  call get_grid_loc(ypar,ygrid,ygrid_nxt,yd)
  call get_grid_loc(zpar,zgrid,zgrid_nxt,zd)
	
	a1=(1.0d0-zd)*psi(xgrid,ygrid,zgrid)+zd*psi(xgrid,ygrid,zgrid_nxt)
	a2=(1.0d0-zd)*psi(xgrid,ygrid_nxt,zgrid)+zd*psi(xgrid,ygrid_nxt,zgrid_nxt)
	b1=(1.0d0-zd)*psi(xgrid_nxt,ygrid,zgrid)+zd*psi(xgrid_nxt,ygrid,zgrid_nxt)
  b2=(1.0d0-zd)*psi(xgrid_nxt,ygrid_nxt,zgrid)+zd*psi(xgrid_nxt,ygrid_nxt,zgrid_nxt)

  c1=(1.0d0-yd)*a1+yd*a2
  c2=(1.0d0-yd)*b1+yd*b2

  psi_interpolated=(1.0d0-xd)*c1+xd*c2

	end subroutine linear_interp_to_offgrid
!--------------------------------------------------------------------------------------------------------------------------

	subroutine get_grid_loc(xpar,grid,grid_nxt,xd)
	implicit none
	real*8::xpar,xd
	integer::grid,grid_nxt
	real*8::Xp_fld

	Xp_fld=modulo(xpar,(2.0d0*pi))	
	grid=floor(Xp_fld/dx)+1
	grid_nxt=(1-grid/n1)*grid+1
	xd=modulo(Xp_fld,dx)/dx

	end subroutine get_grid_loc
!---------------------------------------------------------------------------------------------------------------------

	subroutine dermat_at_particle
	implicit none
	integer::ist,i1
	real*8::interpolated_value,xxp,yyp,zzp
	real*8,dimension(3,3)::DM,DM_sqr,DM_cub
	
  do ist = 1,nstokes
!$OMP PARALLEL DO PRIVATE(i1,xxp,yyp,zzp,interpolated_value,DM,DM_sqr,DM_cub) &
!$OMP SHARED(Xp,Yp,Zp,dx_ux,dx_uy,dx_uz,dy_ux,dy_uy,dy_uz,dz_ux,dz_uy,dz_uz,dx_uxp,dx_uyp,dx_uzp,dy_uxp,dy_uyp,dy_uzp,dz_uxp,dz_uyp,dz_uzp,n1,n2,n3,ist,Qp,Rp,discrip)

    do i1 = 1,nprtcl
      xxp = Xp(i1,ist); yyp = Yp(i1,ist); zzp = Zp(i1,ist)

      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,dx_ux,interpolated_value)
      dx_uxp(i1,ist) = interpolated_value
			DM(1,1) = interpolated_value
      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,dx_uy,interpolated_value)
      dx_uyp(i1,ist) = interpolated_value
			DM(1,2) = interpolated_value
      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,dx_uz,interpolated_value)
      dx_uzp(i1,ist) = interpolated_value
			DM(1,3) = interpolated_value

      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,dy_ux,interpolated_value)
      dy_uxp(i1,ist) = interpolated_value
			DM(2,1) = interpolated_value
      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,dy_uy,interpolated_value)
      dy_uyp(i1,ist) = interpolated_value
			DM(2,2) = interpolated_value
      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,dy_uz,interpolated_value)
      dy_uzp(i1,ist) = interpolated_value
			DM(2,3) = interpolated_value

      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,dz_ux,interpolated_value)
      dz_uxp(i1,ist) = interpolated_value
			DM(3,1) = interpolated_value
      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,dz_uy,interpolated_value)
      dz_uyp(i1,ist) = interpolated_value
			DM(3,2) = interpolated_value
      call linear_interp_to_offgrid(xxp,yyp,zzp,n1,n2,n3,dz_uz,interpolated_value)
      dz_uzp(i1,ist) = interpolated_value
			DM(3,3) = interpolated_value
			
			DM_sqr = matmul(DM,DM)
			DM_cub = matmul(DM,DM_sqr)
			
			Qp(i1,ist) = -0.5d0*(DM_sqr(1,1)+DM_sqr(2,2)+DM_sqr(3,3))
			Rp(i1,ist) = -(1.0d0/3.0d0)*(DM_cub(1,1)+DM_cub(2,2)+DM_cub(3,3))
			discrip(i1,ist) = Qp(i1,ist)**3+(27.0d0/4.0d0)*Rp(i1,ist)**2

    enddo
!$OMP END PARALLEL DO
  enddo
	
	end subroutine dermat_at_particle
!-----------------------------------------------------------------------------------------------------------

	subroutine open_output_particle
	implicit none
  character*80 :: fname1,fname2
  integer::i1,fno,ist
 !Printing tracer details separately
    ist=1
    write(fname1,'(g8.0)') ist
    fno = 1000*ist + 30
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/q1.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/q2.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/q3.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/q4.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/q1dot.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/q2dot.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/q3dot.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/q4dot.out',status='unknown')
    fno=fno+1
  
 !Printing  
  do ist = 2,nstokes
    write(fname1,'(g8.0)')  ist
    fno = 1000*ist + 30
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/xptrack.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/yptrack.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/zptrack.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/vxint.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/vyint.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/vzint.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/uxint.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/uyint.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/uzint.out',status='unknown')
    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dAdtxint.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dAdtyint.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dAdtzint.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/R_int.out',status='unknown')
!    fno = fno + 1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/Delta_int.out',status='unknown')
!    fno = fno + 1
!     open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dx_uxp.out',status='unknown')
!     fno = fno + 1
!     open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dx_uyp.out',status='unknown')
!     fno = fno + 1
!     open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dx_uzp.out',status='unknown')
!     fno = fno + 1
!     open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dy_uxp.out',status='unknown')
!     fno = fno + 1
!     open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dy_uyp.out',status='unknown')
!     fno = fno + 1
!     open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dy_uzp.out',status='unknown')
!     fno = fno + 1
!     open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dz_uxp.out',status='unknown')
!     fno = fno + 1
!     open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dz_uyp.out',status='unknown')
!     fno = fno + 1
!     open(unit=fno,file='tau'//trim(adjustl(fname1))//'/dz_uzp.out',status='unknown')
!     fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/q1.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/q2.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/q3.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/q4.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/q1dot.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/q2dot.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/q3dot.out',status='unknown')
    fno = fno + 1
    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/q4dot.out',status='unknown')
    fno=fno+1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/Egg.out',status='unknown')
!    fno=fno+1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/Qp.out',status='unknown')
!    fno=fno+1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/Rp.out',status='unknown')
!    fno=fno+1
!    open(unit=fno,file='tau'//trim(adjustl(fname1))//'/vort_mag.out',status='unknown')
!    fno=fno+1
  enddo
!!
	end subroutine open_output_particle
!-------------------------------------------------------------------------------------------------

	subroutine write_output_particle
	implicit none
	integer::i1,fno,ist
	character*500 :: cnpar,formp

      write(cnpar,'(g8.0)') nprtcl
      formp = '('//trim(adjustl(cnpar))//'es16.3E2)'
	  !Printing tracer details separately
      ist=1
      fno = 1000*ist + 30
	  write(fno,formp) (q1(i1,ist),i1=1,nprtcl); 
	  fno = fno + 1
      write(fno,formp) (q2(i1,ist),i1=1,nprtcl); 
      fno = fno + 1
      write(fno,formp) (q3(i1,ist),i1=1,nprtcl); 
      fno = fno + 1
      write(fno,formp) (q4(i1,ist),i1=1,nprtcl); 
      fno = fno + 1
      write(fno,formp) (q1dot(i1,ist),i1=1,nprtcl); 
      fno = fno + 1
      write(fno,formp) (q2dot(i1,ist),i1=1,nprtcl); 
      fno = fno + 1
      write(fno,formp) (q3dot(i1,ist),i1=1,nprtcl); 
      fno = fno + 1
      write(fno,formp) (q4dot(i1,ist),i1=1,nprtcl); 
      fno=fno+1
	
  do ist = 2,nstokes
    fno = 1000*ist + 30
    write(fno,formp) (Xp(i1,ist),i1=1,nprtcl);  
    fno = fno + 1
    write(fno,formp) (Yp(i1,ist),i1=1,nprtcl);  
    fno = fno + 1
    write(fno,formp) (Zp(i1,ist),i1=1,nprtcl);  
    fno = fno + 1
     write(fno,formp) (Vxp(i1,ist),i1=1,nprtcl); 
     fno = fno + 1
     write(fno,formp) (Vyp(i1,ist),i1=1,nprtcl); 
     fno = fno + 1
    write(fno,formp) (Vzp(i1,ist),i1=1,nprtcl); 
    fno = fno + 1
    write(fno,formp) (Uxp(i1,ist),i1=1,nprtcl); fno = fno + 1
    write(fno,formp) (Uyp(i1,ist),i1=1,nprtcl); fno = fno + 1
    write(fno,formp) (Uzp(i1,ist),i1=1,nprtcl); 
    fno = fno + 1
!    write(fno,formp) (dAdtxp(i1,ist),i1=1,nprtcl); fno = fno + 1
!    write(fno,formp) (dAdtyp(i1,ist),i1=1,nprtcl); fno = fno + 1
!    write(fno,formp) (dAdtzp(i1,ist),i1=1,nprtcl); fno = fno + 1
!    write(fno,formp) (discrip(i1,ist),i1=1,nprtcl); fno = fno + 1
!     write(fno,formp) (dx_uxp(i1,ist),i1=1,nprtcl); fno = fno + 1
!     write(fno,formp) (dx_uyp(i1,ist),i1=1,nprtcl); fno = fno + 1
!     write(fno,formp) (dx_uzp(i1,ist),i1=1,nprtcl); fno = fno + 1
!     write(fno,formp) (dy_uxp(i1,ist),i1=1,nprtcl); fno = fno + 1
!     write(fno,formp) (dy_uyp(i1,ist),i1=1,nprtcl); fno = fno + 1
!     write(fno,formp) (dy_uzp(i1,ist),i1=1,nprtcl); fno = fno + 1
!     write(fno,formp) (dz_uxp(i1,ist),i1=1,nprtcl); fno = fno + 1
!     write(fno,formp) (dz_uyp(i1,ist),i1=1,nprtcl); fno = fno + 1
!     write(fno,formp) (dz_uzp(i1,ist),i1=1,nprtcl); fno = fno + 1
      write(fno,formp) (q1(i1,ist),i1=1,nprtcl); 
      fno = fno + 1
      write(fno,formp) (q2(i1,ist),i1=1,nprtcl); 
      fno = fno + 1
      write(fno,formp) (q3(i1,ist),i1=1,nprtcl); 
      fno = fno + 1
      write(fno,formp) (q4(i1,ist),i1=1,nprtcl); 
      fno = fno + 1
      write(fno,formp) (q1dot(i1,ist),i1=1,nprtcl); 
      fno = fno + 1
      write(fno,formp) (q2dot(i1,ist),i1=1,nprtcl); 
      fno = fno + 1
      write(fno,formp) (q3dot(i1,ist),i1=1,nprtcl); 
      fno = fno + 1
      write(fno,formp) (q4dot(i1,ist),i1=1,nprtcl); 
      fno=fno+1
!      write(fno,formp) (E(i1,ist,3,3),i1=1,nprtcl);  !E:gg 
!      fno=fno+1
!      write(fno,formp) (Qp(i1,ist),i1=1,nprtcl); 
!      fno=fno+1
!      write(fno,formp) (Rp(i1,ist),i1=1,nprtcl); 
!      fno=fno+1
!      write(fno,formp) (vort_mag(i1,ist),i1=1,nprtcl); 
!      fno=fno+1
  enddo

		
	end subroutine write_output_particle

	end module mod_particle
