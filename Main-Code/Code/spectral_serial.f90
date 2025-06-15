!! This is a code for Direct Numerical Simulation (DNS) of Navier-Stokes
!! equation, in a periodic box. We simulate a divergence-less flow. 
!! This is a spectral code, i.e. the derivative are taken in Fourier Space. 
!!This is a 2/3rd de-aliased code ,this is used to prevent the blow up
!!of the tail b'cz of not resolving the largeK(or small l) modes(length scales).!!The other technique to avoid this is to add a viscousity(vis2*k^2) term to 
!! the viscousity.
!! ---------------------------------------------------------------------- 
	program spectral_serial
	use mod_serial_fluid
	use mod_initial
	use	omplib
	use mod_particle
	implicit none
!! other defination of variables 
	integer :: ispectra,i1,i2,i3,iouter,inner,p
	integer :: itrial1,icall
  real*8 :: En,Omega,sumEkbyk,realk,time !,OMP_get_wtime
  character(1000)::fnum
	character*500::cnpar,formp

!! %-------------------------------------------------------------------
!! %---read initial params ----
	call read_input_params

!---------------------
!$OMP PARALLEL 
!$OMP MASTER
 call get_total_threads(nthreads) 
 if (nthreads .gt. 1) then
    write(*,*) 'Using threads with number of threads = ',nthreads
 else
    write(*,*) 'Not using any threads, nthreads = ',nthreads
 endif
!$OMP END MASTER
!$OMP END PARALLEL
        
	call open_output_files
	if (particl) call open_output_particle

	call initialise_variable
	if (particl) call allocate_particle

	call global_array

	call initial_configuration

	if (particl) call initialize_particle

	call energy_serial 
	En = 0.0d0
	Omega = 0.0d0
	open(unit=110,file='initial_spectra.out',status='unknown')
	do ispectra = 0,nshell
		write(110,*) ispectra,E_Omega(ispectra,1),E_Omega(ispectra,2)
	enddo
	close(110)
	En = sum(E_Omega(:,1))
	Omega = sum(E_Omega(:,2))
!---------------

	iter = 0; cnt = 0      
	t1=OMP_get_wtime()
        write(*,*) 'Entering loop'
	do iouter = 1,maxiter/navg
		count_out=iouter
		do inner = 1,navg

	    iter = iter + 1
            cnt = iter
!       write(*,*) cnt
      if (forcing) call force_serial

      call adbsh_serial

!! first go to the CM frame of the fluid 
			Vk1(1,1,1) = 0.0d0
			Vk2(1,1,1) = 0.0d0
			Vk3(1,1,1) = 0.0d0 !this sets the real part 
!!to be zero, the imaginary part should be zero anyway, otherwise 
!! the real space velocities wont be real. But we still set them to 
!! zero.
			Vk1(2,1,1) = 0.0d0
			Vk2(2,1,1) = 0.0d0
			Vk3(2,1,1) = 0.0d0

!! and after each call the forcing function is enforced. 
			
			E_Omega = 0.0d0
			call energy_serial
			En = 0.0d0
			Omega = 0.0d0
			SumEkbyk=0.0d0
			En = sum(E_Omega(:,1))
			Omega = sum(E_Omega(:,2))
			tlr_micro_scale = dsqrt(5.0d0*En/Omega)
			vkrms = dsqrt(2.0d0*En/3.0d0)
!! The Taylor Micro Scale Reynolds number 
			Rlambda = vkrms*tlr_micro_scale/vis
      c_cfl=umax*dt_by_dx                         
			time = (iter*delta*vkrms)/(2.0d0*pi)	
			write(200,*)time,En,Omega  
      write(201,*)Rlambda,c_cfl
             

		enddo
!! writing the spectra after every navg number of steps
!!------------------------------------------------------


  write(fnum,'(g8.0)') iouter
 open(unit=224,file='spectra/spectra'//trim(adjustl(fnum))//'.out',status='unknown')
		do ispectra= 0,nshell
	write(224,*)ispectra,E_Omega(ispectra,1),E_Omega(ispectra,2)
		enddo
    close(224)

!!------------------------------------------------------
!	open(unit=11,file='vel/Vk'//trim(adjustl(fnum))//'.in',form='unformatted',status='unknown')
!	write(11)(((Vk1(i1,i2,i3),Vk2(i1,i2,i3),Vk3(i1,i2,i3), &
!         i1=1,n1+2),i2=1,n2),i3=1,n3)
!	close(11)
	enddo
	t2=OMP_get_wtime()
     
!! First write down the velocity arrays. 
	open(unit=11,file='Vk.in',form='unformatted',status='unknown')
	write(11)(((Vk1(i1,i2,i3),Vk2(i1,i2,i3),Vk3(i1,i2,i3), &
         i1=1,n1+2),i2=1,n2),i3=1,n3)
	close(11)
!!  and \box(t - \delta t) from file. 
	open(unit=12,file='VXW.in',form='unformatted',status='unknown')
	write(12)(((VWk1(i1,i2,i3),VWk2(i1,i2,i3),VWk3(i1,i2,i3), &
         i1=1,n1+2),i2=1,n2),i3=1,n3)
	close(12)

	En = 0.0d0
	Omega = 0.0d0
	sumEkbyk = 0.0d0
	write(110,*)
	do ispectra = 1,nshell
		En = En + E_Omega(ispectra,1)
		Omega = Omega + E_Omega(ispectra,2)
		sumEkbyk = sumEkbyk + E_Omega(ispectra,1)/dfloat(ispectra)  
		write(110,*) ispectra,E_Omega(ispectra,1),E_Omega(ispectra,2)
	enddo
	close(110)
!! The total energy and dissipation is then used to calculate 
!! the other averaged quantities which characterise the simulation.
!! The energy dissipation rate $\epsilon = 2*\nu \Omega $         
	edr = 2.0*vis*Omega
!! The Taylor Micro scale 
	tlr_micro_scale = dsqrt(5.0d0*En/Omega)
!! The dissipation wave number, 
	dissipation_scale = (edr/(vis**3))**(0.25)
!! The root mean square velocity, 
	vkrms = dsqrt(2.0d0*En/3.0d0)
!! The Taylor Micro Scale Reynolds number 
	Rlambda = vkrms*tlr_micro_scale/vis
!! the integral scale, 
	integral_scale = (3.0d0*pi*sumEkbyk)/(4.0d0*En)
!! and the large eddy turnover time. 
	tau_eddy = (integral_scale/vkrms)/delta
!! Defination of Integral Scale is taken from Ref.(\cite{pope}) page 240, 
!! Eq (6.259) and Eq (6.260). 
!! The parameter characterising the simulation are written down in the
!! file parameters.out 
	open(unit=112,file='parameters.out',status='unknown')
	write(112,*) 'viscosity (nu)    :',vis,' hyperviscosity (\nu_h)  :',vis2
	write(112,*) 'energy diss-rate(epsilon)  :',edr
	write(112,*) 'Mean Energy                 :',En
	write(112,*) 'Mean Enstrophy(Omega)      :',Omega
	write(112,*) 'Box-Length                  :',length	
	write(112,*) 'Taylor-Micro-Scale (lambda):',tlr_micro_scale
	write(112,*) 'Dissipation Scale(k_d)      :',dissipation_scale
	write(112,*) 'rms-velocity(u_{rms})       :',vkrms
	write(112,*) 'Re_Lambda (Re_{lambda})    :',Rlambda
	write(112,*) 'Integral Scale(L)           :',integral_scale
	write(112,*) 'large eddy turnover time (tau_{eddy}) :',tau_eddy
	write(112,*) 'no. of itn for tau_eddy',tau_eddy/delta
	write(112,*) 'averaging done over',tau_eddy/(delta*maxiter),'tau_{eddy} times'
	close(112)
	close(200)
	close(201)
	close(202)
!!	close(13)
!!	close(14)
!! Finally we deallocate the arrays allocated,  
	deallocate(den_state)
	deallocate(time_increment)
	deallocate(Vk1,Vk2,Vk3)
	deallocate(E_Omega)
!! Destroy the plans created (for FFTW only ) 
	call rfftwnd_f77_destroy_plan(pfor)
	call rfftwnd_f77_destroy_plan(pinv)

  t=t2-t1
  print*,"Time taken=",t
  
!! %-------------*********************-------------------------------
!! and end our program 			
	end program spectral_serial 
	
