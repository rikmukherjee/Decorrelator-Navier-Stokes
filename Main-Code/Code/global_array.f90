!!        global_array.f90
!! This subroutine evaluates the global arrays used in the serial 
!! spectral code. This doesnot use any system call thus should
!! remain unchanged from one serial machine to another. 
!! changed for slaved adam-bashforth scheme. 
!! -----------------------------------------------------------------
			subroutine global_array
	use mod_serial_fluid
	implicit none
!! defination of other variables 
	integer :: l,j,k,k1,k2,k3,ireal,iimag,ksqr,mshl,&
               mode_no1,mode_no2,mode_no3,itrial1,jtrial
	real*8 :: rk,rk2,hypvis
!! --the density of states calculations ----
!! --the formula below does the folding in from i to j 
!! where i varies from 1 to n1 and j should be from -n1/2 + 1 
!! to + n1/2 : 
!!          j = (i -1) - n1*(i/(n1hf+1))
!! this has been checked in the program modulo.f in home/res/fortran
!! in the theory machines.  
!! -----------------------------------------
	do k=1,n3
		k3 = (k-1) - n3*(k/(n1hf+1))	
		do j=1,n2
			k2 = (j-1) - n2*(j/(n1hf+1))	
			do l=1,n1hf
				k1=l-1
!! -------------
				ksqr = k1*k1 + k2*k2 + k3*k3
				rk2 = dfloat(ksqr)*factor*factor
				rk = dsqrt(rk2)
				mshl = nint(rk)
				hypvis = vis + vis2*rk2
				time_increment(ksqr)= dexp(-delta*hypvis*rk2)
!! -------------
				if((k1.eq.0).or.(k1.eq.n1h)) then
				den_state(mshl)=den_state(mshl) + 1
				else
				den_state(mshl)=den_state(mshl) + 2
				endif
			enddo
		enddo
	enddo
!! ----now list the forced modes  -----
	mode_no1 = den_state(1)
	mode_no2 = den_state(2)	
	mode_no3 = den_state(3)
	allocate(forced_modes1(mode_no1,4),forced_modes2(mode_no2,4), & 
                                      forced_modes3(mode_no3,4))
	forced_modes1 = 0.0d0
	forced_modes2 = 0.0d0
	forced_modes3 = 0.0d0
!! -------------------------------------------------------
	count1 = 0 
	count2 = 0
	count3 = 0
	do k=1,n3
		k3 = (k-1) - n3*(k/(n1hf+1))	
		do j=1,n2
			k2 = (j-1) - n2*(j/(n1hf+1))	
			do l=1,n1hf
				k1=l-1
!! -------------
				ksqr = k1*k1 + k2*k2 + k3*k3
				rk2 = dfloat(ksqr)*factor*factor
				rk = dsqrt(rk2)
				mshl = nint(rk)
!! ----------------
				if(mshl.eq.1) then
				count1 = count1 + 1
				forced_modes1(count1,1) = l
				forced_modes1(count1,2) = j
				forced_modes1(count1,3) = k
				if((k1.eq.0).or.(k1.eq.n1/2))then
				forced_modes1(count1,4) = 1
				else
				forced_modes1(count1,4) = 2
				endif
				else
				endif

!! ----------------
				if(mshl.eq.2) then
				count2 = count2 + 1
				forced_modes2(count2,1) = l
				forced_modes2(count2,2) = j
				forced_modes2(count2,3) = k
				if((k1.eq.0).or.(k1.eq.n1/2))then
				forced_modes2(count2,4) = 1
				else
				forced_modes2(count2,4) = 2
				endif
				else
				endif
!! ----------------
				if(mshl.eq.3) then
				count3 = count3 + 1
				forced_modes3(count3,1) = l
				forced_modes3(count3,2) = j
				forced_modes3(count3,3) = k
				if((k1.eq.0).or.(k1.eq.n1/2))then
				forced_modes3(count3,4) = 1
				else
				forced_modes3(count3,4) = 2
				endif
				else
				endif
!! -------------
			enddo
		enddo
	enddo
!! -------------checks ----------------------------------------
	open(unit=209,file='density_of_states.out',status='unknown')
	do itrial1 = 0,nshell
		write(209,*)itrial1,den_state(itrial1)
	enddo
	write(209,*)count1,count2,count3
	write(209,*) 'mode 1'
	do itrial1 = 1,count1
		write(209,*) itrial1,(forced_modes1(itrial1,jtrial),jtrial=1,4)
	enddo
	write(209,*) 'mode 2'
	do itrial1 = 1,count2
		write(209,*) itrial1,(forced_modes2(itrial1,jtrial),jtrial=1,4)
	enddo
	write(209,*) 'mode 3'
	do itrial1 = 1,count3
		write(209,*) itrial1,(forced_modes3(itrial1,jtrial),jtrial=1,4)
	enddo
	close(209)	
!! -------------------------------------------------------------------
				
			end subroutine global_array
