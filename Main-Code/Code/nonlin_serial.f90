!!         serial_nonlin.f90
!! This subroutine evaluates the non-linear part of the navier-stokes
!! equation. This is a serial code.This IS the part that contains
!! call to fourier transform subroutines. This SHOULD be changed
!! from one serial machine to another.
!! -----------------------------------------------------------------
			subroutine nonlin_serial
	use mod_serial_fluid
	use mod_particle
	use mod_initial
	implicit none
	real*8 :: cvr1,cvc1,cvr2,cvc2,cvr3,cvc3,w1,w2,w3,rv1,rv2,rv3
	integer ::i1,i2,i3,k1,k2,k3,ireal,iimag,j,i4
	character(1000)::fnum1        
!! ---variables defined for fourier transform -----
!! --- everest  -----------------------------
!	real*8,dimension(n1+15+ 2*(n2+15)+ 2*(n3+15)) :: coeff 
!	integer :: sign,ncube
!! ---alpha machines ------------------------------
!	integer :: status
!! --for fourier transforming in the alpha machines : 
!	include '/usr1/dhruba/include/DXMLDEF.FOR' ! in theory machines
!	include '/home/phd/97/phymitra/include/DXMLDEF.FOR' !for SERC  
!! -----------------------------------------------------------------
!!    ***************************************************************
!!    * This uses           FFTW       to fourier transform         *
!!    * pfor for forward FFT                                        *
!!    * pinv for inverse FFT                                        *
!!    * Plan variables etc are defined in mod_serial_fluid and      *
!!    * assigned in serial_spectral(main)                           *
!!    * see fftw documentation for further details                  *
!!    ***************************************************************
!! step i : evaluate W from Vk and write it in W 
!! --------------------------------------------------------------
!$omp parallel &
!$omp shared(Wk1,Wk2,Wk3,Vk1,Vk2,Vk3,factor,n2,n3,n1hf)&
!$omp private(i1,i2,i3,k1,k2,k3,ireal,iimag,cvr1,cvc1,cvr2,cvc2,cvr3,cvc3) 
!$omp do

	do i3 = 1,n3
		k3 = (i3-1) - n3*(i3/(n1hf+1))	
		do i2 = 1,n2
			k2 = (i2-1) - n2*(i2/(n1hf+1))	
			do i1 =1,n1hf
				k1 = i1 -1
				ireal = 2*i1 -1
				iimag = 2*i1
!! --------------
				cvr1=Vk1(ireal,i2,i3)
				cvc1=Vk1(iimag,i2,i3)
!! -------------
				cvr2=Vk2(ireal,i2,i3)
				cvc2=Vk2(iimag,i2,i3)
!! -------------         
				cvr3=Vk3(ireal,i2,i3)
				cvc3=Vk3(iimag,i2,i3)
!! ------------------------------------------
				Wk1(ireal,i2,i3) = -factor*(k2*cvc3 - k3*cvc2)
				Wk1(iimag,i2,i3) =  factor*(k2*cvr3 - k3*cvr2)
!! 
				Wk2(ireal,i2,i3) = -factor*(k3*cvc1 - k1*cvc3)
				Wk2(iimag,i2,i3) =  factor*(k3*cvr1 - k1*cvr3)
!!
				Wk3(ireal,i2,i3) = -factor*(k1*cvc2 - k2*cvc1)
				Wk3(iimag,i2,i3) =  factor*(k1*cvr2 - k2*cvr1)
!! -------------------------------------------
			enddo
		enddo
	enddo
!$omp end do
!$omp end parallel

	call dermat_serial
!! -----------------------------------------------------------------
!! step ii : inv-fft Wk to real space, do the inv-fft in-place
!! ------------using FFTW ----------
  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,Wk1, 0)
  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,Wk2, 0)
  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,Wk3, 0)
!! -------normalise etc ------------------------------
	do i1=n1+1,n1+2
		Wk1(i1,:,:) = 0.0d0
		Wk2(i1,:,:) = 0.0d0
		Wk3(i1,:,:) = 0.0d0
	enddo
	Wk1 = Wk1*scale
	Wk2 = Wk2*scale
	Wk3 = Wk3*scale
!! -----------------------------------------------------------------     
!! ------------using FFTW ----------
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
!------------- Calculate energy injection rate-----------------------
       umax1=maxval(abs(Vk1))
       umax2=maxval(abs(Vk2))
       umax3=maxval(abs(Vk3))
       umax=max(umax1,umax2,umax3)

!!	if (mod(cnt,10)==0) then
!!	i4 = 1
!!	do i3 = 1, n3,8
!!		do i2 = 1, n2,8
!!			do i1 = 1, n1,8
!!
!!			uxe(i4) = Vk1(i1,i2,i3)			
!!			uye(i4) = Vk2(i1,i2,i3)			
!!			uze(i4) = Vk3(i1,i2,i3)			
!!
!!			dx_uxe(i4) = dx_ux(i1,i2,i3)			
!!			dx_uye(i4) = dx_uy(i1,i2,i3)			
!!			dx_uze(i4) = dx_uz(i1,i2,i3)			
!!			dy_uxe(i4) = dy_ux(i1,i2,i3)			
!!			dy_uye(i4) = dy_uy(i1,i2,i3)			
!!			dy_uze(i4) = dy_uz(i1,i2,i3)			
!!			dz_uxe(i4) = dz_ux(i1,i2,i3)			
!!			dz_uye(i4) = dz_uy(i1,i2,i3)			
!!			dz_uze(i4) = dz_uz(i1,i2,i3)			
!!			i4 = i4 + 1
!!			enddo
!!		enddo
!!	enddo
!	call write_output_euler
!!	endif

	if (particl) then
                if (mod(cnt,10)==0 .OR. cnt==1) call write_output_particle
		call dermat_at_particle
		call evolv_particle
		call get_dAdt
          
	endif

! Writing the real space velocities to the fil after every navg steps
        
       if(mod(cnt,navg)==0) then
         
         write(fnum1,'(g8.0)') count_out
         open(unit=11,file='vel/Vk'//trim(adjustl(fnum1))//'.in',form='unformatted',status='unknown')
         write(11)(((Vk1(i1,i2,i3),Vk2(i1,i2,i3),Vk3(i1,i2,i3), &
        i1=1,n1),i2=1,n2),i3=1,n3)
         close(11)

      endif

!! ------------------------------------------------------------------
!! step iv : find VXOmega in real space , and overwrite it in the 
!! array Wk ------------------------------------------------------ 

!!	write(*,*) sum(Vk1*Vk1+Vk2*Vk2+Vk3*Vk3)/dfloat(n1*n2*n3)
!$omp parallel &
!$omp shared(Wk1,Wk2,Wk3,Vk1,Vk2,Vk3,n1,n2,n3)&
!$omp private(i1,i2,i3,w1,w2,w3,rv1,rv2,rv3) 
!$omp do

	do i3=1,n3
		do i2=1,n2
			do i1=1,n1
!! ---------        
				w1=Wk1(i1,i2,i3)
				w2=Wk2(i1,i2,i3)
				w3=Wk3(i1,i2,i3)
!! ----------
				rv1=Vk1(i1,i2,i3)
				rv2=Vk2(i1,i2,i3)
				rv3=Vk3(i1,i2,i3)
!! ----------
				Wk1(i1,i2,i3) = rv2*w3 - rv3*w2
				Wk2(i1,i2,i3) = rv3*w1 - rv1*w3
				Wk3(i1,i2,i3) = rv1*w2 - rv2*w1
!! ---------------
			enddo
		enddo
	enddo
!$omp end do
!$omp end parallel
!! -----------------------------------------------------------------
!! step v : fft VXOmega to fourier space, in-place .

!! ------------using FFTW ----------
!!       Forward Transform 
  call rfftwnd_f77_threads_one_real_to_complex(nthreads,pfor,Wk1, 0)
  call rfftwnd_f77_threads_one_real_to_complex(nthreads,pfor,Wk2, 0)
  call rfftwnd_f77_threads_one_real_to_complex(nthreads,pfor,Wk3, 0)

	if (forcing) then
! ------------Adding contribution from forcing----------------   
!$omp parallel &
!$omp shared(Wk1,Wk2,Wk3,fk1,fk2,fk3,forced_modes1,count1)&
!$omp private(i1,i2,i3,j,ireal,iimag) 
!$omp do
        do j=1,count1 !shell 1 
                i1 = forced_modes1(j,1)
                ireal = 2*i1-1
                iimag = 2*i1
                i2 = forced_modes1(j,2)
                i3 = forced_modes1(j,3)
        
                Wk1(ireal,i2,i3) = Wk1(ireal,i2,i3) + fk1(ireal,i2,i3)
                Wk1(iimag,i2,i3) = Wk1(iimag,i2,i3) + fk1(iimag,i2,i3)  

                Wk2(ireal,i2,i3) = Wk2(ireal,i2,i3) + fk2(ireal,i2,i3)
                Wk2(iimag,i2,i3) = Wk2(iimag,i2,i3) + fk2(iimag,i2,i3)
!!
                Wk3(ireal,i2,i3) = Wk3(ireal,i2,i3) + fk3(ireal,i2,i3)         
                Wk3(iimag,i2,i3) = Wk3(iimag,i2,i3) + fk3(iimag,i2,i3)  
        enddo
!$omp end do
!$omp end parallel
!------------------------------------------------------------------------------
!$omp parallel &
!$omp shared(Wk1,Wk2,Wk3,fk1,fk2,fk3,forced_modes2,count2)&
!$omp private(i1,i2,i3,j,ireal,iimag) 
!$omp do

        do j=1,count2 !shell 2 
        i1 = forced_modes2(j,1)
                ireal = 2*i1-1
                iimag = 2*i1
                i2 = forced_modes2(j,2)
                i3 = forced_modes2(j,3)  !! Local to processor

                Wk1(ireal,i2,i3) = Wk1(ireal,i2,i3) + fk1(ireal,i2,i3)
                Wk1(iimag,i2,i3) = Wk1(iimag,i2,i3) + fk1(iimag,i2,i3)
                
                Wk2(ireal,i2,i3) = Wk2(ireal,i2,i3) + fk2(ireal,i2,i3)
                Wk2(iimag,i2,i3) = Wk2(iimag,i2,i3) + fk2(iimag,i2,i3)
!!
                Wk3(ireal,i2,i3) = Wk3(ireal,i2,i3) + fk3(ireal,i2,i3)
                Wk3(iimag,i2,i3) = Wk3(iimag,i2,i3) + fk3(iimag,i2,i3)
        enddo
!$omp end do
!$omp end parallel
	endif
                call energy_injection
!! ----------------------------------------------------------------
			end subroutine nonlin_serial

