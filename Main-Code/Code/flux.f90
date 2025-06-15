!!         serial_nonlin.f90
!! This subroutine evaluates the non-linear part of the navier-stokes
!! equation. This is a serial code.This IS the part that contains
!! call to fourier transform subroutines. This SHOULD be changed
!! from one serial machine to another.
!! -----------------------------------------------------------------
			subroutine flux
	use mod_serial_fluid
	use mod_particle
	use mod_initial
	implicit none
	real*8 :: cvr1,cvc1,cvr2,cvc2,cvr3,cvc3,w1,w2,w3,rv1,rv2,rv3
        real*8 :: ptr1,ptc1,ptr2,ptc2,ptr3,ptc3,mode_tr
        real*8 :: tr1,tr2,tr3,tc1,tc2,tc3 
        real*8 :: p11,p22,p33,p12,p23,p31 
        real*8 :: rk2,rk2inv,rk
        real*8 :: one_by_ncube,volume_element,energy_norm
        real*8,allocatable,dimension(:,:,:)::Vk1_dum,Vk2_dum,Vk3_dum
        real*8,allocatable,dimension(:,:,:)::Wk1_dum,Wk2_dum,Wk3_dum
	integer ::i1,i2,i3,k1,k2,k3,ireal,iimag,j,i4
        integer :: msh1,ispectra,ksqr
	character(1000)::fnum1        
allocate(Vk1_dum(n1+2,n2,n3),Vk2_dum(n1+2,n2,n3),Vk3_dum(n1+2,n2,n3))
allocate(Wk1_dum(n1+2,n2,n3),Wk2_dum(n1+2,n2,n3),Wk3_dum(n1+2,n2,n3))
       	one_by_ncube = 1.0d0/dfloat(n1*n2*n3)
	volume_element = one_by_ncube*(length**3.0d0)
	energy_norm = one_by_ncube*one_by_ncube 
        tr_shell = 0.0d0
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
!! Copying it in dummy array
                                Vk1_dum(ireal,i2,i3)=Vk1(ireal,i2,i3)
                                Vk1_dum(iimag,i2,i3)=Vk1(iimag,i2,i3)
                                Vk2_dum(ireal,i2,i3)=Vk2(ireal,i2,i3)
                                Vk2_dum(iimag,i2,i3)=Vk2(iimag,i2,i3)
                                Vk3_dum(ireal,i2,i3)=Vk3(ireal,i2,i3)
                                Vk3_dum(iimag,i2,i3)=Vk3(iimag,i2,i3)
				cvr1=Vk1(ireal,i2,i3)
				cvc1=Vk1(iimag,i2,i3)
!! -------------
				cvr2=Vk2(ireal,i2,i3)
				cvc2=Vk2(iimag,i2,i3)
!! -------------         
				cvr3=Vk3(ireal,i2,i3)
				cvc3=Vk3(iimag,i2,i3)
!! ------------------------------------------
				Wk1_dum(ireal,i2,i3) = -factor*(k2*cvc3 - k3*cvc2)
				Wk1_dum(iimag,i2,i3) =  factor*(k2*cvr3 - k3*cvr2)
!! 
				Wk2_dum(ireal,i2,i3) = -factor*(k3*cvc1 - k1*cvc3)
				Wk2_dum(iimag,i2,i3) =  factor*(k3*cvr1 - k1*cvr3)
!!
				Wk3_dum(ireal,i2,i3) = -factor*(k1*cvc2 - k2*cvc1)
				Wk3_dum(iimag,i2,i3) =  factor*(k1*cvr2 - k2*cvr1)
!! -------------------------------------------
			enddo
		enddo
	enddo
!$omp end do
!$omp end parallel

!! -----------------------------------------------------------------
!! step ii : inv-fft Wk to real space, do the inv-fft in-place
!! ------------using FFTW ----------
  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,Wk1_dum, 0)
  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,Wk2_dum, 0)
  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,Wk3_dum, 0)
!! -------normalise etc ------------------------------
	do i1=n1+1,n1+2
		Wk1_dum(i1,:,:) = 0.0d0
		Wk2_dum(i1,:,:) = 0.0d0
		Wk3_dum(i1,:,:) = 0.0d0
	enddo
	Wk1_dum = Wk1_dum*scale
	Wk2_dum = Wk2_dum*scale
	Wk3_dum = Wk3_dum*scale
!! -----------------------------------------------------------------     
!! ------------using FFTW ----------
  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,Vk1_dum, 0)
  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,Vk2_dum, 0)
  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,Vk3_dum, 0)

!! -------normalise etc ------------------------------
	do i1=n1+1,n1+2
		Vk1_dum(i1,:,:) = 0.0d0
		Vk2_dum(i1,:,:) = 0.0d0
		Vk3_dum(i1,:,:) = 0.0d0
	enddo
	Vk1_dum = Vk1_dum*scale
	Vk2_dum = Vk2_dum*scale
	Vk3_dum = Vk3_dum*scale
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
				w1=Wk1_dum(i1,i2,i3)
				w2=Wk2_dum(i1,i2,i3)
				w3=Wk3_dum(i1,i2,i3)
!! ----------
				rv1=Vk1_dum(i1,i2,i3)
				rv2=Vk2_dum(i1,i2,i3)
				rv3=Vk3_dum(i1,i2,i3)
!! ----------
				Wk1_dum(i1,i2,i3) = rv2*w3 - rv3*w2
				Wk2_dum(i1,i2,i3) = rv3*w1 - rv1*w3
				Wk3_dum(i1,i2,i3) = rv1*w2 - rv2*w1
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
  call rfftwnd_f77_threads_one_real_to_complex(nthreads,pfor,Wk1_dum, 0)
  call rfftwnd_f77_threads_one_real_to_complex(nthreads,pfor,Wk2_dum, 0)
  call rfftwnd_f77_threads_one_real_to_complex(nthreads,pfor,Wk3_dum, 0)
!! -----------------------------------------------------------
	do i3 = 1,n3
		k3 = (i3-1) - n3*(i3/(n1hf+1))	
		do i2 = 1,n2
			k2 = (i2-1) - n2*(i2/(n1hf+1))	
			do i1 =1,n1hf
				k1 = i1 -1
				ksqr = k1*k1 + k2*k2 + k3*k3 
				rk2 = factor*factor*dfloat(ksqr)
                                rk=dsqrt(rk2)
                                msh1=nint(rk)
				ireal = 2*i1-1
				iimag=2*i1
!! -------the projection operator --------------				
				if(rk2.gt.1.0E-4) then
				rk2inv= 1.0d0/dfloat(ksqr)
				p11=1-k1*k1*rk2inv
				p22=1-k2*k2*rk2inv
				p33=1-k3*k3*rk2inv
				p12= -k1*k2*rk2inv
				p23= -k2*k3*rk2inv
				p31= -k3*k1*rk2inv
				else
				p11=2.0/3.0
				p22=p11
				p33=p22
				p12= -1.0/3.0
				p23=p12
				p31=p23
				endif
!! ---------------------------------------
				tr1=Wk1_dum(ireal,i2,i3)
				tc1=Wk1_dum(iimag,i2,i3)
				tr2=Wk2_dum(ireal,i2,i3)
				tc2=Wk2_dum(iimag,i2,i3)
				tr3=Wk3_dum(ireal,i2,i3)
				tc3=Wk3_dum(iimag,i2,i3)
!! ------------------------------------------------
!! APPLYING THE PROJECTION
!! ------------------------------------------------
                                ptr1=  (p11*tr1+ p12*tr2+ p31*tr3)  
				ptc1=  (p11*tc1+ p12*tc2+ p31*tc3) 
				ptr2=  (p12*tr1+ p22*tr2+ p23*tr3)
				ptc2=  (p12*tc1+ p22*tc2+ p23*tc3) 
				ptr3=  (p31*tr1+ p23*tr2+ p33*tr3)
				ptc3=  (p31*tc1+ p23*tc2+ p33*tc3)
!! ---------------------------------------------------
!! CALCULATING THE TRANSFER TERM FOR MODE 'k'
!! ---------------------------------------------------
                                mode_tr =  Vk1(ireal,i2,i3)*ptr1+Vk1(iimag,i2,i3)*ptc1 &
                                        + Vk2(ireal,i2,i3)*ptr2+Vk2(iimag,i2,i3)*ptc2 &
                                        + Vk3(ireal,i2,i3)*ptr3+Vk3(iimag,i2,i3)*ptc3
!! ---------------------------------------------------
                                !! SHELL AVERAGING
!! ---------------------------------------------------
                                if((k1.eq.0).or.(k1.eq.n1h))then
                                tr_shell(msh1,1)=tr_shell(msh1,1)+mode_tr
                                else
                                tr_shell(msh1,1)=tr_shell(msh1,1)+2.0d0*mode_tr
                                endif
			enddo
		enddo
	enddo
        tr_shell = energy_norm*tr_shell
! -------------------------------------------
! FINDING THE FLUX BY SUMMING THE TRANSFER TERMS, WITH NEG SIGN.
! -------------------------------------------
        do ispectra=0,nshell
             tr_shell(ispectra,2)=-sum(tr_shell(:ispectra,1))
        end do        
deallocate(Vk1_dum)
deallocate(Vk2_dum)
deallocate(Vk3_dum)
deallocate(Wk1_dum)
deallocate(Wk2_dum)
deallocate(Wk3_dum)
        end subroutine flux

