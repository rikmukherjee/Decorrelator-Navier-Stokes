!*******************************************************************
!*   *     Subroutine to calculate energy injection rate    *     *!
!*******************************************************************

             subroutine energy_injection


	use mod_serial_fluid
	implicit none

	integer::i1,i2,i3
	real*8::f_dot_v,sum_f_dot_v

  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,fk1, 0)
  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,fk2, 0)
  call rfftwnd_f77_threads_one_complex_to_real(nthreads,pinv,fk3, 0)

	 do i1=n1+1,n1+2
		 fk1(i1,:,:) = 0.0d0
		fk2(i1,:,:) = 0.0d0
		fk3(i1,:,:) = 0.0d0
	enddo
	fk1 = fk1*scale
	fk2 = fk2*scale
	fk3 = fk3*scale

	sum_f_dot_v = 0.0d0
!$omp parallel &
!$omp shared(fk1,fk2,fk3,Vk1,Vk2,Vk3,n1,n2,n3) &
!$omp private(i1,i2,i3,f_dot_v) &
!$omp reduction(+ : sum_f_dot_v)
!$omp do

	do i3 = 1,n3
		do i2 = 1,n2
			do i1 = 1,n1

			f_dot_v = fk1(i1,i2,i3)*Vk1(i1,i2,i3)+ &
                                  fk2(i1,i2,i3)*Vk2(i1,i2,i3)+ &
                                  fk3(i1,i2,i3)*Vk3(i1,i2,i3)

			sum_f_dot_v = sum_f_dot_v+f_dot_v

			enddo
		enddo
	enddo
!$omp end do
!$omp end parallel

	sum_f_dot_v = sum_f_dot_v*scale

	write(202,*)sum_f_dot_v


            end subroutine energy_injection
