      PROGRAM code

      IMPLICIT NONE

      INTEGER i,j,k,l,d,Npart,n,Iread,Ipart,iq,q_max
      REAL*8 q(13)
      REAL*8 dt,a,b
      PARAMETER (d=3,dt=0.008d0,n=4000) ! dt must be the same as for the time series of grad
      PARAMETER (Npart=50000,q_max=13)
      REAL*8 lambda1(Npart),lambda2(Npart),lambda3(Npart)
      REAL*8 last_lambda1(Npart),last_lambda2(Npart),last_lambda3(Npart)
      REAL*8 lambda_max(n),lambda_min(n),lambdaq_max(n,q_max),lambdaq_min(n,q_max)
      CHARACTER*500 :: cn,cnpar,formp,fname,formwrite
      write(cnpar,'(g8.0)') (Npart)
      formp = '('//trim(adjustl(cnpar))//'es16.3E2)'
      write(cn,'(g8.0)') (q_max+1)
      formwrite = '('//trim(adjustl(cn))//'es16.3E2)'

      q(1) = -3.0d0
      q(2) = -2.5d0
      q(3) = -2.0d0
      q(4) = -1.5d0
      q(5) = -1.0d0
      q(6) = -0.5d0
      q(7) = 0.0d0
      q(8) = 0.5d0
      q(9) = 1.0d0
      q(10) = 1.5d0
      q(11) = 2.0d0
      q(12) = 2.5d0
      q(13) = 3.0d0

      lambda1=0.d0
      lambda2=0.d0
      lambda3=0.d0
      last_lambda1=0.d0
      last_lambda2=0.d0
      last_lambda3=0.d0
      lambda_max=0.d0
      lambda_min=0.d0
      lambdaq_max=0.d0
      lambdaq_min=0.d0

open(unit=100,file='Averaged_Lyapunov_D3p00.out',status='unknown')
open(unit=110,file='Averaged_q-Lyapunov-largest_D3p00.out',status='unknown')
open(unit=120,file='Averaged_q-Lyapunov-minimum_D3p00.out',status='unknown')

open(unit=200,file='Lyapunov_D3p00_d1.out',status='old')
open(unit=300,file='Lyapunov_D3p00_d2.out',status='old')
open(unit=400,file='Lyapunov_D3p00_d3.out',status='old')
do k = 2,n
read(200,formp) (lambda1(Iread),Iread=1,Npart)
read(300,formp) (lambda2(Iread),Iread=1,Npart)
read(400,formp) (lambda3(Iread),Iread=1,Npart)
enddo      
close(200)
close(300)
close(400)
last_lambda1 = lambda1
last_lambda2 = lambda2
last_lambda3 = lambda3
open(unit=200,file='Lyapunov_D3p00_d1.out',status='old')
open(unit=300,file='Lyapunov_D3p00_d2.out',status='old')
open(unit=400,file='Lyapunov_D3p00_d3.out',status='old')

      Do k = 2,n
      read(200,formp) (lambda1(Iread),Iread=1,Npart)
      read(300,formp) (lambda2(Iread),Iread=1,Npart)
      read(400,formp) (lambda3(Iread),Iread=1,Npart)
      Do Ipart=1,Npart
      a = max(last_lambda1(Ipart),last_lambda2(Ipart),last_lambda3(Ipart))
      b = min(last_lambda1(Ipart),last_lambda2(Ipart),last_lambda3(Ipart))

      if(a.eq.last_lambda1(Ipart)) then
      lambda_max(k) = lambda_max(k) + lambda1(Ipart)
      do iq = 1,13
      lambdaq_max(k,iq) = lambdaq_max(k,iq) + exp(q(iq)*lambda1(Ipart)*k*dt)
      enddo
      endif      
      if(a.eq.last_lambda2(Ipart)) then
      lambda_max(k) = lambda_max(k) + lambda2(Ipart)
      do iq = 1,13
      lambdaq_max(k,iq) = lambdaq_max(k,iq) + exp(q(iq)*lambda2(Ipart)*k*dt)
      enddo
      endif      
      if(a.eq.last_lambda3(Ipart)) then
      lambda_max(k) = lambda_max(k) + lambda3(Ipart)
      do iq = 1,13
      lambdaq_max(k,iq) = lambdaq_max(k,iq) + exp(q(iq)*lambda3(Ipart)*k*dt)
      enddo
      endif      

      if(b.eq.last_lambda1(Ipart)) then
      lambda_min(k) = lambda_min(k) + lambda1(Ipart)
      do iq = 1,13
      lambdaq_min(k,iq) = lambdaq_min(k,iq) + exp(q(iq)*lambda1(Ipart)*k*dt)
      enddo
      endif      
      if(b.eq.last_lambda2(Ipart)) then
      lambda_min(k) = lambda_min(k) + lambda2(Ipart)
      do iq = 1,13
      lambdaq_min(k,iq) = lambdaq_min(k,iq) + exp(q(iq)*lambda2(Ipart)*k*dt)
      enddo
      endif      
      if(b.eq.last_lambda3(Ipart)) then
      lambda_min(k) = lambda_min(k) + lambda3(Ipart)
      do iq = 1,13
      lambdaq_min(k,iq) = lambdaq_min(k,iq) + exp(q(iq)*lambda3(Ipart)*k*dt)
      enddo
      endif     

      enddo

      do iq = 1,13
      lambdaq_min(k,iq) = lambdaq_min(k,iq)/dble(Npart)
      lambdaq_max(k,iq) = lambdaq_max(k,iq)/dble(Npart)
      enddo
      lambda_max(k) = lambda_max(k)/dble(Npart)
      lambda_min(k) = lambda_min(k)/dble(Npart)
      write(100,*) k*dt,lambda_max(k),lambda_min(k) 
      write(110,formwrite) k*dt,(lambdaq_max(k,iq),iq=1,13)
      write(120,formwrite) k*dt,(lambdaq_min(k,iq),iq=1,13)
      enddo
      close(100)
      close(110)
      close(120)
close(200)
close(300)
close(400)
      STOP
      END
