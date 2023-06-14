program bfss

  implicit none
  integer dim,lambda,N,idum,ngauge,dtau_phi,dtau_theta,isweep,jsweeps,nskip,ntau
  real a,coupling,temperature,mass
  parameter (dim=9)
  parameter (lambda=10)
  parameter (N=10)
  double precision phi(1:dim,1:lambda,1:N,1:N),theta(1:N)
  double precision P_phi(1:dim,1:lambda,1:N,1:N),P_theta(1:N)
  double precision force_phi(1:dim,1:lambda,1:N,1:N),force_theta(1:N)
  double precision action_boson,ham,variationH,accept,reject,acceptance,ncv
  parameter (a=1.0d0,mass=1.0d0,coupling=1.0d0,temperature=1.0d0,dtau_phi=1.0d0,dtau_theta=1.0d0,ntau=5.0d0)
  parameter (ncv=0.0d0,idum=5000,ngauge=0.0d0,jsweeps=100)

!!!!!!!!!!! INPUTS !!!!!!!!!!!
  !  ngauge=1.0d0 ! 0=ungauged,1=gauged !
  ! ntau =5.0d0 ! Number of time leap frog is repeated
  !  dtau_theta=1.0d0 !dtau for theta variable in leap frog (molecular dynamics) 
  !  dtau_phi=1.0d0 ! dtau for phi variable in leap frog (molecular dynamics)
  !  idum=50000 ! seed for random number generator ran2(idum)
  !  temperature=1.0d0 ! fixing the temperature
  !  ncv=0.0d0 ! number of the constraint violation
  ! jsweeps= total number of sweeps, isweep=actual sweep 1<isweep<jsweeps

  isweep=1
  do while (isweep .le. jsweeps)

     call metropolis(idum,mass,coupling,N,dim,lambda,ngauge,acceptance,force_phi,force_theta,&
          &accept,reject,P_phi,P_theta,ncv,ham,action_boson,variationH,Phi,theta)
     print*,ham
     isweep=isweep+1
     
  end do
  


end program bfss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! ACTION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine action(idum,dim,N,lambda,a,coupling,mass,action_boson)
 
  integer d,N,j,lambda,dim,idum,l,i
  real a,y,ran2,y1,y2,y3,coupling,mass
  double precision phi(1:dim,1:lambda,1:N,1:N),rr,angle,theta(1:N),phase(1:N),phase1(1:N)
  double precision S_FP,action_boson,kinetic,uxumx(1:N,1:N)
  complex ii
  
  action_boson=0.0d0 
  ii=cmplx(1,0)
  phi=0.0d0
  
  do d=1,dim
     do l=1,lambda
        do i=1,N
           do j=1,N
              rr=rr+phi(d,l,i,j)*phi(d,l,j,i)
           end do
        end do
     end do
  end do

  do d=1,dim
     do l=1,lambda
        l1=l+1 
        do i=1,N
           do j=1,N
              y=ran2(idum)
              y1=ran2(idum)
              theta(i)=y
              theta(j)=y1
              phase(i)=(cos(y)/dble(lambda))+ii*cmplx(sin(y)/dble(lambda))
              phase1(j)=(cos(y1)/dble(lambda))-ii*cmplx(sin(y1)/dble(lambda))
              uxumx(i,j)=phi(d,l,i,j)*phase(i)*phi(d,l1,j,i)*phase1(j)
           end do
        end do
        do i=1,N
           do j=1,N
              kinetic=kinetic+dble(uxumx(i,j)*uxumx(j,i))
           end do
        end do
        do i=1,N
           do j=1,N
              y2=ran2(idum)
              y3=ran2(idum)
              ! We write now the Faddeev-Popov term
              angle=(y2-y3)/2d0
              S_FP=-0.5d0*log(sin(angle)**2)
           end do
        end do
     end do
  end do
  action_boson=S_fp-(N/a)*kinetic+N*((1/a)*rr+(0.5d0*a*(mass**2))*rr)

  print*,'kinetic term=',kinetic
  print*,'radius is =',rr
  print*,'Faddeev-Popov term =',S_fp
    
end subroutine action


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! FORCE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine force(idum,dim,N,lambda,a,coupling,mass,force_phi,force_theta)
  
  integer dim,N,lambda,d,l,idum
  real a,coupling,mass,ran2,y,y1
  double precision theta(1:N),phi(1:dim,1:lambda,1:N,1:N),force_theta(1:N),force_phi(1:dim,1:lambda,1:N,1:N)
  complex ii,phase(1:N),phase1(1:N),uxdu(1:N,1:N),udxu(1:N,1:N)

  ii=cmplx(0.0d0,1.0d0)

  

  force_theta=0.0d0

  do d=1,dim
     do l=1,lambda
        do i=1,N
           do j=1,N
              y=ran2(idum)
              y1=ran2(idum)
              theta(i)=y
              theta(j)=y1
              phase(i)=(cos(y)/dble(lambda))+ii*cmplx(sin(y)/dble(lambda))
              phase1(j)=(cos(y1)/dble(lambda))-ii*cmplx(sin(y1)/dble(lambda))
              if (l .lt. lambda) then
                 force_theta=force_theta-dble(phase(i)*phi(d,l+1,i,j)*phase1(j)*phi(d,l,i,j))
                 force_theta=force_theta+dble(phase(i)*phi(d,l+1,i,j)*phase1(j)*phi(d,l,i,j))
              else
                 force_theta=force_theta-dble(phase(i)*phi(d,1,i,j)*phase1(j)*phi(d,l,i,j))
                 force_theta=force_theta+dble(phase(i)*phi(d,1,i,j)*phase1(j)*phi(d,l,i,j))
              end if
           end do
        end do
     end do
  end do

  force_theta=force_theta*dble(N)/dble(lambda)

  force_phi=0.0d0

  do d=1,dim
     do l=1,lambda
        do i=1,N
           do j=1,N
              y=ran2(idum)
              y1=ran2(idum)
              theta(i)=y
              theta(j)=y1
              phase(i)=(cos(y)/dble(lambda))+ii*cmplx(sin(y)/dble(lambda))
              phase1(j)=(cos(y1)/dble(lambda))-ii*cmplx(sin(y1)/dble(lambda))
              if (l .lt. lambda) then 
                 uxdu(i,j)=phase(i)*phi(d,l+1,i,j)*phase1(j)
              else
                 uxdu(i,j)=phase(i)*phi(d,1,i,j)*phase1(j)
              end if
              if (l .gt. 1) then
                 uxdu(i,j)=conjg(phase(i))*phi(d,l-1,i,j)*conjg(phase1(j))
              else
                 uxdu(i,j)=conjg(phase(i))*phi(d,1,i,j)*conjg(phase1(j))
              end if
           end do
        end do
        do i=1,N
           do j=1,N
              force_phi(d,l,i,j)=force_phi(d,l,i,j)+2.0d0*phi(d,l,i,j)
              force_phi(d,l,i,j)=force_phi(d,l,i,j)-uxdu(i,j)
           end do
        end do
     end do
  end do

  force_phi(d,l,i,j)=force_phi(d,l,i,j)*dble(N)/dble(lambda)
  

end subroutine force

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! HAMILTONIAN !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine hamiltonian(idum,dim,N,lambda,a,coupling,mass,action_boson,ham)

  integer N,dim,lambda,d,l
  real a,mass,coupling
  double precision ham,P_phi(1:dim,1:lambda,1:N,1:N),P_theta(1:N),action_boson

  call action(idum,dim,N,lambda,a,coupling,mass,action_boson)

  ham=action_boson
  do d=1,dim
     do l=1,lambda
        do i=1,N
           do j=1,N
              ham=ham+dble(P_phi(d,l,i,j)*P_phi(d,l,j,i))*0.5d0
           end do
        end do
     end do
  end do
  if (ngauge .eq. 0) then
     do i=1,N
        ham=ham+dble(P_theta(i)*P_theta(i))*0.5d0
     end do
  end if
  print*,'hamiltonian=',ham

end subroutine hamiltonian

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! GENERATE MOMENTUM !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Generate_momentum(idum,N,dim,lambda,P_theta,P_phi)

  integer N,dim,lambda,d,l,idum
  real ran2,r,y,y1,z,ph,ph_y,ph_y1,ph_z
  double precision P_theta(1:N),P_phi(1:dim,1:lambda,1:N,1:N)



  pi=acos(-1.0d0)

  ! We generate the momentum with a gaussian distribution !
  do i=1,N
     r=sqrt(-2.0d0*log(1.0d0-ran2(idum)))
     ph=2.0d0*pi*ran2(idum)
     P_theta(i)=r*cos(ph)
  end do

  do d=1,dim
     do l=1,lambda
        do i=1,N
           do j=1,N
              y=sqrt(-2.0d0*log(1.0d0-ran2(idum)))
              y1=sqrt(-2.0d0*log(1.0d0-ran2(idum)))
              ph_y=2.0d0*pi*ran2(idum)
              ph_y1=2.0d0*pi*ran2(idum)
              P_phi(d,l,i,j)=y*cos(ph_y)
              P_phi(d,l,j,i)=y1*cos(ph_y1)
           end do
        end do
        
        do i=1,N
           z=sqrt(-2.0d0*log(1.0d0-ran2(idum)))
           ph_z=2.0d0*pi*ran2(idum)
           P_phi(d,l,i,i)=z*cos(ph_z)
        end do
        
     end do
  end do

end subroutine Generate_momentum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! CHECK FOR THETA CONSTRAINTS !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine theta_constraints(N,theta,acceptance)

  integer N,i
  double precision theta(1:N),max,min,pi,acceptance

  max=theta(1)
  min=theta(1)
  do imat=2,N
     if(max.lt.theta(i))then
        max=theta(i)
     else if(min.gt.theta(i))then
        min=theta(i)
     end if
  end do

  pi=2d0*dasin(1d0)
  if(max-min.lt.2d0*pi)then
     acceptance=0
  else
     acceptance=1
  end if
  
end subroutine theta_constraints

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! MOLECULAR DYNAMICS !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine molecular_dynamics(idum,N,dim,lambda,ntau,dtau_phi,dtau_theta,ham,force_phi,force_theta,P_phi,phi,theta,P_theta)

  integer i,j,k,ntau,step,N,dim,lambda,d,idum
  real    a,mass,coupling
  double precision ham_init,ham_fin,dtau_phi,dtau_theta,action
  double precision phi(1:dim,1:lambda,1:N,1:N),ham,theta(1:N),action_boson
  double precision P_phi(1:dim,1:lambda,1:N,1:N),P_theta(1:N)
  double precision force_phi(1:dim,1:lambda,1:N,1:N),force_theta(1:N)
  
  call hamiltonian(idum,dim,N,lambda,a,coupling,mass,action_boson,ham)
  
!!! FIRST LEAP FROG STEP !!!
  do d=1,dim
     do l=1,lambda
        do i=1,N
           do j=1,N
              
              phi(d,l,i,j)=phi(d,l,i,j)+P_phi(d,l,i,j)*0.5d0*dtau_phi
              
           end do
        end do
     end do
  end do
      
  theta=theta*P_theta*0.5d0*dtau_theta
  
!!! SECOND LEAP FROP STEP !!!
  step=1
  do while (step .lt. ntau)
     step=step+1
     call force(idum,dim,N,lambda,a,coupling,mass,force_phi,force_theta)
     do d=1,dim
        do l=1,lambda
           do i=1,N
              do j=1,N
                 P_phi(d,l,i,j)=P_phi(d,l,i,j)+force_phi(d,l,i,j)*0.5d0*dtau_phi
              end do
           end do
        end do
     end do
     do d=1,dim
        do n=1,lambda
           do i=1,N
              do j=1,N
                 phi(d,l,i,j)=phi(d,l,i,j)+P_phi(d,l,i,j)*dtau_phi
              end do
           end do
        end do
     end do
  end do
  

!!! LAST LEAP FROG STEP !!!
  call force(idum,dim,N,lambda,a,coupling,mass,force_phi,force_theta)
  
  P_theta=P_theta-force_theta*dtau_theta
  do d=1,dim
     do l=1,lambda
        do i=1,N
           do j=1,N
              P_phi(d,l,i,j)=P_phi(d,l,i,j)-force_phi(d,l,i,j)*dtau_phi
           end do
        end do
     end do
  end do
  
  do d=1,dim
     do l=1,lambda
        do i=1,N
           do j=1,N
              phi(d,l,i,j)=phi(d,l,i,j)+P_phi(d,l,i,j)*0.5d0*dtau_phi
           end do
        end do
     end do
  end do
      
  theta=theta+P_theta*0.5d0*dtau_theta
      
  call hamiltonian(idum,dim,N,lambda,a,coupling,mass,action_boson,ham)


end subroutine molecular_dynamics

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! VARIATION OF ACTION !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine metropolis(idum,mass,coupling,N,dim,lambda,ngauge,acceptance,force_phi,force_theta,&
     &accept,reject,P_phi,P_theta,ncv,ham,action_boson,variationH,Phi,theta)

  integer N,dim,lambda,idum,ntau
  double precision P_phi(1:dim,1:lambda,1:N,1:N),P_theta(1:N),backup_phi(1:dim,1:lambda,1:N,1:N),backup_theta(1:N)
  double precision action_boson,ham,variationH,accept,reject,ncv,acceptance_new,acceptance
  double precision phi(1:dim,1:lambda,1:N,1:N),theta(1:N),force_phi(1:dim,1:lambda,1:N,1:N),force_theta(1:N)
  double precision dtau_theta,dtau_phi
  real ran2,r,mass

  backup_phi=phi
  backup_theta=theta

  call generate_momentum(idum,N,dim,lambda,P_theta,P_phi)

  if (ngauge .eq. 1) then
     P_theta=0.0d0
  end if

  call hamiltonian(idum,dim,N,lambda,a,coupling,mass,action_boson,ham)
  call molecular_dynamics(idum,N,dim,lambda,ntau,dtau_phi,dtau_theta,ham,force_phi,force_theta,P_phi,phi,theta,P_theta)
  call theta_constraints(N,theta,acceptance)
  variationH=-ham
  call hamiltonian(idum,dim,N,lambda,a,coupling,mass,action_boson,ham)
  variationH=variationH+ham

  if (acceptance .eq. 1) then
     ncv=ncv+1
     acceptance_new=1
  end if
  
  if (acceptance .eq. 0) then
     if (variationH .lt. 0.0d0) then
        accept=accept+1.0d0
     else
        probability=exp(-variationH)
        r=ran2(idum)
        if (r .lt. probability) then
           accept=accept+1.0d0
        else
           phi=backup_phi
           theta=backup_theta
           reject=reject+1.0d0
        end if
     end if
  end if

end subroutine metropolis


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! OBSERVABLES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine observables(N,theta,pol)

  integer N
  double precision theta(1:N),pol

  call polyakov(N,theta,pol)
  

end subroutine observables


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! POLYAKOV LOOP !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine polyakov(N,theta,pol)
  
  integer N
  double precision theta(1:N),pol,re_pol,im_pol

  re_pol=0.0d0
  im_pol=0.0d0

  do i=1,N
     re_pol=re_pol+cos(theta(i))
     im_pol=re_pol+sin(theta(i))
  end do

  re_pol=re_pol/dble(N)
  im_pol=im_pol/dble(N)

  pol=sqrt(re_pol**2.0d0+im_pol**2.0d0)

end subroutine polyakov




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! FUNCTION RAN2(idum) !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function ran2(idum)

    integer idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
    real ran2,AM,EPS,RNMX
    parameter    (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
      NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
    integer idum2,j,k,iv(NTAB),iy
    save iv,iy,idum2
    DATA idum2/123456789/, iv/NTAB*0/, iy/0/

    if (idum.le.0) then
       idum=max(-idum,1)
       idum2=idum
       do 11 j=NTAB+8,1,-1

          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11        continue
          iy=iv(1)
       endif
       k=idum/IQ1
       idum=IA1*(idum-k*IQ1)-k*IR1
       if (idum.lt.0) idum=idum+IM1
       k=idum2/IQ2
       idum2=IA2*(idum2-k*IQ2)-k*IR2
       if (idum2.lt.0) idum2=idum2+IM2
       j=1+iy/NDIV
       iy=iv(j)-idum2
       iv(j)=idum
       if(iy.lt.1)iy=iy+IMM1
       ran2=min(AM*iy,RNMX)
       return
     end function ran2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
