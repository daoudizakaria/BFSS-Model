!********************************************************************
! Main Program: BFSS_Simulation
!  - Allocates arrays, initializes variables, runs the simulation,
!    and logs results to a file.
!********************************************************************
program BFSS_Simulation
  use bfss_parameters
  use bfss_simulation
  implicit none
  integer :: isweep, idum, ngauge, out_unit
  double precision :: ham, accept, reject, polyakov
  double precision, allocatable :: phi(:,:,:,:), P_phi(:,:,:,:)
  double precision, allocatable :: force_phi(:,:,:,:)
  double precision, allocatable :: theta(:), P_theta(:)
  integer :: i
  
  ! Allocate arrays using dimensions from bfss_parameters
  allocate(phi(dim, lambda, N, N))
  allocate(P_phi(dim, lambda, N, N))
  allocate(force_phi(dim, lambda, N, N))
  allocate(theta(N))
  allocate(P_theta(N))

  ! Initialize arrays to zero
  phi       = 0.0_dp
  P_phi     = 0.0_dp
  force_phi = 0.0_dp
  theta     = 0.0_dp
  P_theta   = 0.0_dp

  ! Set simulation parameters
  idum   = 5000
  ngauge = ngauge_default
  accept = 0.0_dp
  reject = 0.0_dp

  ! Open file for logging simulation results
  out_unit = 10
  open(unit=out_unit, file="simulation_results.dat", status="replace", action="write")
  write(out_unit, '(A)') "Sweep   Hamiltonian      Acceptance      Rejection      PolyakovLoop"

  ! Main simulation loop
  do isweep = 1, jsweeps
    call metropolis(idum, phi, theta, P_phi, P_theta, force_phi, force_theta, &
                      ngauge, accept, reject, ham)
                      
    ! Optionally print to console
    print*, "Sweep ", isweep, ": Hamiltonian = ", ham

    ! Compute Polyakov loop observable at current sweep
    call observables(phi, theta, P_phi, P_theta, polyakov)

    ! Log the results to file (formatted output)
    write(out_unit, '(I5,2X,F12.6,2X,F12.6,2X,F12.6,2X,F12.6)') isweep, ham, accept, reject, polyakov
  end do

  ! Print final observables to console
  print*, "Final Polyakov Loop =", polyakov

  ! Close output file
  close(out_unit)

  ! Deallocate arrays
  deallocate(phi, P_phi, force_phi, theta, P_theta)
end program BFSS_Simulation
!********************************************************************
! Module: bfss_parameters
!  - Contains simulation constants and parameters.
!  - Uses iso_fortran_env to define a double‐precision kind.
!********************************************************************
module bfss_parameters
  use iso_fortran_env, only: dp => real64
  implicit none
  integer, parameter :: dim    = 9
  integer, parameter :: lambda = 10
  integer, parameter :: N      = 10
  integer, parameter :: jsweeps = 100
  integer, parameter :: ntau    = 5

  real(dp), parameter :: a         = 1.0_dp
  real(dp), parameter :: mass      = 1.0_dp
  real(dp), parameter :: coupling  = 1.0_dp
  real(dp), parameter :: temperature = 1.0_dp

  real(dp), parameter :: dtau_phi   = 1.0_dp
  real(dp), parameter :: dtau_theta = 1.0_dp

  integer, parameter :: ngauge_default = 0
end module bfss_parameters

!********************************************************************
! Module: bfss_random
!  - Contains the custom random number generator ran2.
!********************************************************************
module bfss_random
  implicit none
  integer, parameter :: IM1   = 2147483563
  integer, parameter :: IM2   = 2147483399
  integer, parameter :: IMM1  = IM1 - 1
  integer, parameter :: IA1   = 40014
  integer, parameter :: IA2   = 40692
  integer, parameter :: IQ1   = 53668
  integer, parameter :: IQ2   = 52774
  integer, parameter :: IR1   = 12211
  integer, parameter :: IR2   = 3791
  integer, parameter :: NTAB  = 32
  integer, parameter :: NDIV  = 1 + IMM1 / NTAB
  double precision, parameter :: AM   = 1.0d0 / IM1
  double precision, parameter :: EPS  = 1.2d-7
  double precision, parameter :: RNMX = 1.0d0 - EPS
contains
  double precision function ran2(idum)
    implicit none
    integer, intent(inout) :: idum
    integer :: j, k, iv(NTAB), iy, idum2
    data idum2 /123456789/
    data iv /NTAB*0/, iy /0/
    
    if (idum <= 0) then
       idum = max(-idum, 1)
       idum2 = idum
       do j = NTAB+8, 1, -1
          k = idum / IQ1
          idum = IA1*(idum - k*IQ1) - k*IR1
          if (idum < 0) idum = idum + IM1
          if (j <= NTAB) iv(j) = idum
       end do
       iy = iv(1)
    end if
    k = idum / IQ1
    idum = IA1*(idum - k*IQ1) - k*IR1
    if (idum < 0) idum = idum + IM1
    k = idum2 / IQ2
    idum2 = IA2*(idum2 - k*IQ2) - k*IR2
    if (idum2 < 0) idum2 = idum2 + IM2
    j = 1 + mod(iy, NDIV)
    iy = iv(j) - idum2
    iv(j) = idum
    if (iy < 1) iy = iy + IMM1
    ran2 = min(AM * iy, RNMX)
  end function ran2
end module bfss_random

!********************************************************************
! Module: bfss_simulation
!  - Contains all subroutines for the simulation: action,
!    force, hamiltonian, momentum generation, molecular dynamics,
!    metropolis update, and observables.
!********************************************************************
module bfss_simulation
  use bfss_parameters
  use bfss_random
  implicit none
contains

  !------------------------------------------------------------------
  ! Subroutine: action
  !  - Computes the bosonic action including kinetic and Faddeev-Popov terms.
  !------------------------------------------------------------------
  subroutine action(idum, phi, theta, action_boson)
    implicit none
    integer, intent(inout) :: idum
    double precision, intent(in) :: phi(dim, lambda, N, N)
    double precision, intent(inout) :: theta(N)
    double precision, intent(out) :: action_boson
    integer :: d, l, i, j, l1
    double precision :: rr, kinetic, S_FP, y, angle, y2, y3
    double precision, dimension(N, N) :: uxumx
    complex(dp) :: ii, phase(N), phase1(N)
    
    rr      = 0.0d0
    kinetic = 0.0d0
    S_FP    = 0.0d0
    ii      = cmplx(1.0d0, 0.0d0, kind=dp)
    
    ! Compute radius term (rr) from field phi
    do d = 1, dim
      do l = 1, lambda
        do i = 1, N
          do j = 1, N
            rr = rr + phi(d, l, i, j) * phi(d, l, j, i)
          end do
        end do
      end do
    end do

    ! Compute kinetic and Faddeev-Popov contributions
    do d = 1, dim
      do l = 1, lambda
        l1 = mod(l, lambda) + 1  ! cyclic index: if l==lambda then l1 becomes 1
        do i = 1, N
          y = ran2(idum)
          theta(i) = y   ! Update for theta; verify with your model
          phase(i) = cmplx(cos(y)/dble(lambda), sin(y)/dble(lambda), kind=dp)
        end do
        do i = 1, N
          do j = 1, N
            y2 = ran2(idum)
            y3 = ran2(idum)
            phase1(j) = cmplx(cos(y2)/dble(lambda), -sin(y2)/dble(lambda), kind=dp)
            uxumx(i, j) = phi(d, l, i, j) * phase(i) * phi(d, l1, j, i) * phase1(j)
          end do
        end do
        do i = 1, N
          do j = 1, N
            kinetic = kinetic + dble(uxumx(i, j) * uxumx(j, i))
          end do
        end do
        do i = 1, N
          do j = 1, N
            y2 = ran2(idum)
            y3 = ran2(idum)
            angle = (y2 - y3) / 2.0d0
            if (sin(angle)**2 > 1.0e-12d0) then
              S_FP = S_FP - 0.5d0 * log(sin(angle)**2)
            end if
          end do
        end do
      end do
    end do

    action_boson = S_FP - (N/a) * kinetic + N * ((1.0d0/a) * rr + 0.5d0 * a * (mass**2) * rr)
  end subroutine action

  !------------------------------------------------------------------
  ! Subroutine: force
  !  - Computes the forces (gradients) with respect to phi and theta.
  !------------------------------------------------------------------
  subroutine force(idum, phi, theta, force_phi, force_theta)
    implicit none
    integer, intent(in) :: idum
    double precision, intent(in) :: phi(dim, lambda, N, N), theta(N)
    double precision, intent(out) :: force_phi(dim, lambda, N, N), force_theta(N)
    integer :: d, l, i, j
    double precision :: y, y1, temp
    complex(dp) :: ii, phase(N), phase1(N)
    
    ii = cmplx(0.0d0, 1.0d0, kind=dp)
    force_theta = 0.0d0

    ! Compute force for theta (example derivative; adjust as needed)
    do d = 1, dim
      do l = 1, lambda
        do i = 1, N
          y = ran2(idum)
          phase(i) = cmplx(cos(y)/dble(lambda), sin(y)/dble(lambda), kind=dp)
        end do
        do i = 1, N
          do j = 1, N
            y1 = ran2(idum)
            phase1(j) = cmplx(cos(y1)/dble(lambda), -sin(y1)/dble(lambda), kind=dp)
            force_theta(i) = force_theta(i) + dble(phase(i)*phi(d, l, i, j)*phase1(j)*phi(d, l, i, j))
          end do
        end do
      end do
    end do
    force_theta = force_theta * dble(N) / dble(lambda)

    ! Compute force for phi as a derivative term (using nearest-neighbor difference)
    force_phi = 0.0d0
    do d = 1, dim
      do l = 1, lambda
        do i = 1, N
          do j = 1, N
            if (l < lambda) then
              temp = phi(d, l+1, i, j)
            else
              temp = phi(d, 1, i, j)
            end if
            force_phi(d, l, i, j) = 2.0d0 * phi(d, l, i, j) - temp
          end do
        end do
      end do
    end do
    force_phi = force_phi * dble(N) / dble(lambda)
  end subroutine force

  !------------------------------------------------------------------
  ! Subroutine: hamiltonian
  !  - Computes the Hamiltonian as the sum of the bosonic action and
  !    kinetic terms from the momenta.
  !------------------------------------------------------------------
  subroutine hamiltonian(idum, phi, theta, P_phi, P_theta, ngauge, action_boson, ham)
    implicit none
    integer, intent(in) :: idum, ngauge
    double precision, intent(in) :: phi(dim, lambda, N, N), theta(N)
    double precision, intent(in) :: P_phi(dim, lambda, N, N), P_theta(N)
    double precision, intent(out) :: action_boson, ham
    integer :: d, l, i, j
    
    call action(idum, phi, theta, action_boson)
    ham = action_boson
    do d = 1, dim
      do l = 1, lambda
        do i = 1, N
          do j = 1, N
            ham = ham + 0.5d0 * P_phi(d, l, i, j) * P_phi(d, l, j, i)
          end do
        end do
      end do
    end do
    if (ngauge == 0) then
      do i = 1, N
        ham = ham + 0.5d0 * P_theta(i)**2
      end do
    end if
  end subroutine hamiltonian

  !------------------------------------------------------------------
  ! Subroutine: Generate_momentum
  !  - Generates momentum variables with a Gaussian distribution.
  !------------------------------------------------------------------
  subroutine Generate_momentum(idum, P_theta, P_phi)
    implicit none
    integer, intent(in) :: idum
    double precision, intent(out) :: P_theta(N), P_phi(dim, lambda, N, N)
    integer :: i, d, l, j
    double precision :: r, ph, y, y1, ph_y, ph_y1, z, ph_z, pi
    pi = acos(-1.0d0)
    do i = 1, N
      r = sqrt(-2.0d0 * log(1.0d0 - ran2(idum)))
      ph = 2.0d0 * pi * ran2(idum)
      P_theta(i) = r * cos(ph)
    end do
    do d = 1, dim
      do l = 1, lambda
        do i = 1, N
          do j = 1, N
            y = sqrt(-2.0d0 * log(1.0d0 - ran2(idum)))
            y1 = sqrt(-2.0d0 * log(1.0d0 - ran2(idum)))
            ph_y = 2.0d0 * pi * ran2(idum)
            ph_y1 = 2.0d0 * pi * ran2(idum)
            P_phi(d, l, i, j) = y * cos(ph_y)
            P_phi(d, l, j, i) = y1 * cos(ph_y1)
          end do
        end do
        do i = 1, N
          z = sqrt(-2.0d0 * log(1.0d0 - ran2(idum)))
          ph_z = 2.0d0 * pi * ran2(idum)
          P_phi(d, l, i, i) = z * cos(ph_z)
        end do
      end do
    end do
  end subroutine Generate_momentum

  !------------------------------------------------------------------
  ! Subroutine: theta_constraints
  !  - Checks whether the theta variables satisfy a periodicity constraint.
  !------------------------------------------------------------------
  subroutine theta_constraints(theta, acceptance)
    implicit none
    double precision, intent(in) :: theta(N)
    double precision, intent(out) :: acceptance
    integer :: i
    double precision :: max_val, min_val, pi_val
    max_val = theta(1)
    min_val = theta(1)
    do i = 2, N
      if (theta(i) > max_val) then
        max_val = theta(i)
      end if
      if (theta(i) < min_val) then
        min_val = theta(i)
      end if
    end do
    pi_val = 2.0d0 * asin(1.0d0)
    if ((max_val - min_val) < 2.0d0 * pi_val) then
      acceptance = 0.0d0
    else
      acceptance = 1.0d0
    end if
  end subroutine theta_constraints

  !------------------------------------------------------------------
  ! Subroutine: molecular_dynamics
  !  - Advances the system using a leapfrog integrator.
  !------------------------------------------------------------------
  subroutine molecular_dynamics(idum, P_phi, phi, theta, P_theta, force_phi, force_theta, ngauge, ham)
    implicit none
    integer, intent(in) :: idum, ngauge
    double precision, intent(inout) :: phi(dim, lambda, N, N), theta(N)
    double precision, intent(inout) :: P_phi(dim, lambda, N, N), P_theta(N)
    double precision, intent(inout) :: force_phi(dim, lambda, N, N), force_theta(N)
    double precision :: action_boson
    integer :: d, l, i, j, step

    call hamiltonian(idum, phi, theta, P_phi, P_theta, ngauge, action_boson, ham)
    ! FIRST LEAPFROG STEP: update fields by half-step
    do d = 1, dim
      do l = 1, lambda
        do i = 1, N
          do j = 1, N
            phi(d, l, i, j) = phi(d, l, i, j) + 0.5d0 * dtau_phi * P_phi(d, l, i, j)
          end do
        end do
      end do
    end do
    theta = theta + 0.5d0 * dtau_theta * P_theta

    ! SECOND LEAPFROG STEPS
    do step = 2, ntau
      call force(idum, phi, theta, force_phi, force_theta)
      do d = 1, dim
        do l = 1, lambda
          do i = 1, N
            do j = 1, N
              P_phi(d, l, i, j) = P_phi(d, l, i, j) + 0.5d0 * dtau_phi * force_phi(d, l, i, j)
            end do
          end do
        end do
      end do
      do d = 1, dim
        do l = 1, lambda
          do i = 1, N
            do j = 1, N
              phi(d, l, i, j) = phi(d, l, i, j) + dtau_phi * P_phi(d, l, i, j)
            end do
          end do
        end do
      end do
    end do

    ! LAST LEAPFROG STEP: update momenta and fields by half-step
    call force(idum, phi, theta, force_phi, force_theta)
    P_theta = P_theta - dtau_theta * force_theta
    do d = 1, dim
      do l = 1, lambda
        do i = 1, N
          do j = 1, N
            P_phi(d, l, i, j) = P_phi(d, l, i, j) - dtau_phi * force_phi(d, l, i, j)
          end do
        end do
      end do
    end do
    do d = 1, dim
      do l = 1, lambda
        do i = 1, N
          do j = 1, N
            phi(d, l, i, j) = phi(d, l, i, j) + 0.5d0 * dtau_phi * P_phi(d, l, i, j)
          end do
        end do
      end do
    end do
    theta = theta + 0.5d0 * dtau_theta * P_theta

    call hamiltonian(idum, phi, theta, P_phi, P_theta, ngauge, action_boson, ham)
  end subroutine molecular_dynamics

  !------------------------------------------------------------------
  ! Subroutine: metropolis
  !  - Implements the Metropolis acceptance/rejection step.
  !------------------------------------------------------------------
  subroutine metropolis(idum, phi, theta, P_phi, P_theta, force_phi, force_theta, ngauge, accept, reject, ham)
    implicit none
    integer, intent(in) :: idum, ngauge
    double precision, intent(inout) :: phi(dim, lambda, N, N), theta(N)
    double precision, intent(inout) :: P_phi(dim, lambda, N, N), P_theta(N)
    double precision, intent(inout) :: force_phi(dim, lambda, N, N), force_theta(N)
    double precision, intent(out) :: ham
    double precision :: action_boson, variationH, acceptance, probability, r
    double precision :: backup_phi(dim, lambda, N, N), backup_theta(N)
    
    ! Backup current fields
    backup_phi = phi
    backup_theta = theta

    call Generate_momentum(idum, P_theta, P_phi)
    if (ngauge == 1) then
      P_theta = 0.0d0
    end if

    call hamiltonian(idum, phi, theta, P_phi, P_theta, ngauge, action_boson, ham)
    call molecular_dynamics(idum, P_phi, phi, theta, P_theta, force_phi, force_theta, ngauge, ham)
    call theta_constraints(theta, acceptance)
    variationH = -ham
    call hamiltonian(idum, phi, theta, P_phi, P_theta, ngauge, action_boson, ham)
    variationH = variationH + ham

    if (acceptance == 1.0d0) then
      ! Constraint violation flag; record if necessary
    else
      if (variationH < 0.0d0) then
        accept = accept + 1.0d0
      else
        probability = exp(-variationH)
        r = ran2(idum)
        if (r < probability) then
          accept = accept + 1.0d0
        else
          phi = backup_phi
          theta = backup_theta
          reject = reject + 1.0d0
        end if
      end if
    end if
  end subroutine metropolis

  !------------------------------------------------------------------
  ! Subroutine: observables
  !  - Computes observables; here, the Polyakov loop is used.
  !------------------------------------------------------------------
  subroutine observables(phi, theta, P_phi, P_theta, polyakov)
    implicit none
    double precision, intent(in) :: phi(dim, lambda, N, N), theta(N)
    double precision, intent(in) :: P_phi(dim, lambda, N, N), P_theta(N)
    double precision, intent(out) :: polyakov
    call polyakov_loop(theta, polyakov)
  end subroutine observables

  !------------------------------------------------------------------
  ! Subroutine: polyakov_loop
  !  - Computes the Polyakov loop observable.
  !------------------------------------------------------------------
  subroutine polyakov_loop(theta, polyakov)
    implicit none
    double precision, intent(in) :: theta(N)
    double precision, intent(out) :: polyakov
    integer :: i
    double precision :: re_pol, im_pol
    re_pol = 0.0d0
    im_pol = 0.0d0
    do i = 1, N
      re_pol = re_pol + cos(theta(i))
      im_pol = im_pol + sin(theta(i))
    end do
    re_pol = re_pol / dble(N)
    im_pol = im_pol / dble(N)
    polyakov = sqrt(re_pol**2 + im_pol**2)
  end subroutine polyakov_loop

end module bfss_simulation







     
