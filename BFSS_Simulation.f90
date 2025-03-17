!********************************************************************
!  BFSS_Simulation.f90
!
!  A consolidated Fortran code that simulates the BFSS (Banks-Fischler-
!  Shenker-Susskind) model. All modules, subroutines, and the main
!  program are in one file for easier compilation.
!********************************************************************

program BFSS_Simulation
  use bfss_parameters
  use bfss_MC
  implicit none

  integer :: isweep, idum, ngauge
  double precision :: ham, accept, reject, pol
  double precision, allocatable :: phi(:,:,:,:), P_phi(:,:,:,:)
  double precision, allocatable :: force_phi(:,:,:,:)
  double precision, allocatable :: theta(:), P_theta(:)
  double precision, allocatable :: force_theta(:)

  ! Allocate arrays
  allocate(phi(dim, lambda, N, N))
  allocate(P_phi(dim, lambda, N, N))
  allocate(force_phi(dim, lambda, N, N))
  allocate(theta(N))
  allocate(P_theta(N))
  allocate(force_theta(N))

  ! Initialize
  phi       = 0.0_dp
  P_phi     = 0.0_dp
  force_phi = 0.0_dp
  theta     = 0.0_dp
  P_theta   = 0.0_dp

  idum   = 5000
  ngauge = ngauge_default
  accept = 0.0_dp
  reject = 0.0_dp

  ! Example: Logging to a file
  open(unit=10, file="simulation_results.dat", status="replace", action="write")
  write(10, '(A)') "Sweep   Hamiltonian      Acceptance      Rejection      PolyakovLoop"

  do isweep = 1, jsweeps
     call metropolis(idum, phi, theta, P_phi, P_theta, force_phi, force_theta, ngauge, &
                     accept, reject, ham)

     ! Compute Polyakov loop
     call observables(phi, theta, P_phi, P_theta, pol)

     ! Print to console
     print*, "Sweep ", isweep, ": H = ", ham*1E-10, "  pol = ", pol

     ! Log to file
     write(10, '(I5,2X,F12.6,2X,F12.6,2X,F12.6,2X,F12.6)') isweep, ham*1E-10, accept, reject, pol
  end do

  close(10)

  ! Deallocate
  deallocate(phi, P_phi, force_phi, theta, P_theta)

end program BFSS_Simulation
