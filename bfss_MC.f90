module bfss_MC
  use bfss_parameters
  use bfss_random
  implicit none

contains

  !------------------------------------------------------------------
  ! Subroutine: action
  !  - Computes the bosonic action including kinetic and Faddeev-Popov terms.
  !  - Modifies 'theta' by assigning random values.
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
        l1 = mod(l, lambda) + 1  ! cyclic index: if l==lambda then l1=1
        do i = 1, N
          y = ran2(idum)
          theta(i) = y
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
            if (sin(angle)**2 > 1.0E-12) then
              S_FP = S_FP - 0.5d0 * log(sin(angle)**2)
            end if
          end do
        end do
      end do
    end do

    action_boson =norm_factor*( S_FP - (N/a) * kinetic + N * ((1.0d0/a)*rr + 0.5d0*a*(mass**2)*rr))
  end subroutine action

  !------------------------------------------------------------------
  ! Subroutine: force
  !  - Computes the forces (gradients) with respect to phi and theta.
  !  - 'phi' and 'theta' are read-only; force_phi and force_theta are computed.
  !  - Changed force_phi and force_theta to INTENT(INOUT) to allow read-modify.
  !------------------------------------------------------------------
  subroutine force(idum, phi, theta, force_phi, force_theta)
    implicit none
    integer, intent(inout) :: idum
    double precision, intent(in) :: phi(dim, lambda, N, N), theta(N)
    double precision, intent(inout) :: force_phi(dim, lambda, N, N), force_theta(N)
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
  !  - Computes the Hamiltonian as the sum of bosonic action and kinetic terms.
  !  - Calls 'action' which modifies 'theta' and 'idum'.
  !------------------------------------------------------------------
  subroutine hamiltonian(idum, phi, theta, P_phi, P_theta, ngauge, action_boson, ham)
    implicit none
    integer, intent(inout) :: idum
    integer, intent(in)    :: ngauge
    double precision, intent(in) :: phi(dim, lambda, N, N)
    double precision, intent(inout) :: theta(N)
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
    ham = norm_factor * ham
  end subroutine hamiltonian

  !------------------------------------------------------------------
  ! Subroutine: Generate_momentum
  !  - Generates momentum variables with a Gaussian distribution.
  !------------------------------------------------------------------
  subroutine Generate_momentum(idum, P_theta, P_phi)
    implicit none
    integer, intent(inout) :: idum
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
            P_phi(d, l, i, j) = y  * cos(ph_y)
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
      if (theta(i) > max_val) max_val = theta(i)
      if (theta(i) < min_val) min_val = theta(i)
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
  !  - Calls 'hamiltonian' (which modifies 'theta', 'idum') and 'force'.
  !------------------------------------------------------------------
  subroutine molecular_dynamics(idum, P_phi, phi, theta, P_theta, force_phi, force_theta, ngauge, ham)
    implicit none
    integer, intent(inout) :: idum
    integer, intent(in)    :: ngauge
    double precision, intent(inout) :: phi(dim, lambda, N, N), theta(N)
    double precision, intent(inout) :: P_phi(dim, lambda, N, N), P_theta(N)
    double precision, intent(inout) :: force_phi(dim, lambda, N, N), force_theta(N)
    double precision :: action_boson, ham
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
  !  - Calls 'hamiltonian' (which modifies 'theta', 'idum') and 'molecular_dynamics'.
  !------------------------------------------------------------------
  subroutine metropolis(idum, phi, theta, P_phi, P_theta, force_phi, force_theta, ngauge, &
                        accept, reject, ham)
    implicit none
    integer, intent(inout) :: idum
    integer, intent(in)    :: ngauge
    double precision, intent(inout) :: phi(dim, lambda, N, N), theta(N)
    double precision, intent(inout) :: P_phi(dim, lambda, N, N), P_theta(N)
    double precision, intent(inout) :: force_phi(dim, lambda, N, N), force_theta(N)
    double precision, intent(inout) :: accept, reject
    double precision, intent(out) :: ham

    double precision :: action_boson, variationH, acceptance, probability, r
    double precision :: backup_phi(dim, lambda, N, N), backup_theta(N)

    ! Backup current fields
    backup_phi = phi
    backup_theta = theta

    ! Generate new momenta
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
      ! Constraint violation: do nothing special (could record if needed)
    else
      if (variationH < 0.0d0) then
        accept = accept + 1.0d0
      else
        probability = exp(-variationH)
        r = ran2(idum)
        if (r < probability) then
          accept = accept + 1.0d0
        else
          ! Revert to old fields
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

end module bfss_MC

