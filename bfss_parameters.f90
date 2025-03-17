!********************************************************************
! Module: bfss_parameters
!  - Contains simulation constants and parameters.
!  - Uses iso_fortran_env to define a double‐precision kind.
!********************************************************************
module bfss_parameters
  use iso_fortran_env, only: dp => real64
  implicit none

  integer, parameter :: dim     = 9
  integer, parameter :: lambda  = 10
  integer, parameter :: N       = 10
  integer, parameter :: jsweeps = 100
  integer, parameter :: ntau    = 2.0

  real(dp), parameter :: a           = 2.0_dp
  real(dp), parameter :: mass        = 0.005_dp
  real(dp), parameter :: coupling    = 1.0_dp
  real(dp), parameter :: temperature = 1.0_dp

  real(dp), parameter :: dtau_phi   = 1.0_dp
  real(dp), parameter :: dtau_theta = 1.0_dp

  integer, parameter :: ngauge_default = 0

  ! Normalization factor = 1 / (dim*lambda*N*N)
  real(dp), parameter :: norm_factor = 1.0_dp / dble(dim*lambda*N*N)
end module bfss_parameters
