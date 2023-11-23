!==============================================================================!
  subroutine Smooths_Mod_Allocate_Smooths(smr, n)
!------------------------------------------------------------------------------!
!>  Allocates arrays for Smooth_Type which are used in grid smoothing.
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Smooths_Type) :: smr  !! smoothing regions
  integer            :: n    !! number of smoothing regions
!==============================================================================!

  smr % n_smooths = n     ! number of smoothing regions

  allocate(smr % iters(n));  smr % iters = 0
  allocate(smr % in_x (n));  smr % in_x  = .false.
  allocate(smr % in_y (n));  smr % in_y  = .false.
  allocate(smr % in_z (n));  smr % in_z  = .false.
  allocate(smr % x_min(n));  smr % x_min = +HUGE
  allocate(smr % y_min(n));  smr % y_min = +HUGE
  allocate(smr % z_min(n));  smr % z_min = +HUGE
  allocate(smr % x_max(n));  smr % x_max = -HUGE
  allocate(smr % y_max(n));  smr % y_max = -HUGE
  allocate(smr % z_max(n));  smr % z_max = -HUGE
  allocate(smr % relax(n));  smr % relax = 1.0   ! is this a smart choice?

  end subroutine
