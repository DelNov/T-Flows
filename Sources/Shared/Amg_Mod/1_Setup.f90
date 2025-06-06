!==============================================================================!
  subroutine setup(amg, nnu, levelx,  &
                   a, u, ia, ja,      &
                   iw, icg, ifg,      &  ! these are "created" from kwork
                   levels,            &
                   iwork, jtr)
!------------------------------------------------------------------------------!
!   Preparation phase of amg1r5 (general part)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(amg_type)  :: amg
  integer          :: nnu, levelx
  double precision :: a(:), u(:)
  integer          :: ia(:), ja(:)
  integer          :: iw(:), icg(:), ifg(:)
  integer          :: levels
  integer          :: iwork(:), jtr(:)
!-----------------------------------[locals]-----------------------------------!
  integer :: i, n_digits
  integer :: digit(AMG_MAX_LEVELS)
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  !------------------------!
  !   Decompose "matrix"   !
  !------------------------!
  call amg % get_integer_digits(amg % matrix, 2, n_digits, digit)
  amg % isym  = digit(1)
  amg % irow0 = digit(2)

  !-------------------------!
  !   Reset time counters   !
  !-------------------------!
  do i = 1, 20
    amg % time(i) = 0.0
  end do

  !-------------------------------------!
  !   Preparation (ignored in timing)   !
  !-------------------------------------!
  amg % imin(1) = 1
  amg % imax(1) = nnu
  call amg % check_matrix_properties(a, ia, ja, icg, ifg)
  if(amg % ierr .gt. 0) return

  !-----------------------------------------------------------------!
  !   Define coarser grids + operators. reset levels if necessary   !
  !-----------------------------------------------------------------!
  call amg % coarsening(levelx,        &
                        a, u, ia, ja,  &  ! store linear system
                        iw, icg, ifg,  &  ! work arrays
                        levels,        &
                        iwork, jtr)

  end subroutine
