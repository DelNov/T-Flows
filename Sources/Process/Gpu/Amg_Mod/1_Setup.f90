!==============================================================================!
  subroutine Setup(Amg, nnu, levelx,  &
                   a, u, ia, ja,      &
                   iw, icg, ifg,      &  ! these are "created" from kwork
                   levels,            &
                   iwork, jtr)
!------------------------------------------------------------------------------!
!   Preparation phase of Amg1r5 (general part)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type) :: Amg
  integer         :: nnu, levelx
  real            :: a(:), u(:)
  integer         :: ia(:), ja(:)
  integer         :: iw(:), icg(:), ifg(:)
  integer         :: levels
  integer         :: iwork(:), jtr(:)
!-----------------------------------[locals]-----------------------------------!
  integer :: i, n_digits
  integer :: digit(AMG_MAX_LEVELS)
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  !------------------------!
  !   Decompose "matrix"   !
  !------------------------!
  call Amg % Get_Integer_Digits(Amg % matrix, 2, n_digits, digit)
  Amg % isym  = digit(1)
  Amg % irow0 = digit(2)

  !-------------------------------------!
  !   Preparation (ignored in timing)   !
  !-------------------------------------!
  Amg % imin(1) = 1
  Amg % imax(1) = nnu
  call Amg % Check_Matrix_Properties(a, ia, ja, icg, ifg)
  if(Amg % ierr .gt. 0) return

  !-----------------------------------------------------------------!
  !   Define coarser grids + operators. reset levels if necessary   !
  !-----------------------------------------------------------------!
  call Amg % Coarsening(levelx,        &
                        a, u, ia, ja,  &  ! store linear system
                        iw, icg, ifg,  &  ! work arrays
                        levels,        &
                        iwork, jtr)

  end subroutine
