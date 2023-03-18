!==============================================================================!
  real function Root_Mean_Square(Native, ni, r)
!------------------------------------------------------------------------------!
!   Calculates root means square of vector r without normalization.            !
!   This non-normalized variant seems to be better option for ACM.             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Native_Type), intent(in) :: Native
  integer,            intent(in) :: ni
  real,               intent(in) :: r(:)  ! this may be only in inner cells
!-----------------------------------[Locals]-----------------------------------!
  real    :: rms
  integer :: i
!==============================================================================!

  ! Compute rms normalizing it with main diagonal in the system matrix
  rms = 0.0
  !$omp parallel do private(i) shared(r) reduction(+ : rms)
  do i = 1, ni
    rms = rms + r(i)**2
  end do
  !$omp end parallel do
  call Comm_Mod_Global_Sum_Real(rms)
  rms = sqrt(rms)

  Root_Mean_Square = rms

  end function
