!==============================================================================!
  real function Root_Mean_Square(ni, r)
!------------------------------------------------------------------------------!
!   Calculates root means square of vector r without normalization.            !
!   This non-normalized variant seems to be better option for ACM.             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Comm_Mod
  use Matrix_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer :: ni
  real    :: r(:)  ! this may be only in inner cells
!-----------------------------------[Locals]-----------------------------------!
  real    :: rms
  integer :: i
!==============================================================================!

  ! Compute rms normalizing it with main diagonal in the system matrix
  rms = 0.0
  do i = 1, ni
    rms = rms + r(i)**2
  end do
  call Comm_Mod_Global_Sum_Real(rms)
  rms = sqrt(rms)

  Root_Mean_Square = rms

  end function
