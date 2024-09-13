!==============================================================================!
  subroutine Insert_Energy_Bc(Process, Grid, Flow)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Process
  type(Grid_Type),  target :: Grid
  type(Field_Type), target :: Flow
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: b(:), fc(:), cond(:)
  real                      :: a12
  integer                   :: reg, s, c, c1, c2
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Insert_Energy_Bc')

  ! Take some aliases
  b    => Flow % Nat % b
  fc   => Flow % Nat % C % fc
  cond => Flow % conductivity

  !-----------------------------------------------------------------------!
  !   Handle boundary conditions on the right-hand side (in the source)   !
  !-----------------------------------------------------------------------!

  !$tf-acc loop begin
  do c = Cells_In_Domain_And_Buffers()  ! all present
    b(c) = 0.0
  end do
  !$tf-acc loop end

  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. WALL .or.  &
       Grid % region % type(reg) .eq. INFLOW) then

      !$tf-acc loop begin
      do s = Faces_In_Region(reg)  ! all present
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)
        a12 = cond(c1) * fc(s)
        b(c1) = b(c1) + a12 * Flow % t % n(c2)
      end do
      !$tf-acc loop end

    end if
  end do

  call Profiler % Stop('Insert_Energy_Bc')

  end subroutine
