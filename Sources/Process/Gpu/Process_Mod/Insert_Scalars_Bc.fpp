!==============================================================================!
  subroutine Insert_Scalars_Bc(Process, Grid, Flow, sc)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Process
  type(Grid_Type),  target :: Grid
  type(Field_Type), target :: Flow
  integer, intent(in)      :: sc
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type),   pointer :: phi
  real, contiguous, pointer :: b(:), fc(:), diff(:)
  real                      :: a12
  integer                   :: reg, s, c, c1, c2
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Insert_Scalars_Bc')

  ! Take some aliases
  b    => Flow % Nat % b
  fc   => Flow % Nat % A % fc
  diff => Flow % diffusivity
  phi  => Flow % scalar(sc)

  !-----------------------------------------------------------------------!
  !   Handle boundary conditions on the right-hand side (in the source)   !
  !-----------------------------------------------------------------------!

  !$tf-acc loop begin
  do c = Cells_In_Domain_And_Buffers()  ! all present
    b(c) = 0.0
  end do
  !$tf-acc loop end

  !$tf-acc loop begin
  do s = Faces_At_Boundaries()  ! all present
    c1 = Grid % faces_c(1,s)    ! inside cell
    c2 = Grid % faces_c(2,s)    ! boundary cell
    if(phi % bnd_cond_type(c2) .eq. WALL    .or.  &
       phi % bnd_cond_type(c2) .eq. INFLOW) then
      a12 = diff(c1) * fc(s)
      b(c1) = b(c1) + a12 * phi % n(c2)
    end if
  end do
  !$tf-acc loop end

  call Profiler % Stop('Insert_Scalars_Bc')

  end subroutine
