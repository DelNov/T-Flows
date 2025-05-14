!==============================================================================!
  subroutine Insert_Variable_Bc(Turb, Grid, Flow, phi, vis_eff)
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Turb_Type)         :: Turb
  type(Grid_Type),  target :: Grid
  type(Field_Type), target :: Flow
  type(Var_Type),   target :: phi
  real                     :: vis_eff(-Grid % n_bnd_cells:Grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: b(:), fc(:)
  real                      :: a12
  integer                   :: reg, s, c, c1, c2
!==============================================================================!

  call Profiler % Start('Insert_Variable_Bc')

  ! Take some aliases
  b    => Flow % Nat % b
  fc   => Flow % Nat % A % fc

  !-----------------------------------------------------------------------!
  !   Handle boundary conditions on the right-hand side (in the source)   !
  !-----------------------------------------------------------------------!

  !$tf-acc loop begin
  do c = Cells_In_Domain_And_Buffers()  ! all present
    b(c) = 0.0
  end do
  !$tf-acc loop end

  if(phi % name .eq. 'VIS') then
    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. WALL   .or.  &
         Grid % region % type(reg) .eq. WALLFL .or.  &
         Grid % region % type(reg) .eq. INFLOW) then

        !$tf-acc loop begin
        do s = Faces_In_Region(reg)  ! all present
          c1 = Grid % faces_c(1,s)
          c2 = Grid % faces_c(2,s)
          a12 = vis_eff(c1) * fc(s)
          b(c1) = b(c1) + a12 * phi % n(c2)
        end do
        !$tf-acc loop end

      end if
    end do
  end if

  call Profiler % Stop('Insert_Variable_Bc')

  end subroutine
