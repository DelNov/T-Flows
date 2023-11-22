!==============================================================================!
  subroutine User_Mod_Beginning_Of_Time_Step(Flow, Turb, Vof, Swarm)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
!------------------------------[Local parameters]------------------------------!
  real, parameter :: U_BULK = 1.0           ! bulk velocity
  real, parameter :: U_MAX  = 1.5 * U_BULK  ! max velocity
  real, parameter :: HALF_W = 2.05          ! half channel width
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: u, v, w, t, phi
  integer                  :: reg, c2
  real                     :: y
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  u    => Flow % u

  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. INFLOW) then
      do c2 = Cells_In_Region(reg)
        y = Grid % yc(c2)
        u % n(c2) = U_MAX * (1.0 - ( (y-HALF_W)/HALF_W )**2)
      end do  ! through boundary cells of the region
    end if    ! if at inflow region
  end do      ! through all boundary regions

  end subroutine

