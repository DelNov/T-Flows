!==============================================================================!
  subroutine User_Mod_Beginning_Of_Time_Step(Flow, Turb, Vof, Swarm,  &
                                             curr_dt, time)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Swarm_Type),    target :: Swarm
  integer, intent(in)         :: curr_dt  ! time step
  real,    intent(in)         :: time     ! physical time
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: u, v, w, t, phi
  integer                  :: c1, c2, s
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  u    => Flow % u
  v    => Flow % u
  w    => Flow % u

  ! Impose slip condition
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2 < 0) then
      if(Grid % Bnd_Cond_Name_At_Cell(c2) .eq. 'TOP') then
        u % n(c2) = u % n(c1)
        v % n(c2) = v % n(c1)
        w % n(c2) = 0.0
      end if
    end if ! c2 < 0
  end do

  end subroutine

