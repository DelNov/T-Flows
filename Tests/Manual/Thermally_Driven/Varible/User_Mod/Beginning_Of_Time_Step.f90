!==============================================================================!
  subroutine User_Mod_Beginning_Of_Time_Step(Flow, Turb, Vof, Swarm)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of time step.                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Swarm_Type),    target :: Swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: u, v, w, t, phi
  integer                  :: c, i
  real                     :: wi, wip
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid

  !------------------------------!
  !   Browse through all cells   !
  !------------------------------!
  do c = Cells_At_Boundaries_In_Domain_And_Buffers()

    ! Browse through all table entries
    do i = 1, N_ITEMS - 1

      ! Did you find the right interval
      if(Flow % t % n(c) >= air_t(i) .and.  &
         Flow % t % n(c) <  air_t(i+1)) then

        ! If so, calculate interpolation factors ...
        wi  = (air_t(i+1) - Flow % t % n(c)) / (air_t(i+1) - air_t(i))
        wip = 1.0 - wi

        ! ... and interpolate physical properties
        Flow % density(c)      = wi * air_rho   (i)  + wip * air_rho   (i+1)
        Flow % viscosity(c)    = wi * air_mu    (i)  + wip * air_mu    (i+1)
        Flow % conductivity(c) = wi * air_lambda(i)  + wip * air_lambda(i+1)
        Flow % capacity(c)     = wi * air_cp    (i)  + wip * air_cp    (i+1)
      end if

    end do
  end do

  end subroutine

