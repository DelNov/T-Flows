!==============================================================================!
  subroutine Vis_T_Hybrid_Les_Rans(Turb)
!------------------------------------------------------------------------------!
!   Calculates SGS stresses and turbulent viscosity for 'LES'.                 !
!------------------------------------------------------------------------------!
!   near(c) is the number of corresponding cell on the nearest wall.
!   In case that, in parallel executions, the subdomain does not have
!   any nearwall cells, the near(c) is zero.
!   near(c) is calculated in NearWallCells.f90, only ones in the beginig
!   of a simulation.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Turb_Type), target :: Turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: t
  integer                   :: c
  real                      :: lf, nc2
!==============================================================================!

  ! Take aliases
  Flow => Turb % pnt_flow
  Grid => Flow % pnt_grid
  t    => Flow % t

  if(Flow % buoyancy .eq. THERMALLY_DRIVEN) then
    call Flow % Grad_Variable(t)
  end if

  do c = Cells_In_Domain_And_Buffers()
    lf = Grid % vol(c) ** ONE_THIRD
    Turb % vis_t_sgs(c) = Flow % density(c)  &
                        * (lf*lf)            &          ! delta^2
                        * Turb % c_dyn(c)    &          ! c_dynamic
                        * Flow % shear(c)
  end do

  if(Flow % buoyancy .eq. THERMALLY_DRIVEN) then
    do c = Cells_In_Domain_And_Buffers()
      nc2 = -Flow % beta * (  Flow % grav_x * t % x(c)   &
                            + Flow % grav_y * t % y(c)   &
                            + Flow % grav_z * t % z(c))
      Turb % vis_t_sgs(c) = Turb % vis_t_sgs(c)  &
             * max((1.0 - 2.5 * nc2 / (Flow % shear(c) + TINY)), 0.0)
    end do
  end if

  call Grid % Exchange_Cells_Real(Turb % vis_t_sgs)

  end subroutine
