!==============================================================================!
  subroutine Turb_Mod_Vis_T_Hybrid_Les_Rans(turb)
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
  type(Turb_Type), target :: turb
!-----------------------------------[Locals]-----------------------------------!
  type(Field_Type), pointer :: flow
  type(Grid_Type),  pointer :: grid
  type(Var_Type),   pointer :: t
  integer                   :: c
  real                      :: lf, pr_t, nc2
!==============================================================================!

  ! Take aliases
  flow => turb % pnt_flow
  grid => flow % pnt_grid
  t    => flow % t


  if(buoyancy) then
    call Field_Mod_Grad_Variable(flow, t)
  end if

  do c = 1, grid % n_cells
    lf = grid % vol(c) ** ONE_THIRD
    turb % vis_t_sgs(c) = flow % density(c)  &
                        * (lf*lf)            &          ! delta^2 
                        * turb % c_dyn(c)    &          ! c_dynamic   
                        * flow % shear(c)
  end do

  if(buoyancy) then
    do c = 1, grid % n_cells
      nc2 = -flow % beta * (  grav_x * t % x(c)   &   
                            + grav_y * t % y(c)   &   
                            + grav_z * t % z(c))  
      turb % vis_t_sgs(c) = turb % vis_t_sgs(c)  &
             * max((1.0 - 2.5 * nc2 / (flow % shear(c) + TINY)), 0.0)
    end do    
  end if

  call Grid_Mod_Exchange_Cells_Real(grid, turb % vis_t_sgs)

  end subroutine
