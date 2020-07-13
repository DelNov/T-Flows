!==============================================================================!
  subroutine Multiphase_Mod_Vof_Mass_Transfer_Rate_In(mult, mass_in)
!------------------------------------------------------------------------------!
!   Computes mass transfer rate due to phase change as an inflow quantity      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  real                          :: mass_in, mass_out
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),      pointer :: grid
  type(Field_Type),     pointer :: flow
  type(Var_Type),       pointer :: vof
  real, contiguous,     pointer :: flux_rate(:)
  integer                       :: s, ss, c, c1, c2, cc1, cc2
  real                          :: epsloc
!==============================================================================!

  grid => mult % pnt_grid
  flow => mult % pnt_flow
  vof  => mult % vof

  epsloc = epsilon(epsloc)

  flux_rate  => mult % flux_rate

  mult % add_mass_in = 0.0

  ! The total added mass
  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    mult % add_mass_in = mult % add_mass_in + flux_rate(c) * grid % vol(c)     &
                                            * flow % density(c)                &
                                            * ( 1.0 / mult % phase_dens(2)     &
                                              - 1.0 / mult % phase_dens(1) )
  end do

  call Comm_Mod_Global_Sum_Real(mult % add_mass_in)

  mass_in  = mass_in + mult % add_mass_in

end subroutine
