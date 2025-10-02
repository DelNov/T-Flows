!==============================================================================!
  subroutine User_Mod_Beginning_Of_Simulation(Flow, Turb, Vof, Swarm)
!------------------------------------------------------------------------------!
!   This function is called at the beginning of simulation.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Swarm_Type),    target :: Swarm
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: i, fu
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid

  !----------------------------------------!
  !   Open file with physical properties   !
  !----------------------------------------!
  call File % Open_For_Reading_Ascii("air_properties_at_1_bar.dat", fu)

  !-----------------------------!
  !   Read all the properties   !
  !-----------------------------!
  do i = 1, N_ITEMS
    call File % Read_Line(fu)

    ! Read desired properties
    read(line % tokens(1), *) air_t(i)
    read(line % tokens(2), *) air_rho(i)
    read(line % tokens(3), *) air_mu(i)
    read(line % tokens(5), *) air_cp(i)
    read(line % tokens(6), *) air_lambda(i)

    ! Fix units where needed (check the values in the table)
    air_cp(i) = air_cp(i) * 1.0e3
    air_mu(i) = air_mu(i) / 1.0e5
  end do

  close(fu)

  if(First_Proc()) then
    print '(a)',        ' #============================================'
    print '(a)',        ' # Output from user function, read properties!'
    print '(a)',        ' #--------------------------------------------'
  end if

  end subroutine
