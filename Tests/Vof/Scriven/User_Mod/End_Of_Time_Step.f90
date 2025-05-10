!==============================================================================!
  subroutine User_Mod_End_Of_Time_Step(Flow, Turb, Vof, Swarm,  &
                                       n_stat_t, n_stat_p)
!------------------------------------------------------------------------------!
!   This function is computing benchmark for rising bubble.                    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: Turb
  type(Vof_Type),   target :: Vof
  type(Swarm_Type), target :: Swarm
  integer                  :: n_stat_t  ! 1st t.s. statistics turbulence
  integer                  :: n_stat_p  ! 1st t.s. statistics particles
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: fun
  type(Front_Type), pointer:: Front
  type(Elem_Type), pointer :: Elem(:)
  integer                  :: c, e, fu
  real                     :: b_volume, surface, ne, rad_volume 
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  fun  => Vof % fun
  Front => Vof % Front
  Elem => Front % Elem

  !-------------------!
  !   Bubble volume   !
  !-------------------!
  b_volume = 0.0
  surface = 0.0
  ne = 0
!  alpha_l = Vof%phase_cond(1)/(Vof%phase_capa(1)*Vof%phase_dens(1))

! do c = Cells_In_Domain()
  do c = Cells_In_Domain()
    b_volume = b_volume + Grid % vol(c) * (1.0 - fun % n(c))
    e  = Vof % Front % elem_in_cell(c)
    if(e > 0) then
      ne = ne + 1
      surface = surface + Elem(e) % area
    end if
  end do

  call Global % Sum_Real(ne)
  call Global % Sum_Real(surface)
  call Global % Sum_Real(b_volume)
  rad_volume = (3.0*b_volume/(4*PI))**(1.0/3.0)

  ! Write to file
  if(First_Proc()) then
    call File % Append_For_Writing_Ascii('bench-data.dat', fu)
    write(fu,'(4(2X,E16.8))') Time % Get_Time(), surface, b_volume, rad_volume
    close(fu)
  end if

  end subroutine
