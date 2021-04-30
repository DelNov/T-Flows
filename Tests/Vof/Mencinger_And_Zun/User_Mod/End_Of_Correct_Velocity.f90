!==============================================================================!
  subroutine User_Mod_End_Of_Correct_Velocity(flow, mult, sol, curr_dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the end of Correct_Velocity function.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  integer, intent(in)           :: curr_dt
  integer, intent(in)           :: ini
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u
  type(Face_Type),   pointer :: v_flux          ! volume flux
  integer                    :: c, s
  character(39)              :: fil1_name =  &
                                'velocity_after_correction_xxxxx_yyy.dat'
  character(39)              :: fil2_name =  &
                                'face_vel_after_correction_xxxxx_yyy.dat'
!==============================================================================!

  ! Take aliases
  grid   => flow % pnt_grid
  v_flux => flow % v_flux
  u      => flow % u

  !---------------------------------------------------!
  !   Write down corrected cell-centered velocities   !
  !---------------------------------------------------!
  write(fil1_name(27:31), '(i5.5)') curr_dt
  write(fil1_name(33:35), '(i3.3)') ini

  open(99, file=fil1_name)
  write(99, '(a)') '# User_Mod_End_Of_Correct_Velocity: 1:x  2:u  3:rho'

  do c = 1, grid % n_cells
    if(Math_Mod_Approx_Real(grid % yc(c), 0.0) .and.  &
       Math_Mod_Approx_Real(grid % zc(c), 0.0)) then
      write(99, '(99es15.5)')  &
        grid % xc(c), u % n(c), flow % density(c)
    end if
  end do

  close(99)

  !---------------------------------------------------!
  !   Write down corrected face-centered velocities   !
  !---------------------------------------------------!
  write(fil2_name(27:31), '(i5.5)') curr_dt
  write(fil2_name(33:35), '(i3.3)') ini

  open(99, file=fil2_name)
  write(99, '(a)') '# User_Mod_End_Of_Correct_Velocity: '  //  &  
                   ' 1:x  2:u @ face  3:v_flux'
  do s = 1, grid % n_faces
    if(Math_Mod_Approx_Real(grid % yf(s), 0.0) .and.  &
       Math_Mod_Approx_Real(grid % zf(s), 0.0) .and.  &
       grid % faces_c(2,s) > 0) then
      write(99, '(99es15.5)')  &
        grid % xf(s), v_flux % n(s) / grid % s(s), v_flux % n(s)
    end if
  end do
  close(99)

  end subroutine
