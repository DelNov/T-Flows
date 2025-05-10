!==============================================================================!
  subroutine User_Mod_End_Of_Correct_Velocity(Flow, Vof, Sol)
!------------------------------------------------------------------------------!
!   This function is called at the end of Correct_Velocity function.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Vof_Type),    target :: Vof
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: u
  type(Face_Type),   pointer :: v_flux          ! volume flux
  integer                    :: c, s
  character(39)              :: fil1_name =  &
                                'velocity_after_correction_xxxxx_yyy.dat'
  character(39)              :: fil2_name =  &
                                'face_vel_after_correction_xxxxx_yyy.dat'
!==============================================================================!

  ! Take aliases
  Grid   => Flow % pnt_grid
  v_flux => Flow % v_flux
  u      => Flow % u

  !---------------------------------------------------!
  !   Write down corrected cell-centered velocities   !
  !---------------------------------------------------!
  write(fil1_name(27:31), '(i5.5)') Time % Curr_Dt()
  write(fil1_name(33:35), '(i3.3)') Iter % Current()

  open(99, file=fil1_name)
  write(99, '(a)') '# User_Mod_End_Of_Correct_Velocity: 1:x  2:u  3:rho'

  do c = Cells_In_Domain_And_Buffers()
    if(Math % Approx_Real(Grid % yc(c), 0.0) .and.  &
       Math % Approx_Real(Grid % zc(c), 0.0)) then
      write(99, '(99es15.5)')  &
        Grid % xc(c), u % n(c), Flow % density(c)
    end if
  end do

  close(99)

  !---------------------------------------------------!
  !   Write down corrected face-centered velocities   !
  !---------------------------------------------------!
  write(fil2_name(27:31), '(i5.5)') Time % Curr_Dt()
  write(fil2_name(33:35), '(i3.3)') Iter % Current()

  open(99, file=fil2_name)
  write(99, '(a)') '# User_Mod_End_Of_Correct_Velocity: '  //  &  
                   ' 1:x  2:u @ face  3:v_flux'
  do s = 1, Grid % n_faces
    if(Math % Approx_Real(Grid % yf(s), 0.0) .and.  &
       Math % Approx_Real(Grid % zf(s), 0.0) .and.  &
       Grid % faces_c(2,s) > 0) then
      write(99, '(99es25.15)')  &
        Grid % xf(s), v_flux % n(s) / Grid % s(s), v_flux % n(s)
    end if
  end do
  close(99)

  end subroutine
