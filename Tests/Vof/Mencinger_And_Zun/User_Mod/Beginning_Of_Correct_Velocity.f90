!==============================================================================!
  subroutine User_Mod_Beginning_Of_Correct_Velocity(Flow, Vof, Sol)
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
  character(40)              :: file_name =  &
                                'face_vel_before_correction_xxxxx_yyy.dat'
!==============================================================================!

  ! Take aliases
  Grid   => Flow % pnt_grid
  v_flux => Flow % v_flux
  u      => Flow % u

  !-----------------------------------------------------------!
  !   Write down face-centered velocities before correction   !
  !-----------------------------------------------------------!
  write(file_name(28:32), '(i5.5)') Time % Curr_Dt()
  write(file_name(34:36), '(i3.3)') Iter % Current()

  open(99, file=file_name)
  write(99, '(a)') '# User_Mod_Beginning_Of_Correct_Velocity: '  //  &  
                   ' 1:x  2:u @ face  3:v_flux @ face'
  do s = 1, Grid % n_faces
    if(Math % Approx_Real(Grid % yf(s), 0.0) .and.  &
       Math % Approx_Real(Grid % zf(s), 0.0) .and.  &
       Grid % faces_c(2,s) > 0) then
      write(99, '(99es18.7)')  &
        Grid % xf(s), v_flux % n(s) / Grid % s(s), v_flux % n(s)
    end if
  end do
  close(99)

  end subroutine
