!==============================================================================!
  subroutine User_Mod_Beginning_Of_Correct_Velocity(Flow, Vof, sol,  &
                                                    curr_dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the end of Correct_Velocity function.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),   target :: Flow
  type(Vof_Type),     target :: Vof
  type(Solver_Type),  target :: sol
  integer, intent(in)        :: curr_dt
  integer, intent(in)        :: ini
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Var_Type),    pointer :: u
  type(Face_Type),   pointer :: v_flux          ! volume flux
  integer                    :: c, s
  character(40)              :: file_name =  &
                                'face_vel_before_correction_xxxxx_yyy.dat'
!==============================================================================!

  ! Take aliases
  grid   => Flow % pnt_grid
  v_flux => Flow % v_flux
  u      => Flow % u

  !-----------------------------------------------------------!
  !   Write down face-centered velocities before correction   !
  !-----------------------------------------------------------!
  write(file_name(28:32), '(i5.5)') curr_dt
  write(file_name(34:36), '(i3.3)') ini

  open(99, file=file_name)
  write(99, '(a)') '# User_Mod_Beginning_Of_Correct_Velocity: '  //  &  
                   ' 1:x  2:u @ face  3:v_flux @ face'
  do s = 1, grid % n_faces
    if(Math_Mod_Approx_Real(grid % yf(s), 0.0) .and.  &
       Math_Mod_Approx_Real(grid % zf(s), 0.0) .and.  &
       grid % faces_c(2,s) > 0) then
      write(99, '(99es18.7)')  &
        grid % xf(s), v_flux % n(s) / grid % s(s), v_flux % n(s)
    end if
  end do
  close(99)

  end subroutine