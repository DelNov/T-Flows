!==============================================================================!
  subroutine User_Mod_End_Of_Compute_Energy(Flow, Turb, Vof, Sol)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Energy function.             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Turb_Type),   target :: Turb
  type(Vof_Type),    target :: Vof
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Front_Type), pointer :: Front
  integer                   :: c1, c2, s, fui, fut
  character(33)             :: profile_name
  logical, save             :: first_entry = .true.
!==============================================================================!

  ! Take aliases
  Grid  => Flow % pnt_grid
  Front => Vof % Front

  !----------------------------------------!
  !   Write out temperature distribution   !
  !----------------------------------------!
  profile_name = 'temperature_distribution_0000.num'
  write(profile_name(26:29), '(i4.4)')  Time % Curr_Dt()
  call File % Open_For_Writing_Ascii(profile_name, fut)

  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)
    if(c2 > 0) then
      if(Math % Approx_Real(Grid % yf(s), 0.0) .and.  &
         Math % Approx_Real(Grid % zf(s), 0.0)) then
        write(fut, '(99(es12.4))') Grid % xc(c1), Flow % t % n(c1)
      end if
      if(Front % intersects_face(s)) then
        if(Math % Approx_Real(Front % ys(s), 0.0) .and.  &
           Math % Approx_Real(Front % zs(s), 0.0)) then
          write(fut, '(99(es12.4))') Front % xs(s), Vof % t_sat
        end if
      end if
    end if
  end do

  close(fut)

  !-------------------------------------!
  !   Write out temperature gradients   !
  !-------------------------------------!
  profile_name = 'temperature_gradients_0000.num'
  write(profile_name(23:26), '(i4.4)')  Time % Curr_Dt()
  call File % Open_For_Writing_Ascii(profile_name, fut)

  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1, s)
    c2 = Grid % faces_c(2, s)
    if(c2 > 0) then
      if(Math % Approx_Real(Grid % yf(s), 0.0) .and.  &
         Math % Approx_Real(Grid % zf(s), 0.0)) then
        write(fut, '(99(es12.4))') Grid % xc(c1), Flow % t % x(c1)
      end if
    end if
  end do

  close(fut)

  !----------------------------------!
  !   Write out interface position   !
  !----------------------------------!
  if(first_entry) then
    call File % Open_For_Writing_Ascii('interface_position.num', fui)
  else
    call File % Append_For_Writing_Ascii('interface_position.num', fui)
  end if

  do s = 1, Grid % n_faces
    if(Front % intersects_face(s)) then

      ! Write down Stefan's solution
      if(Iter % Current() .eq. 1               .and.  &
         Math % Approx_Real(Front % ys(s), 0.0) .and.  &
         Math % Approx_Real(Front % zs(s), 0.0)) then
        write(fui,  '(99(es12.4))')  Time % Curr_Dt() * Flow % dt, Front % xs(s)
      end if

    end if

  end do

  first_entry = .false.

  close(fui)

  end subroutine
