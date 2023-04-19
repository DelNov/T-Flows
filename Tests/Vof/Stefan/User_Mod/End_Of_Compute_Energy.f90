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
  type(Grid_Type), pointer :: Grid
  integer                  :: s, fu
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid

  call File % Append_For_Writing_Ascii('stefans_solution.dat', fu)

  do s = 1, Grid % n_faces

    if(any(Vof % Front % elems_at_face(1:2,s) .ne. 0)) then

      ! Write down Stefan's solution
      if(Iter % Current() .eq. 1               .and.  &
         Math % Approx_Real(Grid % ys(s), 0.0) .and.  &
         Math % Approx_Real(Grid % zs(s), 0.0)) then
        write(fu,  '(99(es12.4))') Time % Curr_Dt() * Flow % dt, Grid % xs(s)
      end if

    end if

  end do

  close(fu)

  end subroutine
