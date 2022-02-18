!==============================================================================!
  subroutine User_Mod_End_Of_Compute_Energy(Flow, turb, Vof, Nat, curr_dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Energy function.             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: turb
  type(Vof_Type),      target :: Vof
  type(Native_Type),   target :: Nat
  integer, intent(in)         :: curr_dt  ! current time step
  integer, intent(in)         :: ini      ! inner iteration
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: s, fu
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid

  call File % Append_For_Writing_Ascii('stefans_solution.dat', fu)

  do s = 1, Grid % n_faces

    if(any(Vof % Front % face_at_elem(1:2,s) .ne. 0)) then

      ! Write down Stefan's solution
      if(ini .eq. 1                            .and.  &
         Math % Approx_Real(grid % ys(s), 0.0) .and.  &
         Math % Approx_Real(grid % zs(s), 0.0)) then
        write(fu,  '(99(es12.4))') curr_dt * Flow % dt, Grid % xs(s)
      end if

    end if

  end do

  close(fu)

  end subroutine
