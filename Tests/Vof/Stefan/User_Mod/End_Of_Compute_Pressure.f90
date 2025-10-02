!==============================================================================!
  subroutine User_Mod_End_Of_Compute_Pressure(Flow, Vof, Sol)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Pressure function.           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Flow
  type(Vof_Type),    target :: Vof
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Var_Type),    pointer :: p, pp
  real                       :: x_ref, p_ref, pp_ref
  integer                    :: c, d
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  p    => Flow % p
  pp   => Flow % pp

  do c = Cells_In_Domain_And_Buffers()

    ! Find relevant cell
    if( Math % Approx_Real(Grid % yc(c), 0.0) .and.  &
        Math % Approx_Real(Grid % zc(c), 0.0) ) then
      x_ref = Grid % xc(c)
      p_ref  = p  % n(c)
      pp_ref = pp % n(c)

      ! Browse through all other cells and homogenize the values
      do d = 1, Grid % n_cells
        if(d .ne. c) then
          if(Math % Approx_Real(Grid % xc(d), x_ref)) then
            p  % n(d) = p_ref
            pp % n(d) = pp_ref
          end if
        end if
      end do
    end if
  end do

  end subroutine
