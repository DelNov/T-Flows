!==============================================================================!
  subroutine User_Mod_End_Of_Compute_Vof(Vof, Sol)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Vof function.                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Vof_Type),    target :: Vof
  type(Solver_Type), target :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: Grid
  type(Field_Type),  pointer :: Flow
  type(Var_Type),    pointer :: fun
  type(Matrix_Type), pointer :: A
  integer                    :: c, d
  real                       :: x_ref, f_ref
!==============================================================================!

  ! Take aliases
  Flow => Vof  % pnt_flow
  Grid => Flow % pnt_grid
  fun  => Vof % fun
  A    => Sol % Nat % A

  do c = Cells_In_Domain_And_Buffers()

    ! Find relevant cell
    if( Math % Approx_Real(Grid % yc(c), 0.0) .and.  &
        Math % Approx_Real(Grid % zc(c), 0.0) ) then
      x_ref = Grid % xc(c)
      f_ref = fun % n(c)

      ! Browse through all other cells and homogenize the values
      do d = 1, Grid % n_cells
        if(d .ne. c) then
          if(Math % Approx_Real(Grid % xc(d), x_ref)) then
            fun % n(d) = f_ref
          end if
        end if
      end do
    end if
  end do

  end subroutine
