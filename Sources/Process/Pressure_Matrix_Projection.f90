!==============================================================================!
  subroutine Pressure_Matrix_Projection(grid, dt)
!------------------------------------------------------------------------------!
!   Forms the pressure system matrix for the fractional step method.           !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Flow_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: dt
!-----------------------------------[Locals]-----------------------------------!
  real    :: a12
  integer :: c, c1, c2, s 
!==============================================================================!

  do c = 1, a % row(grid % n_cells+1)  ! this is number of nozero entries + 1
    a % val(c) = 0.0
  end do

  !-----------------------------!
  !   Calculate system matrix   ! 
  !-----------------------------!
  do s = 1, grid % n_faces

    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 > 0) then
      a12 = dt * f_coef(s)
      a % val(a % pos(1,s)) = -a12
      a % val(a % pos(2,s)) = -a12
      a % val(a % dia(c1))  =              &
      a % val(a % dia(c1))  +  a12
      a % val(a % dia(c2))  =              &
      a % val(a % dia(c2))  +  a12
    end if

  end do ! through faces

  end subroutine
