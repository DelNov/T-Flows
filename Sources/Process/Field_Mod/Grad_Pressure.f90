!==============================================================================!
  subroutine Field_Mod_Grad_Pressure(flow,     &
                                     p,        &
                                     density,  &
                                     grav_x, grav_y, grav_z)
!------------------------------------------------------------------------------!
!   Calculates gradient of pressure of pressure correction.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type) :: flow
  type(Var_Type)   :: p
  real             :: density(-p % pnt_grid % n_bnd_cells:p % pnt_grid % n_cells)
  real             :: grav_x, grav_y, grav_z
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: s, c, c1, c2, iter
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid

  call Grid_Mod_Exchange_Real(grid, p % n)

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. PRESSURE) then
        p % n(c2) = p % n(c1) + dot_product( (/p % x(c1),        &
                                               p % y(c1),        &
                                               p % z(c1)/) ,     &
                                             (/grid % dx(s),     &
                                               grid % dy(s),     &
                                               grid % dz(s)/) )
      end if
    end if
  end do

  call Field_Mod_Grad_Component(flow, p % n, 1, p % x)  ! dp/dx
  call Field_Mod_Grad_Component(flow, p % n, 2, p % y)  ! dp/dy
  call Field_Mod_Grad_Component(flow, p % n, 3, p % z)  ! dp/dz

  end subroutine
