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

  call Comm_Mod_Exchange_Real(grid, p % n)

  do c = 1, grid % n_cells
    p % x(c) = 0.0
    p % y(c) = 0.0
    p % z(c) = 0.0
  end do

  !-------------------------------------------!
  !   First extrapolation to boundary cells   !
  !-------------------------------------------!

  ! User specified gravity vector
  if(sqrt(grav_x ** 2 + grav_y ** 2 + grav_z ** 2) > TINY) then
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 < 0) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. PRESSURE) then
          p % n(c2) = p % n(c1) + density(c1) * (grav_x * grid % dx(s)    &
                                               + grav_y * grid % dy(s)    &
                                               + grav_z * grid % dz(s))
        end if
      end if
    end do

  ! No gravity vector was specified by the user
  else
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 < 0) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. PRESSURE) then
          p % n(c2) = p % n(c1)
        end if
      end if
    end do
  end if

  call Field_Mod_Grad_Component(flow, p % n, 1, p % x)  ! dp/dx
  call Field_Mod_Grad_Component(flow, p % n, 2, p % y)  ! dp/dy
  call Field_Mod_Grad_Component(flow, p % n, 3, p % z)  ! dp/dz

  end subroutine
