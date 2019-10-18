!==============================================================================!
  subroutine Grad_Mod_Pressure(phi)
!------------------------------------------------------------------------------!
!   Calculates gradient of pressure of pressure correction.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Var_Type), target :: phi
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: s, c, c1, c2, iter
!==============================================================================!

  ! Take aliases
  grid  => phi % pnt_grid

  call Comm_Mod_Exchange_Real(grid, phi % n)

  do c = 1, grid % n_cells
    phi % x(c) = 0.0
    phi % y(c) = 0.0
    phi % z(c) = 0.0
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
          phi % n(c2) = phi % n(c1) + density(c1) * (grav_x * grid % dx(s)    &
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
          phi % n(c2) = phi % n(c1)
        end if
      end if
    end do
  end if

  call Grad_Mod_Component(grid, phi % n, 1, phi % x)  ! dp/dx
  call Grad_Mod_Component(grid, phi % n, 2, phi % y)  ! dp/dy
  call Grad_Mod_Component(grid, phi % n, 3, phi % z)  ! dp/dz

  end subroutine
