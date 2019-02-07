!==============================================================================!
  subroutine Grad_Mod_Pressure(phi)
!------------------------------------------------------------------------------!
!   Calculates gradient of pressure of pressure correction.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Var_Mod,  only: Var_Type
  use Grid_Mod, only: Grid_Type, Grid_Mod_Bnd_Cond_Type, PRESSURE
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

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. PRESSURE) then  
        phi % n(c2) = phi % n(c1)
      end if
    end if
  end do

  call Grad_Mod_Component(grid, phi % n, 1, phi % x, .true.)  ! dp/dx
  call Grad_Mod_Component(grid, phi % n, 2, phi % y, .true.)  ! dp/dy
  call Grad_Mod_Component(grid, phi % n, 3, phi % z, .true.)  ! dp/dz

  do iter=1, 1

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 < 0) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. PRESSURE) then
          phi % n(c2) = phi % n(c1)   &
                      + 1.2*( phi % x(c1) * (grid % xc(c2)-grid % xc(c1))  &
                      +       phi % y(c1) * (grid % yc(c2)-grid % yc(c1))  &
                      +       phi % z(c1) * (grid % zc(c2)-grid % zc(c1))  &
                            )
        end if
      end if
    end do

    call Grad_Mod_Component(grid, phi % n, 1, phi % x, .true.)  ! dp/dx
    call Grad_Mod_Component(grid, phi % n, 2, phi % y, .true.)  ! dp/dy
    call Grad_Mod_Component(grid, phi % n, 3, phi % z, .true.)  ! dp/dz

  end do

  end subroutine
