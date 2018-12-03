!==============================================================================!
  subroutine Grad_Mod_Pressure(grid, phi, phi_x, phi_y, phi_z)
!------------------------------------------------------------------------------!
!   Calculates gradient of pressure of pressure correction.
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: phi  (-grid % n_bnd_cells:grid % n_cells),  &
                     phi_x(-grid % n_bnd_cells:grid % n_cells),  &
                     phi_y(-grid % n_bnd_cells:grid % n_cells),  &
                     phi_z(-grid % n_bnd_cells:grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer :: s, c, c1, c2, iter
!==============================================================================!
 
  call Comm_Mod_Exchange_Real(grid, phi)

  do c = 1, grid % n_cells
    phi_x(c)=0.0
    phi_y(c)=0.0
    phi_z(c)=0.0
  end do

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. PRESSURE) then  
        phi(c2) = phi(c1)
      end if
    end if
  end do

  call Grad_Mod_Component(grid, phi, 1, phi_x, .true.)  ! dp/dx
  call Grad_Mod_Component(grid, phi, 2, phi_y, .true.)  ! dp/dy
  call Grad_Mod_Component(grid, phi, 3, phi_z, .true.)  ! dp/dz

  do iter=1, 1

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 < 0) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. PRESSURE) then
          phi(c2) = phi(c1) + 1.2*( phi_x(c1) * (grid % xc(c2)-grid % xc(c1))  &
                            +       phi_y(c1) * (grid % yc(c2)-grid % yc(c1))  &
                            +       phi_z(c1) * (grid % zc(c2)-grid % zc(c1))  &
                            )
        end if
      end if
    end do

    call Grad_Mod_Component(grid, phi, 1, phi_x, .true.)  ! dp/dx
    call Grad_Mod_Component(grid, phi, 2, phi_y, .true.)  ! dp/dy
    call Grad_Mod_Component(grid, phi, 3, phi_z, .true.)  ! dp/dz 

  end do

  end subroutine
