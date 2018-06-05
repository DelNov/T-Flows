!==============================================================================!
  subroutine Calculate_Minimum_Maximum(grid, phi, phi_min, phi_max)
!------------------------------------------------------------------------------!
!                                                                              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Comm_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: phi    (-grid % n_bnd_cells:grid % n_cells) 
  real            :: phi_min(-grid % n_bnd_cells:grid % n_cells) 
  real            :: phi_max(-grid % n_bnd_cells:grid % n_cells) 
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s
!==============================================================================!

  phi_min(:) = phi(:)
  phi_max(:) = phi(:)

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1, s)
    c2 = grid % faces_c(2, s)

    if(c2  < 0) then
      if (Grid_Mod_Bnd_Cond_Type(grid, c2) .ne. BUFFER) then
        cycle ! skip current surface
      end if
    end if
      
    phi_min(c1) = min(phi_min(c1), phi(c2))
    phi_min(c2) = min(phi_min(c2), phi(c1))

    phi_max(c1) = max(phi_max(c1), phi(c2))
    phi_max(c2) = max(phi_max(c2), phi(c1))

  end do

  call Comm_Mod_Exchange(grid, phi_min)
  call Comm_Mod_Exchange(grid, phi_max)

  end subroutine
