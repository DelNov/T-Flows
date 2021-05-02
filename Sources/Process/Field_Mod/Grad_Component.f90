!==============================================================================!
  subroutine Field_Mod_Grad_Component(flow, phi, i, phii)
!------------------------------------------------------------------------------!
!   Calculates gradient of generic variable phi by the least squares method,   !
!   with refershing the buffers.                                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  real                     :: phi ( -flow % pnt_grid % n_bnd_cells  &
                                    :flow % pnt_grid % n_cells)
  integer, intent(in)      :: i
  real,    intent(out)     :: phii( -flow % pnt_grid % n_bnd_cells  &
                                    :flow % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: s, c1, c2
  real                     :: dphi1, dphi2
!-----------------------------[Local parameters]-------------------------------!
  integer, dimension(3,3), parameter :: MAP = reshape((/ 1, 4, 5,  &
                                                         4, 2, 6,  &
                                                         5, 6, 3 /), shape(MAP))
!==============================================================================!

  ! Take alias
  grid => flow % pnt_grid

  ! Refresh buffers
  call Grid_Mod_Exchange_Cells_Real(grid, phi)

  ! Initialize gradients
  phii(1:grid % n_cells) = 0.

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    dphi1 = phi(c2) - phi(c1)
    dphi2 = phi(c2) - phi(c1)

    ! On the boundaries
    if(c2 < 0) then

      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. SYMMETRY) then
        dphi1 = 0.
      end if

      phii(c1) = phii(c1)                                                 &
               + dphi1 * (  flow % grad_c2c(MAP(i,1),c1) * grid % dx(s)   &
                          + flow % grad_c2c(MAP(i,2),c1) * grid % dy(s)   &
                          + flow % grad_c2c(MAP(i,3),c1) * grid % dz(s))
    end if

    ! Inside the domain
    if(c2 > 0) then

      phii(c1) = phii(c1)                                                 &
               + dphi1 * (  flow % grad_c2c(MAP(i,1),c1) * grid % dx(s)   &
                          + flow % grad_c2c(MAP(i,2),c1) * grid % dy(s)   &
                          + flow % grad_c2c(MAP(i,3),c1) * grid % dz(s))

      phii(c2) = phii(c2)                                                 &
               + dphi2 * (  flow % grad_c2c(MAP(i,1),c2) * grid % dx(s)   &
                          + flow % grad_c2c(MAP(i,2),c2) * grid % dy(s)   &
                          + flow % grad_c2c(MAP(i,3),c2) * grid % dz(s))
    end if
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, phii)

  end subroutine
