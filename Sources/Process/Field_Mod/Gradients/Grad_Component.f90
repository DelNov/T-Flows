!==============================================================================!
  subroutine Grad_Component(Flow, phi, i, phii)
!------------------------------------------------------------------------------!
!   Calculates gradient of generic variable phi by the least squares method,   !
!   with refershing the buffers.                                               !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow
  real                      :: phi ( -Flow % pnt_grid % n_bnd_cells  &
                                     :Flow % pnt_grid % n_cells)
  integer, intent(in)       :: i
  real,    intent(out)      :: phii( -Flow % pnt_grid % n_bnd_cells  &
                                     :Flow % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: s, c1, c2, reg
  real                     :: dphi1, dphi2
!-----------------------------[Local parameters]-------------------------------!
  integer, dimension(3,3), parameter :: MAP = reshape((/ 1, 4, 5,  &
                                                         4, 2, 6,  &
                                                         5, 6, 3 /), shape(MAP))
!==============================================================================!

  ! Take alias
  Grid => Flow % pnt_grid

  ! Refresh buffers
  call Grid % Exchange_Cells_Real(phi)

  ! Initialize gradients
  phii(1:Grid % n_cells) = 0.

  ! On the boundaries update only c1 - if not symmetry
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .ne. SYMMETRY) then
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        dphi1 = phi(c2)-phi(c1)

        phii(c1) = phii(c1)                                                 &
                 + dphi1 * (  Flow % grad_c2c(MAP(i,1),c1) * Grid % dx(s)   &
                            + Flow % grad_c2c(MAP(i,2),c1) * Grid % dy(s)   &
                            + Flow % grad_c2c(MAP(i,3),c1) * Grid % dz(s))
      end do
    end if  ! symmetry
  end do

  ! Inside the domain
  do s = Faces_In_Domain_And_At_Buffers()
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    dphi1 = phi(c2)-phi(c1)
    dphi2 = phi(c2)-phi(c1)

    phii(c1) = phii(c1)                                                 &
             + dphi1 * (  Flow % grad_c2c(MAP(i,1),c1) * Grid % dx(s)   &
                        + Flow % grad_c2c(MAP(i,2),c1) * Grid % dy(s)   &
                        + Flow % grad_c2c(MAP(i,3),c1) * Grid % dz(s))

    phii(c2) = phii(c2)                                                 &
             + dphi2 * (  Flow % grad_c2c(MAP(i,1),c2) * Grid % dx(s)   &
                        + Flow % grad_c2c(MAP(i,2),c2) * Grid % dy(s)   &
                        + Flow % grad_c2c(MAP(i,3),c2) * Grid % dz(s))
  end do

  call Grid % Exchange_Cells_Real(phii)

  end subroutine
