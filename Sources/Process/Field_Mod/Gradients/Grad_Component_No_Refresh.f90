!==============================================================================!
  subroutine Grad_Component_No_Refresh(Flow, phi, i, phii)
!------------------------------------------------------------------------------!
!   Calculates gradient of generic variable phi by a least squares method,     !
!   without refershing the buffers.                                            !
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
  integer                  :: s, c1, c2
  real                     :: dphi1, dphi2
!-----------------------------[Local parameters]-------------------------------!
  integer, dimension(3,3), parameter :: MAP = reshape((/ 1, 4, 5,  &
                                                         4, 2, 6,  &
                                                         5, 6, 3 /), shape(MAP))
!==============================================================================!

  ! Take alias
  Grid => Flow % pnt_grid

  ! Initialize gradients
  phii(1:Grid % n_cells) = 0.

  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    dphi1 = phi(c2)-phi(c1)
    dphi2 = phi(c2)-phi(c1)

    ! On the boundaries
    if(c2 < 0) then
      if(Grid % Bnd_Cond_Type(c2) .eq. SYMMETRY) then
        dphi1 = 0.
      end if

      phii(c1) = phii(c1)                                                 &
               + dphi1 * (  Flow % grad_c2c(MAP(i,1),c1) * Grid % dx(s)   &
                          + Flow % grad_c2c(MAP(i,2),c1) * Grid % dy(s)   &
                          + Flow % grad_c2c(MAP(i,3),c1) * Grid % dz(s))
    end if

    ! Inside the domain
    if(c2 > 0) then
      phii(c1) = phii(c1)                                                 &
               + dphi1 * (  Flow % grad_c2c(MAP(i,1),c1) * Grid % dx(s)   &
                          + Flow % grad_c2c(MAP(i,2),c1) * Grid % dy(s)   &
                          + Flow % grad_c2c(MAP(i,3),c1) * Grid % dz(s))

      phii(c2) = phii(c2)                                                 &
               + dphi2 * (  Flow % grad_c2c(MAP(i,1),c2) * Grid % dx(s)   &
                          + Flow % grad_c2c(MAP(i,2),c2) * Grid % dy(s)   &
                          + Flow % grad_c2c(MAP(i,3),c2) * Grid % dz(s))
    end if
  end do

  end subroutine
