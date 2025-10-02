!==============================================================================!
  subroutine Grad_Component(Flow, Grid, phi, i, phii)
!------------------------------------------------------------------------------!
!>  Calculates one gradient component of generic variable phi by the least
!>  squares method, with refreshing the buffers of the variable values before
!>  the calculation and of the calculated gradient components.  (There is a
!>  procedure which does the same, but it doesn't refereshes the buffers.)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), intent(in)    :: Flow  !! parent flow object
  type(Grid_Type),   intent(in)    :: Grid  !! grid object
  real,              intent(inout) :: phi (-Grid % n_bnd_cells:Grid % n_cells)
    !! field whose gradients are being calculated
  integer,           intent(in)    :: i     !! gradient component (1 to 3)
  real,              intent(out)   :: phii(-Grid % n_bnd_cells:Grid % n_cells)
    !! calculated gradient in direction specified by i
!-----------------------------------[Locals]-----------------------------------!
  integer                  :: s, c1, c2, reg
  real                     :: dphi1, dphi2
!-----------------------------[Local parameters]-------------------------------!
  integer, dimension(3,3), parameter :: MAP = reshape((/ 1, 4, 5,  &
                                                         4, 2, 6,  &
                                                         5, 6, 3 /), shape(MAP))
!==============================================================================!

  ! Aret these checks overkill?
  Assert(i > 0)
  Assert(i < 4)

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
