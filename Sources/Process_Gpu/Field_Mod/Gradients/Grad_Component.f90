!==============================================================================!
  subroutine Grad_Component(Flow, Grid, phi, i, phii)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), intent(in)  :: Flow  !! parent flow object
  type(Grid_Type),   intent(in)  :: Grid  !! grid object
  real,              intent(in)  :: phi (-Grid % n_bnd_cells:Grid % n_cells)
  integer,           intent(in)  :: i     !! gradient component (1 to 3)
  real,              intent(out) :: phii(-Grid % n_bnd_cells:Grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s, i_cel
  real    :: dphi, dx, dy, dz, phii_tmp
!-----------------------------[Local parameters]-------------------------------!
  integer, dimension(3,3), parameter :: MAP = reshape((/ 1, 4, 5,  &
                                                         4, 2, 6,  &
                                                         5, 6, 3 /), shape(MAP))
!==============================================================================!

  ! Aret these checks overkill?
  Assert(i .ge. 1)
  Assert(i .le. 3)

  ! Initialize gradients
  phii(:) = 0.0

  ! Try to estimate gradients cell-wise
  ! (face-wise leades to race conditions on GPUs)
  !$acc parallel loop independent
  do c1 = 1, Grid % n_cells - Grid % Comm % n_buff_cells

    phii_tmp = 0.0

    !$acc loop seq
    do i_cel = 1, Grid % cells_n_cells(c1)
      c2 = Grid % cells_c(i_cel, c1)
      s  = Grid % cells_f(i_cel, c1)
      dphi = phi(c2)-phi(c1)
      dx = Grid % dx(s)
      dy = Grid % dy(s)
      dz = Grid % dz(s)
      if(c2 .gt. 0 .and. c1 .gt. c2) then
        dx = -dx
        dy = -dy
        dz = -dz
      end if

      phii_tmp = phii_tmp                                      &
               + dphi * (  Flow % grad_c2c(MAP(i,1),c1) * dx   &
                         + Flow % grad_c2c(MAP(i,2),c1) * dy   &
                         + Flow % grad_c2c(MAP(i,3),c1) * dz)
    end do
    !$acc end loop

    phii(c1) = phii_tmp

  end do
  !$acc end parallel

  call Grid % Exchange_Cells_Real(phii)

  end subroutine
