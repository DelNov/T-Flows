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

  ! Estimate gradients cell-wise (face-wise leades to race conditions on GPUs)
  !$acc parallel loop independent                                &
  !$acc present(grid_cells_n_cells, grid_cells_c, grid_cells_f,  &
  !$acc         grid_dx, grid_dy, grid_dz,                       &
  !$acc         flow_grad_c2c)
  do c1 = Cells_In_Domain()

    phii_tmp = 0.0

    !$acc loop seq
    do i_cel = 1, grid_cells_n_cells(c1)
      c2 = grid_cells_c(i_cel, c1)
      s  = grid_cells_f(i_cel, c1)
      dphi = phi(c2)-phi(c1)
      dx = grid_dx(s)
      dy = grid_dy(s)
      dz = grid_dz(s)
      if(c2 .gt. 0 .and. c1 .gt. c2) then
        dx = -dx
        dy = -dy
        dz = -dz
      end if

      phii_tmp = phii_tmp                                    &
               + dphi * (  flow_grad_c2c(MAP(i,1),c1) * dx   &
                         + flow_grad_c2c(MAP(i,2),c1) * dy   &
                         + flow_grad_c2c(MAP(i,3),c1) * dz)
    end do
    !$acc end loop

    phii(c1) = phii_tmp

  end do
  !$acc end parallel

  call Grid % Exchange_Cells_Real(phii)

  end subroutine
