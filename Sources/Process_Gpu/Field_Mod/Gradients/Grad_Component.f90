!==============================================================================!
  subroutine Grad_Component(Flow, Grid, phi, i, phii, boundary_updated)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), intent(in), target  :: Flow  !! parent flow object
  type(Grid_Type),   intent(in), target  :: Grid  !! grid object
  real                                   :: phi (-Grid % n_bnd_cells &
                                                 :Grid % n_cells)
  integer,           intent(in)          :: i     !! gradient component (1 to 3)
  real,              intent(out)         :: phii(-Grid % n_bnd_cells &
                                                 :Grid % n_cells)
  logical, optional, intent(in)          :: boundary_updated
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s, i_cel, reg
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

  ! Copy values to symmetry boundaries
  if(present(boundary_updated)) then
    if(.not. boundary_updated) then
      do reg = Boundary_Regions()
        if(Grid % region % type(reg) .eq. SYMMETRY) then

          !$acc parallel loop  &
          !$acc present(  &
          !$acc   grid_region_f_face,  &
          !$acc   grid_region_l_face,  &
          !$acc   grid_faces_c,  &
          !$acc   phi   &
          !$acc )
          do s = grid_region_f_face(reg), grid_region_l_face(reg)  ! all present
            c1 = grid_faces_c(1,s)   ! inside cell
            c2 = grid_faces_c(2,s)   ! boundary cell
            phi(c2) = phi(c1)
          end do
          !$acc end parallel

        end if
      end do
    end if
  end if

  ! Estimate gradients cell-wise (face-wise leades to race conditions on GPUs)
  !$acc parallel loop independent  &
  !$acc present(  &
  !$acc   grid_region_f_cell,  &
  !$acc   grid_region_l_cell,  &
  !$acc   phi,  &
  !$acc   phii,  &
  !$acc   grid_cells_n_cells,  &
  !$acc   grid_cells_c,  &
  !$acc   grid_cells_f,  &
  !$acc   flow_grad_c2c,  &
  !$acc   map,  &
  !$acc   grid_dx,  &
  !$acc   grid_dy,  &
  !$acc   grid_dz   &
  !$acc )
  do c1 = grid_region_f_cell(grid_n_regions), grid_region_l_cell(grid_n_regions)  ! all present

    phii_tmp = 0.0

  !$acc loop seq
    do i_cel = 1, grid_cells_n_cells(c1)
      c2 = grid_cells_c(i_cel, c1)
      s  = grid_cells_f(i_cel, c1)
      dphi = phi(c2) - phi(c1)
      if(c2 .gt. 0 .and. c1 .gt. c2) then
        dphi = -dphi
      end if

      phii_tmp = phii_tmp                                                &
               + dphi * (  flow_grad_c2c(map(i,1),c1) * grid_dx(s)   &
                         + flow_grad_c2c(map(i,2),c1) * grid_dy(s)   &
                         + flow_grad_c2c(map(i,3),c1) * grid_dz(s))
    end do
  !$acc end loop

    phii(c1) = phii_tmp

  end do
  !$acc end parallel

  call Grid % Exchange_Cells_Real(phii)

  end subroutine
