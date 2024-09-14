!==============================================================================!
  subroutine Grad_Component(Flow, Grid, phi, i, phii)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), intent(in), target  :: Flow  !! parent flow object
  type(Grid_Type),   intent(in), target  :: Grid  !! grid object
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
  !$tf-acc loop begin
  do c1 = Cells_In_Domain()  ! all present

    phii_tmp = 0.0

    do i_cel = 1, Grid % cells_n_cells(c1)
      c2 = Grid % cells_c(i_cel, c1)
      s  = Grid % cells_f(i_cel, c1)
      dphi = phi(c2) - phi(c1)
      if(c2 .gt. 0 .and. c1 .gt. c2) then
        dphi = -dphi
      end if

      phii_tmp = phii_tmp                                                &
               + dphi * (  Flow % grad_c2c(map(i,1),c1) * Grid % dx(s)   &
                         + Flow % grad_c2c(map(i,2),c1) * Grid % dy(s)   &
                         + Flow % grad_c2c(map(i,3),c1) * Grid % dz(s))
    end do

    phii(c1) = phii_tmp

  end do
  !$tf-acc loop end

  call Grid % Exchange_Cells_Real(phii)

  end subroutine
