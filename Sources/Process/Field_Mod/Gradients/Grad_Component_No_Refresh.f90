!==============================================================================!
  subroutine Grad_Component_No_Refresh(Flow, phi, i, phii)
!------------------------------------------------------------------------------!
!   Calculates gradient of generic variable phi by a least squares method,     !
!   without refershing the buffers.                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow
  real,    intent(in)       :: phi ( -Flow % pnt_grid % n_bnd_cells  &
                                     :Flow % pnt_grid % n_cells)
  integer, intent(in)       :: i
  real,    intent(out)      :: phii( -Flow % pnt_grid % n_bnd_cells  &
                                     :Flow % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),     pointer :: Grid
  real,    contiguous, pointer :: dx(:), dy(:), dz(:), grad_c2c(:,:)
  integer, contiguous, pointer :: faces_c(:,:)
  integer                      :: s, c1, c2, reg
  real                         :: dphi1, dphi2
!-----------------------------[Local parameters]-------------------------------!
  integer, dimension(3,3), parameter :: MAP = reshape((/ 1, 4, 5,  &
                                                         4, 2, 6,  &
                                                         5, 6, 3 /), shape(MAP))
!==============================================================================!

  ! Take alias
  ! OpenMP doesn't unerstand Fortran's members (%), that's why ...
  ! ... aliases for faces_c, grad_c2c, dx, dy and dz are needed
  Grid     => Flow % pnt_grid
  dx       => Grid % dx
  dy       => Grid % dy
  dz       => Grid % dz
  grad_c2c => Flow % grad_c2c
  faces_c  => Grid % faces_c

  ! Initialize gradients
  phii(1:Grid % n_cells) = 0.

  ! On the boundaries update only c1 - if not symmetry
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .ne. SYMMETRY) then
      !$omp parallel do                                       &
      !$omp private(s, c1, c2, dphi1)                         &
      !$omp shared(faces_c, grad_c2c, dx, dy, dz, phi, phii)
      do s = Faces_In_Region(reg)
        c1 = faces_c(1,s)
        c2 = faces_c(2,s)

        dphi1 = phi(c2)-phi(c1)

        phii(c1) = phii(c1)                                   &
                 + dphi1 * (  grad_c2c(MAP(i,1),c1) * dx(s)   &
                            + grad_c2c(MAP(i,2),c1) * dy(s)   &
                            + grad_c2c(MAP(i,3),c1) * dz(s))
      end do
      !$omp end parallel do
    end if  ! symmetry
  end do

  ! Inside the domain
  !$omp parallel do                                       &
  !$omp private(s, c1, c2, dphi1, dphi2)                  &
  !$omp shared(faces_c, grad_c2c, dx, dy, dz, phi, phii)
  do s = Faces_In_Domain_And_At_Buffers()
    c1 = faces_c(1,s)
    c2 = faces_c(2,s)

    dphi1 = phi(c2)-phi(c1)
    dphi2 = phi(c2)-phi(c1)

    phii(c1) = phii(c1)                                   &
             + dphi1 * (  grad_c2c(MAP(i,1),c1) * dx(s)   &
                        + grad_c2c(MAP(i,2),c1) * dy(s)   &
                        + grad_c2c(MAP(i,3),c1) * dz(s))

    phii(c2) = phii(c2)                                   &
             + dphi2 * (  grad_c2c(MAP(i,1),c2) * dx(s)   &
                        + grad_c2c(MAP(i,2),c2) * dy(s)   &
                        + grad_c2c(MAP(i,3),c2) * dz(s))
  end do
  !$omp end parallel do

  end subroutine
