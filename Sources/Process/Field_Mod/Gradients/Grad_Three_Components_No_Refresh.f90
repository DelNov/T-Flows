!==============================================================================!
  subroutine Grad_Three_Components_No_Refresh(Flow, phi, phi_x, phi_y, phi_z)
!------------------------------------------------------------------------------!
!   Calculates gradient of generic variable phi by a least squares method,     !
!   without refershing the buffers.                                            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow
  real,    intent(in)       :: phi ( -Flow % pnt_grid % n_bnd_cells  &
                                     :Flow % pnt_grid % n_cells)
  real,    intent(out)      :: phi_x( -Flow % pnt_grid % n_bnd_cells  &
                                      :Flow % pnt_grid % n_cells)
  real,    intent(out)      :: phi_y( -Flow % pnt_grid % n_bnd_cells  &
                                      :Flow % pnt_grid % n_cells)
  real,    intent(out)      :: phi_z( -Flow % pnt_grid % n_bnd_cells  &
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
  phi_x(1:Grid % n_cells) = 0.
  phi_y(1:Grid % n_cells) = 0.
  phi_z(1:Grid % n_cells) = 0.

  ! On the boundaries update only c1 - if not symmetry
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .ne. SYMMETRY) then
      !$omp parallel do                                                      &
      !$omp private(s, c1, c2, dphi1)                                        &
      !$omp shared(faces_c, grad_c2c, dx, dy, dz, phi, phi_x, phi_y, phi_z)
      do s = Faces_In_Region(reg)
        c1 = faces_c(1,s)
        c2 = faces_c(2,s)

        dphi1 = phi(c2)-phi(c1)

        phi_x(c1) = phi_x(c1)                                  &
                  + dphi1 * (  grad_c2c(MAP(1,1),c1) * dx(s)   &
                             + grad_c2c(MAP(1,2),c1) * dy(s)   &
                             + grad_c2c(MAP(1,3),c1) * dz(s))
        phi_y(c1) = phi_y(c1)                                  &
                  + dphi1 * (  grad_c2c(MAP(2,1),c1) * dx(s)   &
                             + grad_c2c(MAP(2,2),c1) * dy(s)   &
                             + grad_c2c(MAP(2,3),c1) * dz(s))
        phi_z(c1) = phi_z(c1)                                  &
                  + dphi1 * (  grad_c2c(MAP(3,1),c1) * dx(s)   &
                             + grad_c2c(MAP(3,2),c1) * dy(s)   &
                             + grad_c2c(MAP(3,3),c1) * dz(s))
      end do
      !$omp end parallel do
    end if  ! symmetry
  end do

  ! Inside the domain
  !$omp parallel do                                                      &
  !$omp private(s, c1, c2, dphi1, dphi2)                                 &
  !$omp shared(faces_c, grad_c2c, dx, dy, dz, phi, phi_x, phi_y, phi_z)
  do s = Faces_In_Domain()
    c1 = faces_c(1,s)
    c2 = faces_c(2,s)

    dphi1 = phi(c2)-phi(c1)
    dphi2 = phi(c2)-phi(c1)

    phi_x(c1) = phi_x(c1)                                  &
              + dphi1 * (  grad_c2c(MAP(1,1),c1) * dx(s)   &
                         + grad_c2c(MAP(1,2),c1) * dy(s)   &
                         + grad_c2c(MAP(1,3),c1) * dz(s))
    phi_y(c1) = phi_y(c1)                                  &
              + dphi1 * (  grad_c2c(MAP(2,1),c1) * dx(s)   &
                         + grad_c2c(MAP(2,2),c1) * dy(s)   &
                         + grad_c2c(MAP(2,3),c1) * dz(s))
    phi_z(c1) = phi_z(c1)                                  &
              + dphi1 * (  grad_c2c(MAP(3,1),c1) * dx(s)   &
                         + grad_c2c(MAP(3,2),c1) * dy(s)   &
                         + grad_c2c(MAP(3,3),c1) * dz(s))

    phi_x(c2) = phi_x(c2)                                  &
              + dphi2 * (  grad_c2c(MAP(1,1),c2) * dx(s)   &
                         + grad_c2c(MAP(1,2),c2) * dy(s)   &
                         + grad_c2c(MAP(1,3),c2) * dz(s))
    phi_y(c2) = phi_y(c2)                                  &
              + dphi2 * (  grad_c2c(MAP(2,1),c2) * dx(s)   &
                         + grad_c2c(MAP(2,2),c2) * dy(s)   &
                         + grad_c2c(MAP(2,3),c2) * dz(s))
    phi_z(c2) = phi_z(c2)                                  &
              + dphi2 * (  grad_c2c(MAP(3,1),c2) * dx(s)   &
                         + grad_c2c(MAP(3,2),c2) * dy(s)   &
                         + grad_c2c(MAP(3,3),c2) * dz(s))
  end do
  !$omp end parallel do

  end subroutine
