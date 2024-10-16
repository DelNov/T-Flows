!==============================================================================!
  subroutine Grad_Component_Faces_To_Cells(Flow, phic, phif, i, phii)
!------------------------------------------------------------------------------!
!   Calculates gradient of a cell-based cell array phic by the least squares   !
!   method from facial values in face-based array phif.  This function is very !
!   similar to more usual cell-based version for calculation of gradients,     !
!   and was introduced as an addition, or alternative, to Gauss based method   !
!   of gradient computation.  It is considered experimental, and I am not sure !
!   if it will ever be used but could prove useful when implementing force-    !
!   based face interpolations if Gaussian method proves weak or unstable.      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow
  real                      :: phic( -Flow % pnt_grid % n_bnd_cells  &
                                     :Flow % pnt_grid % n_cells)
  real                      :: phif(1:Flow % pnt_grid % n_faces)
  integer, intent(in)       :: i
  real,    intent(out)      :: phii( -Flow % pnt_grid % n_bnd_cells  &
                                     :Flow % pnt_grid % n_cells)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: s, c1, c2
  real                     :: dphi1, dphi2
  real                     :: dx_c1, dy_c1, dz_c1, dx_c2, dy_c2, dz_c2
!-----------------------------[Local parameters]-------------------------------!
  integer, dimension(3,3), parameter :: MAP = reshape((/ 1, 4, 5,  &
                                                         4, 2, 6,  &
                                                         5, 6, 3 /), shape(MAP))
!==============================================================================!

  ! Take alias
  Grid => Flow % pnt_grid

  ! Refresh buffers
  call Grid_Mod_Exchange_Cells_Real(Grid, phic)

  ! Initialize gradients
  phii(1:Grid % n_cells) = 0.

  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    dphi1 = phif(s) - phic(c1)
    dphi2 = phif(s) - phic(c2)

    dx_c1 = Grid % xf(s) - Grid % xc(c1)
    dy_c1 = Grid % yf(s) - Grid % yc(c1)
    dz_c1 = Grid % zf(s) - Grid % zc(c1)
    dx_c2 = Grid % xf(s) - Grid % xc(c2)
    dy_c2 = Grid % yf(s) - Grid % yc(c2)
    dz_c2 = Grid % zf(s) - Grid % zc(c2)

    ! On the boundaries
    if(c2 < 0) then

      phii(c1) = phii(c1)                                          &
               + dphi1 * (  Flow % grad_f2c(MAP(i,1),c1) * dx_c1   &
                          + Flow % grad_f2c(MAP(i,2),c1) * dy_c1   &
                          + Flow % grad_f2c(MAP(i,3),c1) * dz_c1)
    end if

    ! Inside the domain
    if(c2 > 0) then

      phii(c1) = phii(c1)                                          &
               + dphi1 * (  Flow % grad_f2c(MAP(i,1),c1) * dx_c1   &
                          + Flow % grad_f2c(MAP(i,2),c1) * dy_c1   &
                          + Flow % grad_f2c(MAP(i,3),c1) * dz_c1)

      phii(c2) = phii(c2)                                          &
               + dphi2 * (  Flow % grad_f2c(MAP(i,1),c2) * dx_c2   &
                          + Flow % grad_f2c(MAP(i,2),c2) * dy_c2   &
                          + Flow % grad_f2c(MAP(i,3),c2) * dz_c2)
    end if
  end do

  call Grid_Mod_Exchange_Cells_Real(Grid, phii)

  end subroutine
