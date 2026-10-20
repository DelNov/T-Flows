!==============================================================================!
  subroutine Grad_Least_Pressure(Flow, p)
!------------------------------------------------------------------------------!
!>  This subroutine calculates gradient of pressure (or pressure correction)
!>  using the least squares method.
!------------------------------------------------------------------------------!
!   Functionality                                                              !
!                                                                              !
!   * Buffer refresh: The subroutine starts by refreshing the buffers for the  !
!     pressure (or pressure correction) variable p % n.                        !
!   * Iterative calculation: It then enters an iterative loop (up to a         !
!     maximum number of iterations defined by Flow % least_miter):             !
!     - Boundary extrapolation: For each boundary region which is not of the   !
!       PRESSURE type, it extrapolates the pressure values to boundary cells   !
!       based on the neighboring cells' pressure and gradient values.          !
!     - Gradient calculation: Individual gradient components (p % x, p % y,    !
!       and p % z for dp/dx, dp/dy, dp/dz respectively) are calculated for     !
!       each cell in the domain without refreshing buffers.                    !
!     - Buffer refresh for gradient components: After calculating the gradient !
!       components, the subroutine refreshes the buffers for each component.   !
!   * Balanced with body forces (Flow % balanced_p_grad, pressure only):       !
!     the gradient is computed as grad(p) = f_c + G(dp - f_s . d), where f_c   !
!     is the cell-centered body force, G the least squares operator, dp the    !
!     pressure difference over a face and f_s . d the pressure difference      !
!     body forces sustain over that face at rest; f_s is the face force if     !
!     Gu's correction is used and linear interpolation of cell forces if not.  !
!     At boundaries, only the deviation of pressure gradient from the body     !
!     force is extrapolated.  This way, grad(p) equals f_c in hydrostatic      !
!     equilibrium, also in cells next to boundaries, and spurious velocities   !
!     in quiescent, stratified fluid vanish on orthogonal grids and are much   !
!     reduced on non-orthogonal ones (most with Gu's correction).  Without     !
!     body forces, it is the same as the standard method.                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow
  type(Var_Type)            :: p
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: s, c, c1, c2, iter, reg
  logical                  :: balanced
  real                     :: fs, fx_s, fy_s, fz_s
  real, contiguous, pointer :: p_res(:)
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid

  ! Refresh buffers for variable
  call Grid % Exchange_Cells_Real(p % n)

  balanced = Flow % balanced_p_grad .and. p % name .eq. 'P'

  !---------------------------!
  !   The standard approach   !
  !---------------------------!
  if(.not. balanced) then

    do iter = 1, Flow % least_miter

      ! Extrapolation to boundaries
      do reg = Boundary_Regions()
        if(Grid % region % type(reg) .ne. AMBIENT) then
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)

            p % n(c2) = p % n(c1) + p % x(c1) * Grid % dx(s)  &
                                  + p % y(c1) * Grid % dy(s)  &
                                  + p % z(c1) * Grid % dz(s)
          end do  ! faces
        else
          ! This was used for checking
          ! do s = Faces_In_Region(reg)
          !   c2 = Grid % faces_c(2,s)
          ! end do
        end if    ! pressure, ambient or other
      end do      ! regions

      ! Compute individual gradients without refreshing buffers
      call Flow % Grad_Component_No_Refresh(Grid, p % n, 1, p % x)  ! dp/dx
      call Flow % Grad_Component_No_Refresh(Grid, p % n, 2, p % y)  ! dp/dy
      call Flow % Grad_Component_No_Refresh(Grid, p % n, 3, p % z)  ! dp/dz

      ! Refresh buffers for gradient components
      call Grid % Exchange_Cells_Real(p % x)
      call Grid % Exchange_Cells_Real(p % y)
      call Grid % Exchange_Cells_Real(p % z)

    end do

  !--------------------------------------!
  !   Approach balanced with the body    !
  !   forces (see the comment on top)    !
  !--------------------------------------!
  else

    call Work % Connect_Real_Face(p_res)

    do iter = 1, Flow % least_miter

      ! Extrapolation to boundaries: hydrostatic part with the face force,
      ! and only the deviation of the gradient from the body force with the
      ! gradient (the face force at boundary is the only one available)
      do reg = Boundary_Regions()
        if(Grid % region % type(reg) .ne. AMBIENT) then
          do s = Faces_In_Region(reg)
            c1 = Grid % faces_c(1,s)
            c2 = Grid % faces_c(2,s)

            p % n(c2) = p % n(c1)                                         &
                      + (p % x(c1) - Flow % cell_fx(c1)) * Grid % dx(s)   &
                      + (p % y(c1) - Flow % cell_fy(c1)) * Grid % dy(s)   &
                      + (p % z(c1) - Flow % cell_fz(c1)) * Grid % dz(s)   &
                      + Flow % face_fx(s) * Grid % dx(s)                  &
                      + Flow % face_fy(s) * Grid % dy(s)                  &
                      + Flow % face_fz(s) * Grid % dz(s)
          end do  ! faces
        end if
      end do      ! regions

      ! Pressure differences not sustained by the body forces.  Note
      ! that Grad_Component_No_Refresh works with differences over faces
      ! (phi(c2)-phi(c1)), so they are stored here in a face array and
      ! the gradient is computed from them in Grad_Residual below
      do s = 1, Grid % n_faces
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        ! Force sustaining pressure difference at rest: linear interpolation
        ! of cell forces without Gu's correction, face force with it (and at
        ! boundaries, where the face force is the only one available)
        if(c2 > 0 .and. .not. Flow % gu_correction) then
          fs   = Grid % f(s)
          fx_s = fs * Flow % cell_fx(c1) + (1.0 - fs) * Flow % cell_fx(c2)
          fy_s = fs * Flow % cell_fy(c1) + (1.0 - fs) * Flow % cell_fy(c2)
          fz_s = fs * Flow % cell_fz(c1) + (1.0 - fs) * Flow % cell_fz(c2)
        else
          fx_s = Flow % face_fx(s)
          fy_s = Flow % face_fy(s)
          fz_s = Flow % face_fz(s)
        end if

        p_res(s) = p % n(c2) - p % n(c1)        &
                 - (  fx_s * Grid % dx(s)       &
                    + fy_s * Grid % dy(s)       &
                    + fz_s * Grid % dz(s) )
      end do

      ! Gradient is body force plus least squares of the residuals
      call Grad_Residual(1, p % x, Flow % cell_fx)
      call Grad_Residual(2, p % y, Flow % cell_fy)
      call Grad_Residual(3, p % z, Flow % cell_fz)

      ! Refresh buffers for gradient components
      call Grid % Exchange_Cells_Real(p % x)
      call Grid % Exchange_Cells_Real(p % y)
      call Grid % Exchange_Cells_Real(p % z)

    end do

    call Work % Disconnect_Real_Face(p_res)

  end if

  contains

  !============================================================================!
    subroutine Grad_Residual(i, phii, cell_fi)
  !----------------------------------------------------------------------------!
  !   Least squares gradient (component i) from differences over faces in      !
  !   p_res, plus the body force.  Mirrors Grad_Component_No_Refresh.          !
  !----------------------------------------------------------------------------!
    implicit none
  !--------------------------------[Arguments]---------------------------------!
    integer, intent(in) :: i
    real                :: phii   (-Grid % n_bnd_cells:Grid % n_cells)
    real                :: cell_fi(-Grid % n_bnd_cells:Grid % n_cells)
  !----------------------------------[Locals]----------------------------------!
    real, contiguous, pointer :: g(:,:)
    integer, dimension(3,3), parameter :: MAP = reshape((/ 1, 4, 5,    &
                                                           4, 2, 6,    &
                                                           5, 6, 3 /), &
                                                        shape(MAP))
  !============================================================================!

    g => Flow % grad_c2c

    do c = Cells_In_Domain()
      phii(c) = cell_fi(c)
    end do

    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .ne. SYMMETRY) then
        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1,s)
          phii(c1) = phii(c1)                                  &
                   + p_res(s) * (  g(MAP(i,1),c1) * Grid % dx(s)  &
                                 + g(MAP(i,2),c1) * Grid % dy(s)  &
                                 + g(MAP(i,3),c1) * Grid % dz(s))
        end do
      end if
    end do

    do s = Faces_In_Domain_And_At_Buffers()
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      phii(c1) = phii(c1)                                  &
               + p_res(s) * (  g(MAP(i,1),c1) * Grid % dx(s)  &
                             + g(MAP(i,2),c1) * Grid % dy(s)  &
                             + g(MAP(i,3),c1) * Grid % dz(s))

      phii(c2) = phii(c2)                                  &
               + p_res(s) * (  g(MAP(i,1),c2) * Grid % dx(s)  &
                             + g(MAP(i,2),c2) * Grid % dy(s)  &
                             + g(MAP(i,3),c2) * Grid % dz(s))
    end do

    end subroutine

  end subroutine
