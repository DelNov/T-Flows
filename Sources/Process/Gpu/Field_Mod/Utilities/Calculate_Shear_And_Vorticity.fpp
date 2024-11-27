!==============================================================================!
  subroutine Calculate_Shear_And_Vorticity(Flow, Grid)
!------------------------------------------------------------------------------!
!   Computes the magnitude of the Flow % shear stress.                         !
!------------------------------------------------------------------------------!
!   Shear of the velocity vield is a tensor define as:                         !
!                                                                              !
!   s_ij  = 1/2 ( dui/dxj + duj/dxi )                                          !
!                                                                              !
!   Shear's magnitude is computed as:                                          !
!                                                                              !
!   Flow % shear = sqrt( 2 * s_ij * s_ij )                                     !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !
!   Rotation of a velocity field is defined as:                                !
!                                                                              !
!         1         1                                                          !
!   rot = - ∇ x u = - ((dw/dy-dv/dz) i + (du/dz-dw/dx) j + (dv/dx-du/dy) k)    !
!         2         2                                                          !
!                                                                              !
!   Vorticity is twice the rotation vector, hence                              !
!                                                                              !
!                                                                              !
!   Flow % vort = ∇ x u = (dw/dy-dv/dz) i + (du/dz-dw/dx) j + (dv/dx-du/dy) k  !
!                                                                              !
!                = v_x i + v_y j + v_z k                                       !
!                                                                              !
!   Vorticity magnitude would be the magnitude of Flow % vorticity vector      !
!                                                                              !
!   |Flow % vort| = sqrt(v_x^2 + v_y^2 + v_z^2)                                !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Field_Type), target :: Flow
  type(Grid_Type),   target :: Grid  !! grid object is need for GPU
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: u_x(:), u_y(:), u_z(:),  &
                               v_x(:), v_y(:), v_z(:),  &
                               w_x(:), w_y(:), w_z(:)
  integer                   :: c, c1, c2, s, reg, run
!==============================================================================!

  !----------------------------------------!
  !   Copy values to symmetry boundaries   !
  !   (Probably not the most consistent)   !
  !----------------------------------------!
  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. SYMMETRY) then

      !$tf-acc loop begin
      do s = Faces_In_Region(reg)  ! all present
        c1 = Grid % faces_c(1,s)   ! inside cell
        c2 = Grid % faces_c(2,s)   ! boundary cell
        Flow % u % n(c2) = Flow % u % n(c1)
        Flow % v % n(c2) = Flow % v % n(c1)
        Flow % w % n(c2) = Flow % w % n(c1)
      end do
      !$tf-acc loop end

    end if
  end do

  !----------------------------------------!
  !   Compute Flow % shear in three runs   !
  !----------------------------------------!
  do run = 1, 3

    !---------------------------!
    !   Run #1: u_x, w_y, v_z   !
    !---------------------------!
    if(run .eq. 1) then
      u_x => Flow % phi_x
      w_y => Flow % phi_y
      v_z => Flow % phi_z
      call Flow % Grad_Component(Grid, Flow % u % n, 1, u_x)
      call Flow % Grad_Component(Grid, Flow % w % n, 2, w_y)
      call Flow % Grad_Component(Grid, Flow % v % n, 3, v_z)

      !$tf-acc loop begin
      do c = Cells_In_Domain()  ! all present
        Flow % shear(c) = u_x(c)**2 + .5 * (v_z(c)+w_y(c))**2
        Flow % vort (c) =           - .5 * (v_z(c)-w_y(c))**2
      end do
      !$tf-acc loop end

    !---------------------------!
    !   Run #2: w_x, v_y, u_z   !
    !---------------------------!
    else if(run .eq. 2) then
      w_x => Flow % phi_x
      v_y => Flow % phi_y
      u_z => Flow % phi_z
      call Flow % Grad_Component(Grid, Flow % w % n, 1, w_x)
      call Flow % Grad_Component(Grid, Flow % v % n, 2, v_y)
      call Flow % Grad_Component(Grid, Flow % u % n, 3, u_z)

      !$tf-acc loop begin
      do c = Cells_In_Domain()  ! all present
        Flow % shear(c) = Flow % shear(c) + v_y(c)**2 + .5 * (u_z(c)+w_x(c))**2
        Flow % vort (c) = Flow % vort (c)             - .5 * (u_z(c)-w_x(c))**2
      end do
      !$tf-acc loop end

    !---------------------------!
    !   Run #3: v_x, u_y, w_z   !
    !---------------------------!
    else if(run .eq. 3) then
      v_x => Flow % phi_x
      u_y => Flow % phi_y
      w_z => Flow % phi_z
      call Flow % Grad_Component(Grid, Flow % v % n, 1, v_x)
      call Flow % Grad_Component(Grid, Flow % u % n, 2, u_y)
      call Flow % Grad_Component(Grid, Flow % w % n, 3, w_z)

      !$tf-acc loop begin
      do c = Cells_In_Domain()  ! all present
        Flow % shear(c) = Flow % shear(c) + w_z(c)**2 + .5 * (v_x(c)+u_y(c))**2
        Flow % vort (c) = Flow % vort (c)             - .5 * (v_x(c)-u_y(c))**2
      end do
      !$tf-acc loop end

    end if

  end do

  !$tf-acc loop begin
  do c = Cells_In_Domain()  ! all present
    Flow % shear(c) = sqrt(2.0 * Flow % shear(c))
    Flow % vort (c) = sqrt(2.0 * abs(Flow % vort(c)))
  end do
  !$tf-acc loop end


# if T_FLOWS_DEBUG == 1
    call Grid % Save_Debug_Vtu("shear",                  &
                               scalar_name="shear"    ,  &
                               scalar_cell=Flow % shear)
    call Grid % Save_Debug_Vtu("vorticity",              &
                               scalar_name="vorticity",  &
                               scalar_cell=Flow % vort)
# endif

  end subroutine
