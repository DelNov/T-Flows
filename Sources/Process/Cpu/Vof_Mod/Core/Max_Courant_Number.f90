!==============================================================================!
  subroutine Max_Courant_Number(Vof, dt, interf, courant_max)
!------------------------------------------------------------------------------!
!   Computes the Maximum Courant Number at cells. The argument interf helps    !
!   selecting if calculation will be performed close the interface, which      !
!   in turn will modify the time step if necessary and perform a simple        !
!   stepping approach. If interf = 0, then calculation is made everywhere      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Vof_Type), target :: Vof
  real                    :: dt            ! time step
  integer, intent(in)     :: interf
  real                    :: courant_max
!--------------------------------[Locals]--------------------------------------!
  type(Field_Type), pointer :: Flow
  type(Grid_Type),  pointer :: Grid
  type(Var_Type),   pointer :: fun
  type(Face_Type),  pointer :: v_flux
  real, contiguous, pointer :: c_d(:)
  integer                   :: c, c1, c2, s, reg
  real                      :: fun_dist
!==============================================================================!

  ! Take aliases
  Flow   => Vof % pnt_flow
  v_flux => Flow % v_flux
  Grid   => Flow % pnt_grid
  fun    => Vof % fun
  c_d    => Vof % c_d

  courant_max = -HUGE

  ! Initialize
  c_d(-Vof % pnt_grid % n_bnd_cells:Vof % pnt_grid % n_cells) = 0.0

  if(interf == 1) then

    ! Interior faces
    do s = Faces_In_Domain_And_At_Buffers()
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      ! For c1
      fun_dist = min(max(fun % n(c1), 0.0), 1.0)
      fun_dist = (1.0 - fun_dist) ** 2 * fun_dist ** 2 * 16.0

      c_d(c1) = c_d(c1) + fun_dist  &
              * max(-v_flux % n(s) * dt / Grid % vol(c1), 0.0)

      ! For c2
      fun_dist = min(max(fun % n(c2), 0.0), 1.0)

      fun_dist = (1.0 - fun_dist) ** 2 * fun_dist ** 2 * 16.0

      c_d(c2) = c_d(c2) + fun_dist  &
              * max( v_flux % n(s) * dt / Grid % vol(c2), 0.0)
    end do

    ! Boundary faces
    do reg = Boundary_Regions()
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)

        ! For c1
        fun_dist = min(max(fun % n(c1), 0.0), 1.0)
        fun_dist = (1.0 - fun_dist) ** 2 * fun_dist ** 2 * 16.0

        c_d(c1) = c_d(c1) + fun_dist  &
                * max(-v_flux % n(s) * dt / Grid % vol(c1), 0.0)
      end do  ! faces
    end do    ! regions

    ! if(Vof % phase_Change) then
    !   do c = Cells_In_Domain_And_Buffers()
    !     fun_dist = min(max(fun % n(c1), 0.0),1.0)
    !     fun_dist = (1.0 - fun_dist) ** 2 * fun_dist ** 2 * 16.0
    !     c_d(c) = c_d(c) + fun_dist * Vof % flux_rate(c)    &
    !                                / Flow % density_f(s) * dt
    !   end do
    ! end if

    call Grid % Exchange_Cells_Real(c_d)

    do c = Cells_In_Domain_And_Buffers()
      courant_max = max(c_d(c), courant_max)
    end do
    call Global % Max_Real(courant_max)

  else  ! interf = 0

    ! Interior faces
    do s = Faces_In_Domain_And_At_Buffers()
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      c_d(c1) = c_d(c1) + max(-v_flux % n(s) * dt / Grid % vol(c1), 0.0)
      c_d(c2) = c_d(c2) + max( v_flux % n(s) * dt / Grid % vol(c2), 0.0)
    end do

    ! Boundary faces
    do reg = Boundary_Regions()
      do s = Faces_In_Region(reg)
        c1 = Grid % faces_c(1,s)
        c_d(c1) = c_d(c1) + max(-v_flux % n(s) * dt / Grid % vol(c1), 0.0)
      end do  ! faces
    end do    ! regions

    ! if(Vof % phase_Change) then
    !   do c = Cells_In_Domain_And_Buffers()
    !     c_d(c) = c_d(c) + Vof % flux_rate(c) / Flow % density_f(s) * dt
    !   end do
    ! end if

    call Grid % Exchange_Cells_Real(c_d)

  end if  ! if interf == 1

  end subroutine
