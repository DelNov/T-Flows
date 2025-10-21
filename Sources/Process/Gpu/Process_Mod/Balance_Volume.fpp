!==============================================================================!
  subroutine Balance_Volume(Process, Flow, Grid)
!------------------------------------------------------------------------------!
!>  First part of inserting volume source for pressure-Poisson equation
!------------------------------------------------------------------------------!
  implicit none
!------------------------------------------------------------------------------!
  class(Process_Type)      :: Process
  type(Field_Type), target :: Flow
  type(Grid_Type),  target :: Grid
!-----------------------------------[Locals]-----------------------------------!
  real, contiguous, pointer :: b(:)
  real                      :: vol_in, vol_out, area_in, area_out, fac
  integer                   :: s, c2, c, reg
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  call Profiler % Start('Balance_Volume')

  ! Take some aliases
  ! GPU version doesn't work if you use directly Flow % whatever_variable
  ! These aliases are really needed, not just some gimmick to shorten the code
  b => Flow % Nat % b

  ! Check if you have pressure gradients at hand and then set aliases properly
  Assert(Flow % stores_gradients_of .eq. 'P')

  ! Nullify the volume source
  !$tf-acc loop begin
  do c = Cells_In_Domain_And_Buffers()  ! all present
    b(c) = 0.0
  end do
  !$tf-acc loop end

  !------------------------------------------!
  !                                          !
  !   Calculate / update volume flow rates   !
  !                                          !
  !------------------------------------------!

  !----------------------------------------------------!
  !   Calculate volume fluxes through boundary faces   !
  !----------------------------------------------------!
  do reg = Boundary_Regions()

    if(Grid % region % type(reg) .eq. SYMMETRY) then
      !$tf-acc loop begin
      do s = Faces_In_Region(reg)
        Flow % v_flux % n(s) = 0.0
      end do
      !$tf-acc loop end
    else
      !$tf-acc loop begin
      do s = Faces_In_Region(reg)
        c2 = Grid % faces_c(2,s)  ! boundary cell
        Flow % v_flux % n(s) = Flow % u % n(c2) * Grid % sx(s)  &
                             + Flow % v % n(c2) * Grid % sy(s)  &
                             + Flow % w % n(c2) * Grid % sz(s)
      end do  ! faces
      !$tf-acc loop end
    end if    ! boundary region type

  end do      ! regions

  !-------------------------------------------!
  !   Calcululate inflow and outflow volume   !
  !   flow rates and inlet and outlet areas   !
  !-------------------------------------------!
  vol_in   = 0.0
  vol_out  = 0.0
  area_in  = 0.0
  area_out = 0.0

  do reg = Boundary_Regions()
    if(Grid % region % type(reg) .eq. INFLOW) then

      !$tf-acc loop begin
      do s = Faces_In_Region(reg)
        area_in = area_in + Grid % s(s)
        vol_in  = vol_in  - Flow % v_flux % n(s)
      end do
      !$tf-acc loop end

    end if

    if(Grid % region % type(reg) .eq. PRESSURE .or.  &
       Grid % region % type(reg) .eq. OUTFLOW  .or.  &
       Grid % region % type(reg) .eq. CONVECT) then

      !$tf-acc loop begin
      do s = Faces_In_Region(reg)
        area_out = area_out + Grid % s(s)
        vol_out  = vol_out  + Flow % v_flux % n(s)
      end do
      !$tf-acc loop end

    end if
  end do

  call Global % Sum_Reals(vol_in,   &
                          vol_out,  &
                          area_in,  &
                          area_out)

  !-----------------------------------------------------------------!
  !   If there is volume imbalance in the source for the pressure   !
  !   (correction) equation, pressure will converge poorly, if at   !
  !   all. Thus, correct the outlet fluxes to enforce the balance   !
  !-----------------------------------------------------------------!

  if(.not. Math % Approx_Real(vol_in, vol_out, FEMTO)) then

    ! You should compute the "fac" now ...
    fac = vol_in / (vol_out + TINY)

    ! ... and correct all velocities
    vol_out = 0.0

    do reg = Boundary_Regions()
      if(Grid % region % type(reg) .eq. OUTFLOW   .or.  &
         Grid % region % type(reg) .eq. CONVECT   .or.  &
         Grid % region % type(reg) .eq. PRESSURE) then

        !$tf-acc loop begin
        do s = Faces_In_Region(reg)  ! all present
          c2 = Grid % faces_c(2,s)

          ! Update velocity components ...
          Flow % u % n(c2) = Flow % u % n(c2) * fac
          Flow % v % n(c2) = Flow % v % n(c2) * fac
          Flow % w % n(c2) = Flow % w % n(c2) * fac

          ! ... volume flux itself ...
          Flow % v_flux % n(s) = Flow % u % n(c2) * Grid % sx(s)    &
                               + Flow % v % n(c2) * Grid % sy(s)    &
                               + Flow % w % n(c2) * Grid % sz(s)

          ! ... and bulk volume out
          vol_out = vol_out + Flow % v_flux % n(s)
        end do
        !$tf-acc loop end

      end if
    end do
  end if

  ! Holy mackrele: summ it up over all processors
  call Global % Sum_Real(vol_out)  ! not checked

  call Profiler % Stop('Balance_Volume')

  end subroutine
