!==============================================================================!
  subroutine Balance_Volume(Process, Flow, Vof)
!------------------------------------------------------------------------------!
!   Modifies fluxes and velocities at outflows to conserve the volume.         !
!   It is called from Compute_Pressure because there we need strict            !
!   volume conservation to achieve convergence of pressure linear solver.      !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
!  use Const_Mod
!  use Comm_Mod
!  use Grid_Mod
!  use Field_Mod
!  use Var_Mod,   only: Var_Type
!  use Face_Mod,  only: Face_Type
!  use Bulk_Mod,  only: Bulk_Type
!  use Math_Mod
!  use Vof_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Process_Type)         :: Process
  type(Field_Type),    target :: Flow
  type(Vof_Type),      target :: Vof
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w
  type(Face_Type), pointer :: v_flux
  integer                  :: s, c1, c2
  real                     :: fac, vol_outflow, area_outflow
!==============================================================================!

  ! Take aliases
  Grid   => Flow % pnt_grid
  bulk   => Flow % bulk
  v_flux => Flow % v_flux
  call Flow % Alias_Momentum(u, v, w)

  !--------------------------------------------------------!
  !                                                        !
  !   Update fluxes at boundaries with latest velocities   !
  !     (These might not obey the volume conservation)     !
  !                                                        !
  !--------------------------------------------------------!
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)

    ! Side is on the boundary
    if(c2 < 0) then

      if(Grid % Bnd_Cond_Type(c2) .eq. SYMMETRY) then
        v_flux % n(s) = 0.0
      else
        v_flux % n(s) = ( u % n(c2) * Grid % sx(s)     &
                        + v % n(c2) * Grid % sy(s)     &
                        + w % n(c2) * Grid % sz(s) )
      end if
    end if

  end do

  !--------------------------------------!
  !                                      !
  !   Added volume due to phase change   !
  !                                      !
  !--------------------------------------!
  if(Flow % mass_transfer) then
    call Vof % Mass_Transfer_Added_Volume(bulk % vol_src)
  end if

  !----------------------------------------------------------!
  !                                                          !
  !   Check if anything comes from real outflow boundaries   !
  !   at all.  For example, at the beginning of simulation   !
  !                                                          !
  !----------------------------------------------------------!
  vol_outflow  = 0.0
  area_outflow = 0.0
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2 < 0 .and. Grid % comm % cell_proc(c1) .eq. this_proc) then
      if(Grid % Bnd_Cond_Type(c2) .eq. PRESSURE .or.  &
         Grid % Bnd_Cond_Type(c2) .eq. OUTFLOW  .or.  &
         Grid % Bnd_Cond_Type(c2) .eq. CONVECT) then

        ! Integrate volume flux only where something does come out
        if(v_flux % n(s) > 0.0) then
          vol_outflow  = vol_outflow  + v_flux % n(s)
        end if

        ! Integrate area everywhere where outflow is specified
        area_outflow = area_outflow + Grid % s(s)
      end if
    end if
  end do
  call Comm_Mod_Global_Sum_Real(vol_outflow)
  call Comm_Mod_Global_Sum_Real(area_outflow)

  !------------------------------------------------------!
  !                                                      !
  !   Calculate all boundary volume fluxes; in and out   !
  !                                                      !
  !------------------------------------------------------!
  bulk % vol_in   = 0.0
  bulk % vol_out  = 0.0
  bulk % area_in  = 0.0
  bulk % area_out = 0.0

  !---------------------------------------------------------------!
  !   This works better if domain features pressure outflow ...   !
  !---------------------------------------------------------------!
  if(Flow % has_pressure) then
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      ! Volume flux at the boundary face
      if(c2 < 0 .and. Grid % comm % cell_proc(c1) .eq. this_proc) then

        if(Grid % Bnd_Cond_Type(c2) .eq. INFLOW   .or.  &
           Grid % Bnd_Cond_Type(c2) .eq. OUTFLOW  .or.  &
           Grid % Bnd_Cond_Type(c2) .eq. CONVECT  .or.  &
           Grid % Bnd_Cond_Type(c2) .eq. PRESSURE) then
          if(v_flux % n(s) > 0.0) then
            bulk % vol_out  = bulk % vol_out  + v_flux % n(s)
            bulk % area_out = bulk % area_out + Grid % s(s)
          else
            bulk % vol_in  = bulk % vol_in  - v_flux % n(s)
            bulk % area_in = bulk % area_in + Grid % s(s)
          end if
        end if

      end if
    end do

  !-----------------------------------------------------------------------!
  !   ... but this variant works better for all other types of outflows   !
  !-----------------------------------------------------------------------!
  else
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)

      ! Volume flux at the boundary face
      if(c2 < 0 .and. Grid % comm % cell_proc(c1) .eq. this_proc) then

        if(Grid % Bnd_Cond_Type(c2) .eq. INFLOW) then
          bulk % vol_in  = bulk % vol_in  - v_flux % n(s)
          bulk % area_in = bulk % area_in + Grid % s(s)
        end if

        if(Grid % Bnd_Cond_Type(c2) .eq. OUTFLOW  .or.  &
           Grid % Bnd_Cond_Type(c2) .eq. CONVECT  .or.  &
           Grid % Bnd_Cond_Type(c2) .eq. PRESSURE) then
          bulk % vol_out  = bulk % vol_out  + v_flux % n(s)
          bulk % area_out = bulk % area_out + Grid % s(s)
        end if

      end if
    end do
  end if

  call Comm_Mod_Global_Sum_Real(bulk % vol_in)
  call Comm_Mod_Global_Sum_Real(bulk % vol_out)
  call Comm_Mod_Global_Sum_Real(bulk % area_in)
  call Comm_Mod_Global_Sum_Real(bulk % area_out)

  ! Avoid divisions by zero for the cases without any fluid motion
  fac = 1.0
  if(bulk % vol_out .gt. FEMTO) fac = bulk % vol_in / (bulk % vol_out)
  bulk % area_out = max(bulk % area_out, FEMTO)

  !----------------------------------------------------------------------!
  !                                                                      !
  !   Correct velocities and volume fluxes at outlet to balance volume   !
  !    (Here you should distinguish two cases: if anything was coming    !
  !     out of real outflows, like in the first time step, or not.)      !
  !                                                                      !
  !----------------------------------------------------------------------!

  !---------------------------------------------------!
  !   Nothing comes from real outflows, no point in   !
  !   scaling with "fac", just add bulk corrections   !
  !---------------------------------------------------!
  if(Math % Approx_Real(vol_outflow, 0.0, tol=FEMTO)) then

    bulk % vol_out = 0.0
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      if(c2 < 0 .and. Grid % comm % cell_proc(c1) .eq. this_proc) then
        if(Grid % Bnd_Cond_Type(c2) .eq. OUTFLOW  .or.  &
           Grid % Bnd_Cond_Type(c2) .eq. CONVECT  .or.  &
           Grid % Bnd_Cond_Type(c2) .eq. PRESSURE) then

          ! Update velocity components ...
          u % n(c2) = (bulk % vol_in + bulk % vol_src) / area_outflow  &
                    * Grid % sx(s) / Grid % s(s)
          v % n(c2) = (bulk % vol_in + bulk % vol_src) / area_outflow  &
                    * Grid % sy(s) / Grid % s(s)
          w % n(c2) = (bulk % vol_in + bulk % vol_src) / area_outflow  &
                    * Grid % sz(s) / Grid % s(s)

          ! ... volume flux itself ...
          v_flux % n(s) = u % n(c2)*Grid % sx(s)    &
                        + v % n(c2)*Grid % sy(s)    &
                        + w % n(c2)*Grid % sz(s)

          ! ... and bulk volume out
          bulk % vol_out = bulk % vol_out + v_flux % n(s)
        end if
      end if
    end do

    ! Holy mackrele: summ it up over all processors
    call Comm_Mod_Global_Sum_Real(bulk % vol_out)

  !-------------------------------------------------!
  !   Something is coming out from real outflows,   !
  !   presumably not in the first iteration of      !
  !   the first time step, do scaling with "fac"    !
  !-------------------------------------------------!
  else

    ! You should compute the "fac" now ...
    fac = bulk % vol_in / (bulk % vol_out + TINY)

    ! ... and correct all velocities
    bulk % vol_out = 0.0
    do s = 1, Grid % n_faces
      c1 = Grid % faces_c(1,s)
      c2 = Grid % faces_c(2,s)
      if(c2 < 0) then
        if(Grid % Bnd_Cond_Type(c2) .eq. OUTFLOW  .or.  &
           Grid % Bnd_Cond_Type(c2) .eq. CONVECT  .or.  &
           Grid % Bnd_Cond_Type(c2) .eq. PRESSURE) then

          ! Update velocity components ...
          u % n(c2) = u % n(c2) * fac + bulk % vol_src / bulk % area_out  &
                    * Grid % sx(s) / Grid % s(s)
          v % n(c2) = v % n(c2) * fac + bulk % vol_src / bulk % area_out  &
                    * Grid % sy(s) / Grid % s(s)
          w % n(c2) = w % n(c2) * fac + bulk % vol_src / bulk % area_out  &
                    * Grid % sz(s) / Grid % s(s)

          ! ... volume flux itself ...
          v_flux % n(s) = u % n(c2)*Grid % sx(s)    &
                        + v % n(c2)*Grid % sy(s)    &
                        + w % n(c2)*Grid % sz(s)

          ! ... and bulk volume out
          bulk % vol_out = bulk % vol_out + v_flux % n(s)
        end if
      end if
    end do

    ! Holy mackrele: summ it up over all processors
    call Comm_Mod_Global_Sum_Real(bulk % vol_out)  ! not checked

  end if

  ! You were modifying the velocity components -> refresh their buffers
  call Grid % Exchange_Cells_Real(u % n)
  call Grid % Exchange_Cells_Real(v % n)
  call Grid % Exchange_Cells_Real(w % n)

  end subroutine
