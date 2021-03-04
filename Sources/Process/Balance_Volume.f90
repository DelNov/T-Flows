!==============================================================================!
  subroutine Balance_Volume(flow, mult)
!------------------------------------------------------------------------------!
!   Modifies the fluxes at outflow boundaries to conserve the volume.          !
!   This function modifies velocities and volume fluxes at outflows.           !
!   It should be called with fresh values at boundaries, thus after the        !
!   call to Update_Boundary_Values.                                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Comm_Mod
  use Grid_Mod,        only: Grid_Type, Grid_Mod_Bnd_Cond_Type,  &
                             INFLOW, OUTFLOW, CONVECT, PRESSURE, SYMMETRY
  use Field_Mod,       only: Field_Type, Field_Mod_Alias_Momentum
  use Var_Mod,         only: Var_Type
  use Face_Mod,        only: Face_Type
  use Bulk_Mod,        only: Bulk_Type
  use Multiphase_Mod,  only: Multiphase_Type,  &
                             Multiphase_Mod_Vof_Mass_Transfer_Added_Volume
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Multiphase_Type), target :: mult
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w
  type(Face_Type), pointer :: v_flux
  integer                  :: s, e, c, c1, c2
  real                     :: fac, area_in, area_out
!==============================================================================!

  ! Take aliases
  grid   => flow % pnt_grid
  bulk   => flow % bulk
  v_flux => flow % v_flux
  call Field_Mod_Alias_Momentum(flow, u, v, w)

  !----------------------------------!
  !   Calculate all boundary fluxe   !
  !    and the inflow volume flux    !
  !----------------------------------!
  bulk % vol_in = 0.0
  area_in       = 0.0
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 < 0 .and. grid % comm % cell_proc(c1) .eq. this_proc) then
      v_flux % n(s) = u % n(c2) * grid % sx(s)    &
                    + v % n(c2) * grid % sy(s)    &
                    + w % n(c2) * grid % sz(s)

      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
        bulk % vol_in = bulk % vol_in - v_flux % n(s)
        area_in = area_in + grid % s(s)
      end if

      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE  &
         .and. v_flux % n(s) < 0.0) then
        bulk % vol_in = bulk % vol_in - v_flux % n(s)
        area_in = area_in + grid % s(s)
      end if

      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT  &
         .and. v_flux % n(s) < 0.0) then
        bulk % vol_in = bulk % vol_in - v_flux % n(s)
        area_in = area_in + grid % s(s)
      end if

      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. SYMMETRY) then
        v_flux % n(s) = 0.0
      end if

    end if
  end do

  call Comm_Mod_Global_Sum_Real(bulk % vol_in)  ! not checked
  call Comm_Mod_Global_Sum_Real(area_in)

  !-------------------------------------------!
  !   Additional volume due to phase change   !
  !-------------------------------------------!
  if(flow % mass_transfer) then
    call Multiphase_Mod_Vof_Mass_Transfer_Added_Volume(mult, bulk % vol_src)
  end if

  !---------------------------------------!
  !   Calculate the outflow mass fluxes   !
  !     then correct it to satisfy the    ! 
  !          overall mass balance         !
  !---------------------------------------!
  bulk % vol_out = 0.0
  area_out = 0.0
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 < 0 .and. grid % comm % cell_proc(c1) .eq. this_proc) then

      if(Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. OUTFLOW) then
        bulk % vol_out = bulk % vol_out + v_flux % n(s)
        area_out = area_out + grid % s(s)
      end if

      if(Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. CONVECT  &
         .and. v_flux % n(s) > 0.0) then
        bulk % vol_out = bulk % vol_out + v_flux % n(s)
        area_out = area_out + grid % s(s)
      end if

      if(Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. PRESSURE  &
         .and. v_flux % n(s) >0.0) then
        bulk % vol_out = bulk % vol_out + v_flux % n(s)
        area_out = area_out + grid % s(s)
      end if

    end if
  end do

  call Comm_Mod_Global_Sum_Real(bulk % vol_out)  ! not checked
  call Comm_Mod_Global_Sum_Real(area_out)

  ! Avoid divisions by zero for the cases without any fluid motion
  fac = 1.0
  if(bulk % vol_out .gt. FEMTO) fac = bulk % vol_in / (bulk % vol_out)
  area_out = max(area_out, FEMTO)

  !----------------------------------------------------------------------!
  !   Correct velocities and volume fluxes at outlet to balance volume   !
  !----------------------------------------------------------------------!

  fac = bulk % vol_in / (bulk % vol_out + TINY)

  bulk % vol_out = 0.0
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2  < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. OUTFLOW  .or.  &
         Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. CONVECT  .or.  &
         Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. PRESSURE) then
        u % n(c2) = u % n(c2) * fac  &
                  + bulk % vol_src / area_out * grid % sx(s) / grid % s(s)
        v % n(c2) = v % n(c2) * fac  &
                  + bulk % vol_src / area_out * grid % sy(s) / grid % s(s)
        w % n(c2) = w % n(c2) * fac  &
                  + bulk % vol_src / area_out * grid % sz(s) / grid % s(s)
        v_flux % n(s) = u % n(c2)*grid % sx(s)    &
                      + v % n(c2)*grid % sy(s)    &
                      + w % n(c2)*grid % sz(s)
        bulk % vol_out = bulk % vol_out + v_flux % n(s)
      end if
    end if
  end do

  end subroutine
