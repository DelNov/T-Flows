!==============================================================================!
  subroutine Balance_Mass(flow, mult)
!------------------------------------------------------------------------------!
!   Modifies the fluxes at outflow boundaries to conserve the mass.            ! 
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Comm_Mod
  use Grid_Mod,        only: Grid_Type, Grid_Mod_Bnd_Cond_Type,  &
                             INFLOW, OUTFLOW, CONVECT, PRESSURE
  use Field_Mod,       only: Field_Type, Field_Mod_Alias_Momentum
  use Var_Mod,         only: Var_Type
  use Face_Mod,        only: Face_Type
  use Bulk_Mod,        only: Bulk_Type
  use Multiphase_Mod,  only: Multiphase_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Multiphase_Type), target :: mult
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),      pointer :: grid
  type(Bulk_Type),      pointer :: bulk
  type(Var_Type),       pointer :: u, v, w
  type(Face_Type),      pointer :: m_flux
  integer                       :: s, c, c1, c2
  real                          :: fac
!==============================================================================!

  ! Take aliases
  grid   => flow % pnt_grid
  bulk   => flow % bulk
  m_flux => flow % m_flux
  call Field_Mod_Alias_Momentum(flow, u, v, w)

  !--------------------------------------!
  !   Calculate the inflow mass fluxes   !
  !--------------------------------------!
  bulk % mass_in = 0.0
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 < 0) then
      m_flux % n(s) = flow % density_f(s) * ( u % n(c2)*grid % sx(s)    &
                                            + v % n(c2)*grid % sy(s)    &
                                            + w % n(c2)*grid % sz(s) )

      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
        bulk % mass_in = bulk % mass_in - m_flux % n(s)
      end if

      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE  &
         .and. m_flux % n(s) < 0.0) then
        bulk % mass_in = bulk % mass_in - m_flux % n(s)
      end if

      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT  &
         .and. m_flux % n(s) < 0.0) then
        bulk % mass_in = bulk % mass_in - m_flux % n(s)
      end if
    end if
  end do

  ! Mass source:
  if(mult % phase_change) then
    do c = 1, grid % n_cells - grid % comm % n_buff_cells
      bulk % mass_in = bulk % mass_in                                      &
                     + mult % flux_rate(c) * grid % vol(c)                 &
                                           * flow % density(c)             &
                                           * ( 1.0 / mult % phase_dens(1)  &
                                             - 1.0 / mult % phase_dens(2) )
    end do
  end if

  !do c = 1, grid % n_cells
  !  if (grid % xc(c) > 0.0001 .and. grid % xc(c) < 0.00012) then
  !  !  bulk % mass_in = bulk % mass_in + mult % flux_rate(c)
  !  !if (grid % xc(c) > 3.87 .and. grid % xc(c) < 3.92) then
  !    bulk % mass_in = bulk % mass_in + 1e-10
  !  end if
  !end do
  call Comm_Mod_Global_Sum_Real(bulk % mass_in)

  !---------------------------------------!
  !   Calculate the outflow mass fluxes   !
  !     then correct it to satisfy the    ! 
  !          overall mass balance         !
  !---------------------------------------!
  bulk % mass_out = 0.0
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2  < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. OUTFLOW) then
        u % n(c2) = u % n(c1)
        v % n(c2) = v % n(c1)
        w % n(c2) = w % n(c1)
        m_flux % n(s) = flow % density_f(s) * ( u % n(c2)*grid % sx(s)    &
                                              + v % n(c2)*grid % sy(s)    &
                                              + w % n(c2)*grid % sz(s) )  
        bulk % mass_out = bulk % mass_out + m_flux % n(s)
      end if

      if(Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. CONVECT  &
         .and. m_flux % n(s) > 0.0) then
        m_flux % n(s) = flow % density_f(s) * ( u % n(c2)*grid % sx(s)    &
                                              + v % n(c2)*grid % sy(s)    &
                                              + w % n(c2)*grid % sz(s) )
        bulk % mass_out = bulk % mass_out + m_flux % n(s)
      end if

      m_flux % n(s) = flow % density_f(s) * ( u % n(c2)*grid % sx(s)    &
                                            + v % n(c2)*grid % sy(s)    &
                                            + w % n(c2)*grid % sz(s) )

      if(Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. PRESSURE  &
         .and. m_flux % n(s) >0.0) then
        bulk % mass_out = bulk % mass_out + m_flux % n(s)
      end if

    end if
  end do
  call Comm_Mod_Global_Sum_Real(bulk % mass_out)  ! not checked

  fac = bulk % mass_in / (bulk % mass_out + TINY)
  !write(*,*) fac, bulk % mass_in , bulk % mass_out
  bulk % mass_out = 0.0
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2  < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. OUTFLOW  .or.  &
         Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. CONVECT  .or.  &
         Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. PRESSURE) then
        u % n(c2) = u % n(c2) * fac
        v % n(c2) = v % n(c2) * fac
        w % n(c2) = w % n(c2) * fac
        m_flux % n(s) = m_flux % n(s) * fac
        bulk % mass_out = bulk % mass_out + m_flux % n(s)
      end if
    end if
  end do

  if(mult % phase_change) then
    do s = 1, grid % n_bnd_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. OUTFLOW) then
        m_flux % n(s) = m_flux % n(s) / flow % density_f(s)
      end if
    end do
  end if
  end subroutine
