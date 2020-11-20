!==============================================================================!
  subroutine Multiphase_Mod_Vof_Mass_Transfer(mult, sol)
!------------------------------------------------------------------------------!
!   Computes mass transfer due to phase change                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),      pointer :: grid
  type(Field_Type),     pointer :: flow
  type(Var_Type),       pointer :: vof, t
  type(Matrix_Type),    pointer :: a
  real,    contiguous,  pointer :: b(:)
  real,    contiguous,  pointer :: qci(:)
  real,    contiguous,  pointer :: flux_rate(:)
  integer, contiguous,  pointer :: ic(:)
  integer                       :: s, ss, c, c1, c2, cc1, cc2
  integer                       :: n, i_fac, i_fac2, i_nod
  integer                       :: count_c, int_c, cat1, cat2
  real                          :: signo, fw
  real                          :: dotprod1, dotprod2, dotprod3, dotprod4
  real                          :: tx_f, ty_f, tz_f, t_flux, some_flux
!==============================================================================!

  grid => mult % pnt_grid
  flow => mult % pnt_flow
  vof  => mult % vof
  t    => flow % t
  a    => sol % a
  b    => sol % b % val

  some_flux = 1.0e+02

  if( .not. allocated(mult % qci)) then
    allocate(mult % qci      (-grid % n_bnd_cells:grid % n_cells))
    allocate(mult % ic       (-grid % n_bnd_cells:grid % n_cells))
    allocate(mult % flux_rate(-grid % n_bnd_cells:grid % n_cells))
  end if
  ic         => mult % ic
  qci        => mult % qci
  flux_rate  => mult % flux_rate

  qci = 0.0
  ic  = 0
  flux_rate = 0.0

  ! cvs can be of two categories cat = 0, it is not considered, cat = 1,
  ! belongs to the interface, cat = 2, does not belong to interface but it is
  ! neighbour to a cat = 1 and it is inside the saturated phase

  !call Field_Mod_Grad_Variable(flow, t)

  ! classify cells as CI and CIs using interior faces
  do s = grid % n_bnd_faces + 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(vof % n(c1) > PICO .and. vof % n(c1) < 1.0 - PICO) then
      ic(c1) = 1
      if(vof % n(c2) > PICO .and. vof % n(c2) < 1.0 - PICO) then
        ic(c2) = 1
      else if(vof % n(c2) > 1.0 - PICO) then
        ic(c2) = 2
      end if
    end if

    if(vof % n(c2) > PICO .and. vof % n(c2) < 1.0 - PICO) then
      ic(c2) = 1
      if(vof % n(c1) > PICO .and. vof % n(c1) < 1.0 - PICO) then
        ic(c1) = 1
      else if(vof % n(c1) > 1.0 - PICO) then
        ic(c1) = 2
      end if
    end if

  end do
  ! Check all cells that contain the interface and calculate heat imported

  !! Boundary faces
  !do s = 1, grid % n_bnd_faces
  !  c1 = grid % faces_c(1,s)
  !  c2 = grid % faces_c(2,s)

  !  if(Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. SYMMETRY) then

  !    tx_f = t % x(c1) * grid % sx(s)
  !    ty_f = t % y(c1) * grid % sy(s)
  !    tz_f = t % z(c1) * grid % sz(s)
  !    dotprod1 = -(tx_f + ty_f + tz_f)

  !    tx_f = t % x(c1) * grid % dx(s)
  !    ty_f = t % y(c1) * grid % dy(s)
  !    tz_f = t % z(c1) * grid % dz(s)
  !    dotprod2 = -(tx_f + ty_f + tz_f)

  !    t_flux = dotprod1 + a % fc(s) * (t % n(c2) - t % n(c1) - dotprod2)
  !    qci(c1) = qci(c1) + max(0.0,- mult % phase_cond(2) * t_flux)

  !  end if
  !end do

  !! Interior faces
  !do s = grid % n_bnd_faces + 1, grid % n_faces
  !  c1 = grid % faces_c(1,s)
  !  c2 = grid % faces_c(2,s)
  !  fw = grid % fw(s)

  !  tx_f = (fw * t % x(c1) + (1.0 - fw) * t % x(c2) ) * grid % sx(s)
  !  ty_f = (fw * t % y(c1) + (1.0 - fw) * t % y(c2) ) * grid % sy(s)
  !  tz_f = (fw * t % z(c1) + (1.0 - fw) * t % z(c2) ) * grid % sz(s)
  !  dotprod1 = -(tx_f + ty_f + tz_f)
  !  dotprod3 = (tx_f + ty_f + tz_f)

  !  tx_f = (fw * t % x(c1) + (1.0 - fw) * t % x(c2) ) * grid % dx(s)
  !  ty_f = (fw * t % y(c1) + (1.0 - fw) * t % y(c2) ) * grid % dy(s)
  !  tz_f = (fw * t % z(c1) + (1.0 - fw) * t % z(c2) ) * grid % dz(s)
  !  dotprod2 = -(tx_f + ty_f + tz_f)
  !  dotprod4 = (tx_f + ty_f + tz_f)

  !  t_flux = dotprod1 + a % fc(s) * (t % n(c2) - t % n(c1) - dotprod2)
  !  qci(c1) = qci(c1) + max(0.0,- mult % phase_cond(2) * t_flux)

  !  t_flux = dotprod3 + a % fc(s) * (t % n(c2) - t % n(c1) - dotprod4)
  !  qci(c2) = qci(c2) + max(0.0,- mult % phase_cond(2) * t_flux)

  !end do

  ! Computing flux rate
  do c = 1, grid % n_cells
    if(ic(c) > 0) then
      if(ic(c) == 1) then
        !flux_rate(c) = qci(c) / (flow % latent_heat * grid % vol(c))
        flux_rate(c) = some_flux
        !write(*,*) c, grid % xc(c), vof % n(c)
      end if
    else ! cells category 0
      qci(c) = 0.0
    end if
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, qci)
  call Grid_Mod_Exchange_Cells_Int (grid, ic )

  ! Interior faces
  do s = grid % n_bnd_faces + 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(ic(c1) == 1) then
      if(ic(c2) == 2) then
        flux_rate(c1) = some_flux
        !flux_rate(c1) = flux_rate(c1)                                     &
        !              + qci(c2) / (flow % latent_heat * grid % vol(c2))
      end if
    end if

    if(ic(c2) == 1) then
      if(ic(c1) == 2) then
        flux_rate(c2) = some_flux
        !flux_rate(c2) = flux_rate(c2)                                     &
        !              + qci(c1) / (flow % latent_heat * grid % vol(c1))
      end if
    end if
  end do

  ! Zeroing cells type 2
  do c = 1, grid % n_cells
    if(ic(c) == 2) then
      qci(c) = 0.0
    end if
  end do

  call Grid_Mod_Exchange_Cells_Real(grid, flux_rate)

  end subroutine
