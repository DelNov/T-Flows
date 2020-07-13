!==============================================================================!
  subroutine Multiphase_Mod_Vof_Open_Boundary(flow, mult)
!------------------------------------------------------------------------------!
!   Interpolate values at open boundary. This is donde using                   !
!   https://doi.org/10.1080/10407790.2017.1338077, eq. 26                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  integer                       :: ini    ! current temporal iteration
  real                          :: mass_res
!-----------------------------------[Locals]-----------------------------------!
  type(Matrix_Type),    pointer :: a
  type(Grid_Type),      pointer :: grid
  type(Var_Type),       pointer :: v_comp
  real, contiguous,     pointer :: b(:)
  integer                       :: s, ss, c, c1, c2, c11, c22, i_fac
  integer                       :: i, corr_steps
  real                          :: phi_f, fs, signo, sumx, sumy, sumz
  real                          :: a_open
!==============================================================================!

  grid   => flow % pnt_grid

  ! Browse through velocity components
  a_open = 0.0
  do i = 1, 3

    select case(i)
      case(1)
        v_comp => flow % u
      case(2)
        v_comp => flow % v
      case(3)
        v_comp => flow % w
    end select

    do s = 1, grid % n_bnd_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OPENBC) then
        sumx = 0.0; sumy = 0.0; sumz = 0.0
        if (i == 1) then
          a_open = a_open + grid % s(s)
        end if
        do i_fac = 1, grid % cells_n_faces(c1)
          ss = grid % cells_f(i_fac, c1)
          if (ss .ne. s) then
            fs = grid % fw(ss)
            c11 = grid % faces_c(1,ss)
            c22 = grid % faces_c(2,ss)
            if (c11 == c1) then
              signo = 1.0
              phi_f = fs * v_comp % n(c11) + (1.0 - fs) * v_comp % n(c22)
            else
              signo = -1.0
              phi_f = fs * v_comp % n(c22) + (1.0 - fs) * v_comp % n(c11)
            end if
            sumx = sumx + 1.0 / grid % vol(c1) * signo * phi_f * grid % sx(ss)
            sumy = sumy + 1.0 / grid % vol(c1) * signo * phi_f * grid % sy(ss)
            sumz = sumz + 1.0 / grid % vol(c1) * signo * phi_f * grid % sz(ss)
          end if
        end do
        v_comp % n(c2) = v_comp % n(c1)                       &
                       + dot_product( (/sumx, sumy, sumz/),   &
                                      (/grid % dx(s),         &
                                        grid % dy(s),         &
                                        grid % dz(s)/) )
        v_comp % n(c2) = v_comp % n(c2)                       &
                       / (1.0 - 1.0 / grid % vol(c1)          &
                       * dot_product( (/grid % sx(s),         &
                                        grid % sy(s),         &
                                        grid % sz(s)/),       &
                                      (/grid % dx(s),         &
                                        grid % dy(s),         &
                                        grid % dz(s)/) ) )
      end if
    end do
  end do

  call Comm_Mod_Global_Sum_Real(a_open)

  ! correct values using fluxes
  if (mult % phase_change) then

    ! The total added mass
    ! the sign needs to be negative in case of evaporation, so it acts like a
    ! synthetic inflow boundary

    mult % add_mass_in = 0.0
    do c = 1, grid % n_cells - grid % comm % n_buff_cells
      mult % add_mass_in = mult % add_mass_in - mult % flux_rate(c)            &
                                              * grid % vol(c)                  &
                                              * flow % density(c)              &
                                              * ( 1.0 / mult % phase_dens(2)   &
                                                - 1.0 / mult % phase_dens(1) )
    end do
  else
    mult % add_mass_in = 0.0
    do s = 1, grid % n_bnd_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
        flow % m_flux % n(s) = flow % density_f(s)                  &
                             * dot_product( (/flow % u % n(c2),     &
                                              flow % v % n(c2),     &
                                              flow % w % n(c2)/),   &
                                            (/grid % sx(s),         &
                                              grid % sy(s),         &
                                              grid % sz(s)/) )
        mult % add_mass_in  = mult % add_mass_in + flow % m_flux % n(s)
      end if
    end do
  end if

  mult % add_mass_out = 0.0
  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OPENBC) then
      mult % add_mass_out = mult % add_mass_out                  &
                          + flow % density_f(s)                  &
                          * dot_product( (/flow % u % n(c2),     &
                                           flow % v % n(c2),     &
                                           flow % w % n(c2)/),   &
                                         (/grid % sx(s),         &
                                           grid % sy(s),         &
                                           grid % sz(s)/) )
    end if
  end do
  call Comm_Mod_Global_Sum_Real(mult % add_mass_in)
  call Comm_Mod_Global_Sum_Real(mult % add_mass_out)

  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OPENBC) then
      flow % m_flux % n(s) = flow % m_flux % n(s) - grid % s(s) / a_open   &
                           * (mult % add_mass_in - mult % add_mass_out)
      !write(*,*) flow % m_flux % n(s), mult % add_mass_in, mult % add_mass_out
      !write(*,*) flow % u % n(c2), flow % v % n(c2), flow % w % n(c2)
    end if
  end do

  mult % add_mass_out = 0.0
  do s = 1, grid % n_bnd_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OPENBC) then
      mult % add_mass_out = mult % add_mass_out + flow % m_flux % n(s)
    end if
  end do
  write(*,*) mult % add_mass_in, mult % add_mass_out
  end subroutine
