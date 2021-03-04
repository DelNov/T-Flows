!==============================================================================!
  subroutine Correct_Velocity(flow, turb, mult, sol, ini)
!------------------------------------------------------------------------------!
!   Corrects the velocities, and mass (or volume) fluxes on cell faces.        !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Solver_Type),     target :: sol
  integer                       :: ini
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Bulk_Type),   pointer :: bulk
  type(Var_Type),    pointer :: u, v, w, p, pp
  type(Face_Type),   pointer :: v_flux          ! volume flux
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  real,              pointer :: u_relax
  integer                    :: c, c1, c2, s, i_fac
  real                       :: cfl_t, pe_t, dens_f, visc_f, dt, wght
!==============================================================================!

  call Cpu_Timer_Mod_Start('Correct_Velocity')

  ! Take aliases
  grid    => flow % pnt_grid
  bulk    => flow % bulk
  v_flux  => flow % v_flux
  p       => flow % p
  pp      => flow % pp
  dt      =  flow % dt
  a       => sol % a
  b       => sol % b % val
  u_relax => flow % u_rel_corr

  call Field_Mod_Alias_Momentum(flow, u, v, w)

  ! User function
  call User_Mod_Beginning_Of_Correct_Velocity(flow, mult, sol, ini)

  !-----------------------------------------!
  !   Correct velocities and fluxes with    !
  !    periodic part of the pressure to     !
  !    obtain divergence free velocity      !
  !-----------------------------------------!

  ! Normal correction
  if(.not. flow % mass_transfer) then
    do c = 1, grid % n_cells
      u % n(c) = u % n(c) - u_relax * pp % x(c) * grid % vol(c) / a % sav(c)
      v % n(c) = v % n(c) - u_relax * pp % y(c) * grid % vol(c) / a % sav(c)
      w % n(c) = w % n(c) - u_relax * pp % z(c) * grid % vol(c) / a % sav(c)
    end do

  ! In case with mass transfer, dodge cells at the surface
  else
    do c = 1, grid % n_cells
      if(mult % cell_at_elem(c) .eq. 0) then
        u % n(c) = u % n(c) - u_relax * pp % x(c) * grid % vol(c) / a % sav(c)
        v % n(c) = v % n(c) - u_relax * pp % y(c) * grid % vol(c) / a % sav(c)
        w % n(c) = w % n(c) - u_relax * pp % z(c) * grid % vol(c) / a % sav(c)
      end if
    end do

    call Multiphase_Mod_Vof_Grad_Variable_With_Jump(mult, flow % u)
    call Multiphase_Mod_Vof_Grad_Variable_With_Jump(mult, flow % v)
    call Multiphase_Mod_Vof_Grad_Variable_With_Jump(mult, flow % w)

    do c = 1, grid % n_cells
      if(mult % cell_at_elem(c) .ne. 0) then

        u % n(c) = 0.0
        v % n(c) = 0.0
        w % n(c) = 0.0

        ! Browse through faces of the current cell
        do i_fac = 1, grid % cells_n_faces(c)
          s = grid % cells_f(i_fac, c)
          c1 = grid % faces_c(1,s)
          c2 = grid % faces_c(2,s)

          if(c2 > 0) then

            ! Cell 2 is at the front
            if(mult % cell_at_elem(c2) .ne. 0 .and.  &
               mult % cell_at_elem(c1) .eq. 0) then
              wght = (1.0 - abs(mult % vof % n(c1) - mult % vof % n(c2)))
              u % n(c2) = u % n(c2) + wght * u % n(c1)
              v % n(c2) = v % n(c2) + wght * v % n(c1)
              w % n(c2) = w % n(c2) + wght * w % n(c1)
            end if

            ! Cell 1 is at the front
            if(mult % cell_at_elem(c1) .ne. 0 .and.  &
               mult % cell_at_elem(c2) .eq. 0) then
              wght = (1.0 - abs(mult % vof % n(c1) - mult % vof % n(c2)))
              u % n(c1) = u % n(c1) + wght * u % n(c2)
              v % n(c1) = v % n(c1) + wght * v % n(c2)
              w % n(c1) = w % n(c1) + wght * w % n(c2)
            end if

          end if

        end do
      end if  ! cell at elem
    end do    ! cells
  end if      ! if flow % mass_transfer

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2  < 0) then
      if( (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE) ) then
        u % n(c2) = u % n(c1)
        v % n(c2) = v % n(c1)
        w % n(c2) = w % n(c1)
      end if
    end if
  end do

  !----------------------------------------------------------!
  !   Compute velocity gradients, taking jump into account   !
  !----------------------------------------------------------!

  call Update_Boundary_Values(flow, turb, mult)

  ! If there is a jump in velocities, call specialized gradient calculation
  if(mult % model .eq. VOLUME_OF_FLUID .and. flow % mass_transfer) then
    call Multiphase_Mod_Vof_Grad_Variable_With_Jump(mult, flow % u)
    call Multiphase_Mod_Vof_Grad_Variable_With_Jump(mult, flow % v)
    call Multiphase_Mod_Vof_Grad_Variable_With_Jump(mult, flow % w)

  ! No jumps, call usual routines
  else
    call Field_Mod_Grad_Variable(flow, flow % u)
    call Field_Mod_Grad_Variable(flow, flow % v)
    call Field_Mod_Grad_Variable(flow, flow % w)
  end if

  !----------------------------------------------------------------!
  !   Look at the following equation and you will understand why   !
  !   is the matrix for pressure corrections in SIMPLE algorithm   !
  !   formed from the coefficients of the velocity matrix.         !
  !   pp      [kg/ms^2]                                            !
  !   a % val [m^4s/kg]                                            !
  !----------------------------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 > 0) then
                                     !<--- this is correction --->!
      v_flux % n(s) = v_flux % n(s) + ( pp % n(c2) - pp % n(c1) )   &
                                       * a % val(a % pos(1,s))
    end if
  end do

  !------------------------------!
  !   Calculate the CFL number   !
  !     and the Peclet number    !
  !------------------------------!
  flow % cfl_max = 0.0
  flow % pe_max  = 0.0
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    dens_f =        grid % fw(s)  * flow % density  (c1)   &
           + (1.0 - grid % fw(s)) * flow % density  (c2)
    visc_f =        grid % fw(s)  * flow % viscosity(c1)   &
           + (1.0 - grid % fw(s)) * flow % viscosity(c2)
    if(c2 > 0) then
      cfl_t = abs( dt * v_flux % n(s) /          &
                   ( a % fc(s) *                 &
                   (  grid % dx(s)*grid % dx(s)  &
                    + grid % dy(s)*grid % dy(s)  &
                    + grid % dz(s)*grid % dz(s)) ) )
      pe_t    = abs( v_flux % n(s) / a % fc(s) / (visc_f / dens_f + TINY) )
      flow % cfl_max = max( flow % cfl_max, cfl_t )
      flow % pe_max  = max( flow % pe_max,  pe_t  )
    end if
  end do
  call Comm_Mod_Global_Max_Real(flow % cfl_max)
  call Comm_Mod_Global_Max_Real(flow % pe_max)

  if (flow % p_m_coupling == SIMPLE) then
    call Info_Mod_Iter_Fill_At(1, 5, 'dum', -1, flow % vol_res)
  else
    if (flow % i_corr == flow % n_piso_corrections) then
      call Info_Mod_Iter_Fill_At(1, 5, 'dum', -1, flow % vol_res)
    end if
  end if

  ! User function
  call User_Mod_End_Of_Correct_Velocity(flow, mult, sol, ini)

  call Cpu_Timer_Mod_Stop('Correct_Velocity')

  end subroutine
