!==============================================================================!
  real function Correct_Velocity(flow, sol, dt, ini)
!------------------------------------------------------------------------------!
!   Corrects the velocities, and mass fluxes on the cell faces.                !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Field_Mod,    only: Field_Type, viscosity, density
  use Comm_Mod
  use Les_Mod
  use Grid_Mod,     only: Grid_Type
  use Bulk_Mod,     only: Bulk_Type
  use Info_Mod
  use Solver_Mod,   only: Solver_Type
  use Matrix_Mod,   only: Matrix_Type
  use Control_Mod
  use Numerics_Mod
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: flow
  type(Solver_Type), target :: sol
  real                      :: dt
  integer                   :: ini
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),   pointer :: grid
  type(Bulk_Type),   pointer :: bulk
  type(Var_Type),    pointer :: u, v, w, p, pp
  real,              pointer :: flux(:)
  type(Matrix_Type), pointer :: a
  real,              pointer :: b(:)
  integer                    :: c, c1, c2, s
  real                       :: cfl_max, pe_max
  real                       :: cfl_t, pe_t, mass_err
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  bulk => flow % bulk
  flux => flow % flux
  u    => flow % u
  v    => flow % v
  w    => flow % w
  p    => flow % p
  pp   => flow % pp
  a    => sol % a
  b    => sol % b % val

  ! User function
  call User_Mod_Beginning_Of_Correct_Velocity(flow, dt, ini)

  !-----------------------------------------!
  !   Correct velocities and fluxes with    !
  !    periodic part of the pressure to     !
  !    obtain divergence free velocity      !
  !- - - - - - - - - - - - - - - - - - - -  !
  !   For SOLIDs, px, py and pz are zero    !
  !   so this loop will not correct SOLID   !
  !   velocities.                           !
  !-----------------------------------------!
  do c = 1, grid % n_cells
    u % n(c) = u % n(c) - pp % x(c) * grid % vol(c) / a % sav(c)
    v % n(c) = v % n(c) - pp % y(c) * grid % vol(c) / a % sav(c)
    w % n(c) = w % n(c) - pp % z(c) * grid % vol(c) / a % sav(c)
  end do

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

  !-------------------------------------------------------------------!
  !   Look at the following equation and you will understand why      !
  !   is the matrix for pressure corrections in SIMPLE algorythm      !
  !   formed from the coefficients of the velocity matrix.            !
  !   Moreover, it should also be clear that pressure correction      !
  !   matrix must be formed from underrelaxed velocity coefficients   !
  !-------------------------------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 > 0) then
      flux(s) = flux(s) + (pp % n(c2) - pp % n(c1))*a % val(a % pos(1,s))
    end if               !                                               !
  end do                 !<------------ this is correction ------------->!

  !-------------------------------------!
  !    Calculate the max mass error     !
  !   with the new (corrected) fluxes   !
  !-------------------------------------!
  do c = 1, grid % n_cells
    b(c) = 0.0 
  end do

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    b(c1) = b(c1) - flux(s)
    if(c2 > 0) then
      b(c2) = b(c2) + flux(s)
    end if
  end do

  do c = 1, grid % n_cells
    b(c) = b(c) / (grid % vol(c) * density / dt)
  end do

  mass_err = 0.0
  do c = 1, grid % n_cells - grid % comm % n_buff_cells
    mass_err = max(mass_err, abs(b(c)))
  end do
  call Comm_Mod_Global_Max_Real(mass_err)

  !------------------------------!
  !   Calculate the CFL number   !
  !     and the Peclet number    !
  !------------------------------!
  cfl_max = 0.0
  pe_max  = 0.0
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 > 0) then
      cfl_t = abs( dt * flux(s) / density /      &
                   ( a % fc(s) *                 &
                   (  grid % dx(s)*grid % dx(s)  &
                    + grid % dy(s)*grid % dy(s)  &
                    + grid % dz(s)*grid % dz(s)) ) )
      pe_t    = abs( flux(s) / a % fc(s) / (viscosity / density + TINY) )
      cfl_max = max( cfl_max, cfl_t ) 
      pe_max  = max( pe_max,  pe_t  ) 
    end if
  end do
  call Comm_Mod_Global_Max_Real(cfl_max)
  call Comm_Mod_Global_Max_Real(pe_max)

  call Info_Mod_Iter_Fill_At(1, 5, 'dum', -1, mass_err)
  call Info_Mod_Bulk_Fill(cfl_max,          &
                          pe_max,           &
                          bulk % flux_x,    &
                          bulk % flux_y,    &
                          bulk % flux_z,    &
                          bulk % p_drop_x,  &
                          bulk % p_drop_y,  &
                          bulk % p_drop_z)

  Correct_Velocity = mass_err ! /(velmax+TINY)

  ! User function
  call User_Mod_End_Of_Correct_Velocity(flow, dt, ini)

  end function
