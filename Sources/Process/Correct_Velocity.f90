!==============================================================================!
  real function Correct_Velocity(grid, dt, ini)
!------------------------------------------------------------------------------!
!   Corrects the velocities, and mass fluxes on the cell faces.                !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Comm_Mod
  use Les_Mod
  use Grid_Mod,     only: Grid_Type
  use Bulk_Mod
  use Info_Mod
  use Control_Mod
  use User_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: dt
  integer         :: ini
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, c1, c2, s, m
  real              :: cfl_max(256), pe_max(256)
  real              :: cfl_t, pe_t, mass_err
  character(len=80) :: coupling
!==============================================================================!

  ! User function
  call User_Mod_Beginning_Of_Correct_Velocity(grid, dt, ini)

  call Control_Mod_Pressure_Momentum_Coupling(coupling)

  !-----------------------------------------!
  !   Correct velocities and fluxes with    !
  !    periodic part of the pressure to     !
  !    obtain divergence free velocity      !
  !- - - - - - - - - - - - - - - - - - - -  !
  !   For SOLIDs, px, py and pz are zero    !
  !   so this loop will not correct SOLID   !
  !   velocities.                           !
  !-----------------------------------------!
  if(coupling .eq. 'PROJECTION') then
    do c = 1, grid % n_cells
      u % n(c) = u % n(c) - p % x(c) * grid % vol(c) / a % sav(c)
      v % n(c) = v % n(c) - p % y(c) * grid % vol(c) / a % sav(c)
      w % n(c) = w % n(c) - p % z(c) * grid % vol(c) / a % sav(c)
    end do 
  else ! coupling is 'SIMPLE'
    do c = 1, grid % n_cells
      u % n(c) = u % n(c) - p % x(c) * grid % vol(c) / a % sav(c)
      v % n(c) = v % n(c) - p % y(c) * grid % vol(c) / a % sav(c)
      w % n(c) = w % n(c) - p % z(c) * grid % vol(c) / a % sav(c)
    end do 
  end if

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
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  !
  !   Note that for FLUID-SOLID interaction, FLUX correctin is zero   !
  !   because a % val(a % pos(1,s)) is also zero.                     !  
  !   What will happen with parallel version ... only god knows.      !
  !-------------------------------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2 > 0 .or.  &
       c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. BUFFER) then
      if(c2  > 0) then
        flux(s)=flux(s)+(pp % n(c2) - pp % n(c1))*a % val(a % pos(1,s))
      else 
        flux(s)=flux(s)+(pp % n(c2) - pp % n(c1))*a % bou(c2)
      end if
    end if             !                                          !
  end do               !<---------- this is correction ---------->!

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
    if(c2 > 0 .or.  &
       c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. BUFFER) then
      b(c1)=b(c1)-flux(s)
      if(c2  > 0) b(c2)=b(c2)+flux(s)
    else
      b(c1) = b(c1)-flux(s)
    end if
  end do

  do c = 1, grid % n_cells
    b(c) = b(c) / (grid % vol(c) * density / dt)
  end do

  mass_err=0.0
  do c = 1, grid % n_cells
    mass_err=max(mass_err, abs(b(c)))
  end do
  call Comm_Mod_Global_Max_Real(mass_err)

  !------------------------------!
  !   Calculate the CFL number   !
  !     and the Peclet number    !
  !------------------------------!
  m = 1
  cfl_max(m) = 0.0
  pe_max(m)  = 0.0
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if( (grid % material(c1) .eq. m) .or.  &
        (grid % material(c2) .eq. m) ) then
      if(c2 > 0 .or.   &
         c2 < 0.and.Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. BUFFER) then
        cfl_t = abs( dt * flux(s) / density /      &
                     ( f_coef(s) *                 &
                     (  grid % dx(s)*grid % dx(s)  &
                      + grid % dy(s)*grid % dy(s)  &
                      + grid % dz(s)*grid % dz(s)) ) )
        pe_t  = abs( flux(s) / f_coef(s) / (viscosity / density + TINY) )
        cfl_max(m) = max( cfl_max(m), cfl_t ) 
        pe_max(m)  = max( pe_max(m),  pe_t  ) 
      end if
    end if
  end do
  call Comm_Mod_Global_Max_Real(cfl_max(m))
  call Comm_Mod_Global_Max_Real(pe_max(m))

  call Info_Mod_Iter_Fill_At(1, 2, 'dum', -1, mass_err)
  call Info_Mod_Bulk_Fill(cfl_max(m),          &
                          pe_max(m),           &                            
                          bulk(m) % flux_x,    &
                          bulk(m) % flux_y,    &
                          bulk(m) % flux_z,    &
                          bulk(m) % p_drop_x,  &
                          bulk(m) % p_drop_y,  &
                          bulk(m) % p_drop_z)

  Correct_Velocity = mass_err ! /(velmax+TINY)

  ! User function
  call User_Mod_End_Of_Correct_Velocity(grid, dt, ini)

  end function
