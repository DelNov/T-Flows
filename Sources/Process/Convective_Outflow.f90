!==============================================================================!
  subroutine Convective_Outflow(grid, dt)
!------------------------------------------------------------------------------!
!   Extrapoloate variables on the boundaries where needed.                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Rans_Mod
  use Grid_Mod
  use Grad_Mod
  use Bulk_Mod
  use Control_Mod
  use Work_Mod, only: t_x => r_cell_01,  &
                      t_y => r_cell_02,  &
                      t_z => r_cell_03           
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  real            :: dt
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s
!==============================================================================!

  call Bulk_Mod_Compute_Fluxes(grid, bulk, flux)

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    ! On the boundary perform the extrapolation
    if(c2  < 0) then
      if( (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) ) then
        u % n(c2) = u % n(c2)   &
                  - ( bulk % u * u % x(c1)         & 
                    + bulk % v * u % y(c1)         &
                    + bulk % w * u % z(c1) ) * dt
        v % n(c2) = v % n(c2)  &
                  - ( bulk % u * v % x(c1)         & 
                    + bulk % v * v % y(c1)         &
                    + bulk % w * v % z(c1) ) * dt
        w % n(c2) = w % n(c2)  &
                  - ( bulk % u * w % x(c1)         & 
                    + bulk % v * w % y(c1)         &
                    + bulk % w * w % z(c1) ) * dt
      end if
    end if
  end do

  if(heat_transfer .eq. YES) then
    call Grad_Mod_For_Phi(grid, t % n, 1, t_x, .true.)     ! dT/dx
    call Grad_Mod_For_Phi(grid, t % n, 2, t_y, .true.)     ! dT/dy
    call Grad_Mod_For_Phi(grid, t % n, 3, t_z, .true.)     ! dT/dz
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if(c2  < 0) then
        if( (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) ) then
          t % n(c2) = t % n(c2)   &
                    - ( bulk % u * t_x(c1)        & 
                      + bulk % v * t_y(c1)        &
                      + bulk % w * t_z(c1) ) * dt
        end if
      end if
    end do
  end if

  end subroutine
