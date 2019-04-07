!==============================================================================!
  subroutine Convective_Outflow(flow, dt)
!------------------------------------------------------------------------------!
!   Extrapoloate variables on the boundaries where needed.                     !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Field_Mod,   only: Field_Type, heat_transfer, density
  use Turb_Mod
  use Grid_Mod,    only: Grid_Type
  use Bulk_Mod,    only: Bulk_Type
  use Grad_Mod
  use Bulk_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  real                     :: dt
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, t
  real,            pointer :: flux(:)
  integer                  :: c1, c2, s
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  bulk => flow % bulk
  flux => flow % flux
  u    => flow % u
  v    => flow % v
  w    => flow % w
  t    => flow % t

  call Bulk_Mod_Calculate_Fluxes(grid, bulk, flux)

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

  if(heat_transfer) then

    ! Temperature gradients might have been computed and
    ! stored already in t % x, t % y and t % z, check it
    call Grad_Mod_Variable(t, .true.)

    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)

      ! On the boundary perform the extrapolation
      if(c2  < 0) then
        if( (Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT) ) then
          t % n(c2) = t % n(c2)   &
                    - ( bulk % u * t % x(c1)         &
                      + bulk % v * t % y(c1)         &
                      + bulk % w * t % z(c1) ) * dt
        end if
      end if
    end do
  end if

  end subroutine
