!==============================================================================!
  subroutine Calculate_Mass_Flow_Rate(grid)
!------------------------------------------------------------------------------!
!   Calculate mass flow rate at cell faces based on velocities only.           !
!----------------------------------[Modules]-----------------------------------!
  use Flow_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: s, c1, c2  
  real              :: fs
!==============================================================================!

  !-------------------------------------------------!
  !   Calculate the mass fluxes on the cell faces   !
  !-------------------------------------------------!
  do s=1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    fs = grid % f(s)

    ! Face is inside the domain
    if( c2 > 0 .or.  &
        c2 < 0 .and. Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. BUFFER) then 
 
      ! Extract the "centred" pressure terms from cell velocities
      flux(s) = density *                                                 &
              (   (fs * u % n(c1) + (1.0-fs) * u % n(c2)) * grid % sx(s)  &
                + (fs * v % n(c1) + (1.0-fs) * v % n(c2)) * grid % sy(s)  &
                + (fs * w % n(c1) + (1.0-fs) * w % n(c2)) * grid % sz(s) )
 
    ! Face is on the boundary
    else
      flux(s) = density *                   &
              (   u % n(c1) * grid % sx(s)  &
                + v % n(c1) * grid % sy(s)  &
                + w % n(c1) * grid % sz(s) )
    end if

  end do

  end subroutine
