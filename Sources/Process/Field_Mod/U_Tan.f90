!==============================================================================!
  real function Field_Mod_U_Tan(flow, s)
!------------------------------------------------------------------------------!
!   Computes tangential velocity component at the cell near a wall.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  integer                  :: s     ! wall face
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  integer                  :: c1, c2
  real                     :: u_tot_sq, u_nor_sq, u_nor
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid

  ! Take cells around the face
  ! (Note that c2 should be less than zero!
  !  Also note that there is no ckecking yet!)
  c1 = grid % faces_c(1, s)
  c2 = grid % faces_c(2, s)

  ! Compute tangential velocity component
  u_tot_sq = flow % u % n(c1) * flow % u % n(c1)  &
           + flow % v % n(c1) * flow % v % n(c1)  &
           + flow % w % n(c1) * flow % w % n(c1)
  u_nor  = ( flow % u % n(c1) * grid % sx(s)     &
           + flow % v % n(c1) * grid % sy(s)     &
           + flow % w % n(c1) * grid % sz(s) )   &
         / grid % s(s)
  u_nor_sq = u_nor**2

  if( u_tot_sq  > u_nor_sq) then
    Field_Mod_U_Tan = sqrt(u_tot_sq - u_nor_sq)
  else
    Field_Mod_U_Tan = TINY
  end if

  end function
