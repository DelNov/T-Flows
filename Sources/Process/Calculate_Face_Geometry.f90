!==============================================================================!
  subroutine Calculate_Face_Geometry(grid)
!------------------------------------------------------------------------------!
!   Calculates additional geometrical quantities.                              !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Flow_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c1, c2, s
  real    :: xc1, yc1, zc1, xc2, yc2, zc2
!==============================================================================!

  !----------------------------------------------!
  !   Calculate total surface of the cell face   ! 
  !----------------------------------------------!
  do s = 1, grid % n_faces
    f_coef(s) = (  grid % sx(s)*grid % sx(s)  &
                 + grid % sy(s)*grid % sy(s)  &
                 + grid % sz(s)*grid % sz(s) )  
  end do

  !-------------------------------------------------------!
  !   Calculate the distance between neighbouring cells   !
  !-------------------------------------------------------!
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    xc1 = grid % xc(c1)
    yc1 = grid % yc(c1)
    zc1 = grid % zc(c1) 

    xc2 = grid % xc(c2) + grid % dx(s)
    yc2 = grid % yc(c2) + grid % dy(s)
    zc2 = grid % zc(c2) + grid % dz(s)

    grid % dx(s) = xc2-xc1
    grid % dy(s) = yc2-yc1
    grid % dz(s) = zc2-zc1

    f_coef(s) = f_coef(s)                    &
             / (  grid % dx(s)*grid % sx(s)  &
                + grid % dy(s)*grid % sy(s)  &
                + grid % dz(s)*grid % sz(s)) 
  end do  ! faces 

  !--------------------------------------------!
  !   Calculate weight factors for the faces   !
  !--------------------------------------------!
  do s = 1, grid % n_faces    
    c1 = grid % faces_c(1,s) 
    c2 = grid % faces_c(2,s)
                          
    ! Inside the flow, it has usual value: phi_f = fw * phi_1 + (1-fw) * phi_2
    fw(s) = grid % f(s)    
                         
    ! Close to the wall, however, there is inversion. It takes
    ! the value from inside as the representative for the face.
    if(c2 < 0) then
      if (Grid_Mod_Bnd_Cond_Type(grid,c2) .ne. BUFFER) then
        fw(s) = 1.
      end if
    end if           
  end do            


  !------------------------------!
  !   Find the near-wall cells   !
  !------------------------------!
  grid % cell_near_wall = .false.

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
        grid % cell_near_wall(c1) = .true.  
      end if
    end if 

  end do  ! faces 

  end subroutine
