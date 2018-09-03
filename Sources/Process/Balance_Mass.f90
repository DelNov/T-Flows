!==============================================================================!
  subroutine Balance_Mass(grid)
!------------------------------------------------------------------------------!
!   Modifies the fluxes at outflow boundaries to conserve the mass.            ! 
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Comm_Mod
  use Grid_Mod
  use Bulk_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: m, s, c1, c2
  real    :: fac(256) 
!==============================================================================!

  !--------------------------------------!
  !   Calculate the inflow mass fluxes   !
  !--------------------------------------!
  do m = 1, grid % n_materials
    bulk(m) % mass_in = 0.0
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 < 0) then
        flux(s) = density*( u % n(c2)*grid % sx(s) + &
                            v % n(c2)*grid % sy(s) + &
                            w % n(c2)*grid % sz(s) )
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
          if(grid % material(c1) .eq. m) then 
            bulk(m) % mass_in = bulk(m) % mass_in - flux(s)
          end if
        end if
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE .and. flux(s)<0.) then
          if(grid % material(c1) .eq. m) then
            bulk(m) % mass_in = bulk(m) % mass_in - flux(s)
          end if
        end if
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT .and. flux(s)<0.) then
          if(grid % material(c1) .eq. m) then
            bulk(m) % mass_in = bulk(m) % mass_in - flux(s)
          end if
        end if
      end if
    end do
    call Comm_Mod_Global_Sum_Real(bulk(m) % mass_in)
  end do                    

  !---------------------------------------!
  !   Calculate the outflow mass fluxes   !
  !     then correct it to satisfy the    ! 
  !          overall mass balance         !
  !---------------------------------------!
  do m = 1, grid % n_materials
    bulk(m) % mass_out = 0.0
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2  < 0) then
        if(Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. OUTFLOW) then
          u % n(c2) = u % n(c1)
          v % n(c2) = v % n(c1)
          w % n(c2) = w % n(c1)
          flux(s) = density * ( u % n(c2)*grid % sx(s) +  & 
                                v % n(c2)*grid % sy(s) +  &
                                w % n(c2)*grid % sz(s) )
          if(grid % material(c1) .eq. m) then
            bulk(m) % mass_out = bulk(m) % mass_out + flux(s)
          end if
        end if
        if(Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. CONVECT .and. flux(s)>0.) then
          flux(s) = density * ( u % n(c2)*grid % sx(s) +  & 
                                v % n(c2)*grid % sy(s) +  &
                                w % n(c2)*grid % sz(s) )
          if(grid % material(c1) .eq. m) then
            bulk(m) % mass_out = bulk(m) % mass_out + flux(s)
          end if
        end if

        flux(s) = density * ( u % n(c2)*grid % sx(s) +  & 
                              v % n(c2)*grid % sy(s) +  &
                              w % n(c2)*grid % sz(s) )
        if(Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. PRESSURE .and. flux(s)>0.) then
          if(grid % material(c1) .eq. m) then
            bulk(m) % mass_out = bulk(m) % mass_out + flux(s)
          end if
        end if

      end if
    end do
    call Comm_Mod_Global_Sum_Real(bulk(m) % mass_out)  ! not checked
  end do

  do m = 1, grid % n_materials
    fac(m) = bulk(m) % mass_in/(bulk(m) % mass_out+TINY)
  end do

  do m = 1, grid % n_materials
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2  < 0) then
        if(Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. OUTFLOW  .or.  &
           Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. CONVECT  .or.  &
           Grid_Mod_Bnd_Cond_Type(grid, c2) .eq. PRESSURE) then
          if(grid % material(c1) .eq. m) then
            u % n(c2) = u % n(c2) * fac(m)
            v % n(c2) = v % n(c2) * fac(m)
            w % n(c2) = w % n(c2) * fac(m)
            flux(s) = flux(s)*fac(m) 
            bulk(m) % mass_out = bulk(m) % mass_out + flux(s)
          end if
        end if
      end if
    end do
  end do

  end subroutine
