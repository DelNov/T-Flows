!==============================================================================!
  subroutine Calculate_Mean(grid, n0, n1)   
!------------------------------------------------------------------------------!
!   Calculates time averaged velocity and velocity fluctuations.               !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Les_Mod
  use Rans_Mod
  use Grid_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  integer         :: n0, n1
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, n
!==============================================================================!

  if(turbulence_statistics .eq. NO) return

  n = n1 - n0

  if(n > -1) then

    do c = -grid % n_bnd_cells, grid % n_cells

      !---------------------------------!
      !   Scale-resolving simulations   ! 
      !---------------------------------!
      if(turbulence_model .eq. SMAGORINSKY .or.  &
         turbulence_model .eq. DYNAMIC     .or.  &
         turbulence_model .eq. WALE        .or.  &
         turbulence_model .eq. DES_SPALART .or.  &
         turbulence_model .eq. DNS) then

        u % mean(c) = (u % mean(c) * (1.*n) + u % n(c)) / (1.*(n+1))
        v % mean(c) = (v % mean(c) * (1.*n) + v % n(c)) / (1.*(n+1))
        w % mean(c) = (w % mean(c) * (1.*n) + w % n(c)) / (1.*(n+1))
        p % mean(c) = (p % mean(c) * (1.*n) + p % n(c)) / (1.*(n+1))

        uu % mean(c) = (uu % mean(c)*(1.*n) + u % n(c) * u % n(c)) / (1.*(n+1))
        vv % mean(c) = (vv % mean(c)*(1.*n) + v % n(c) * v % n(c)) / (1.*(n+1))
        ww % mean(c) = (ww % mean(c)*(1.*n) + w % n(c) * w % n(c)) / (1.*(n+1))

        uv % mean(c) = (uv % mean(c)*(1.*n) + u % n(c) * v % n(c)) / (1.*(n+1))
        uw % mean(c) = (uw % mean(c)*(1.*n) + u % n(c) * w % n(c)) / (1.*(n+1))
        vw % mean(c) = (vw % mean(c)*(1.*n) + v % n(c) * w % n(c)) / (1.*(n+1))
      end if

      !-----------------!
      !   K-eps model   !
      !-----------------!
      if(turbulence_model .eq. K_EPS) then           
        kin % mean(c) = (kin % mean(c) * (1.*n) + kin % n(c)) / (1.*(n+1))
        eps % mean(c) = (eps % mean(c) * (1.*n) + eps % n(c)) / (1.*(n+1))
      end if

      !------------------!
      !   K-eps-zeta-f   !
      !------------------!
      if(turbulence_model .eq. K_EPS_ZETA_F) then           
        kin  % mean(c) = (kin  % mean(c) * (1.*n) + kin  % n(c)) / (1.*(n+1))
        eps  % mean(c) = (eps  % mean(c) * (1.*n) + eps  % n(c)) / (1.*(n+1))
        zeta % mean(c) = (zeta % mean(c) * (1.*n) + zeta % n(c)) / (1.*(n+1))
        f22  % mean(c) = (f22  % mean(c) * (1.*n) + f22  % n(c)) / (1.*(n+1))
      end if

      !----------------------------!
      !   Reynolds stress models   !
      !----------------------------!
      if(turbulence_model .eq. HANJALIC_JAKIRLIC .or.  &
         turbulence_model .eq. REYNOLDS_STRESS) then           
        uu  % mean(c) = (uu  % mean(c) * (1.*n) + uu  % n(c)) / (1.*(n+1))
        vv  % mean(c) = (vv  % mean(c) * (1.*n) + vv  % n(c)) / (1.*(n+1))
        ww  % mean(c) = (ww  % mean(c) * (1.*n) + ww  % n(c)) / (1.*(n+1))
        uv  % mean(c) = (uv  % mean(c) * (1.*n) + uv  % n(c)) / (1.*(n+1))
        uw  % mean(c) = (uw  % mean(c) * (1.*n) + uw  % n(c)) / (1.*(n+1))
        vw  % mean(c) = (vw  % mean(c) * (1.*n) + vw  % n(c)) / (1.*(n+1))
        kin % mean(c) = (kin % mean(c) * (1.*n) + kin % n(c)) / (1.*(n+1))
        eps % mean(c) = (eps % mean(c) * (1.*n) + eps % n(c)) / (1.*(n+1))
        if(turbulence_model .eq. REYNOLDS_STRESS) then
          f22  % mean(c) = (f22  % mean(c) * (1.*n) + f22  % n(c)) / (1.*(n+1))
        end if
      end if

      !---------------------------------!
      !   Temperature and heat fluxes   !
      !---------------------------------!
      if(heat_transfer == YES) then
        t  % mean(c) = (t % mean(c) * (1.*n) + t % n(c) ) / (1.*(n+1))
        tt % mean(c) = (tt % mean(c)*(1.*n) + t % n(c) * t % n(c) ) / (1.*(n+1))
        ut % mean(c) = (ut % mean(c)*(1.*n) + u % n(c) * t % n(c) ) / (1.*(n+1))
        vt % mean(c) = (vt % mean(c)*(1.*n) + v % n(c) * t % n(c) ) / (1.*(n+1))
        wt % mean(c) = (wt % mean(c)*(1.*n) + w % n(c) * t % n(c) ) / (1.*(n+1))
      end if

      !------------------------------!
      !   User scalars are missing   !
      !------------------------------!
    end do 
  end if

  end subroutine
