!==============================================================================!
  subroutine Calculate_Sgs(grid)
!------------------------------------------------------------------------------!
!   Calculates SGS stresses and turbulent viscosity for 'LES'.                 !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Comm_Mod
  use Les_Mod
  use Rans_Mod
  use Grid_Mod
  use Grad_Mod
  use Control_Mod
  use Work_Mod, only: t_x => r_cell_01,  &
                      t_y => r_cell_02,  &
                      t_z => r_cell_03
!------------------------------------------------------------------------------!
!   Near(c) is the number of corresponding cell on the nearest wall.           !
!   In case that, in parallel executions, the subdomain does not have          !
!   any nearwall cells, the nearest_wall_cell(c) is zero.                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer :: c, s, c1, c2 
  real    :: nx, ny, nz
  real    :: cs
  real    :: lf, u_tau_l, u_f 
  real    :: u_tot, u_nor, u_tan, a_pow, b_pow, nu, dely
  real    :: nc2
!==============================================================================!
  
  !---------------!
  !               !
  !   SGS terms   !
  !               !
  !---------------!
  if(buoyancy) then
    call Grad_Mod_For_Phi(grid, t % n, 1, t_x, .true.)  ! dT/dx
    call Grad_Mod_For_Phi(grid, t % n, 2, t_y, .true.)  ! dT/dy
    call Grad_Mod_For_Phi(grid, t % n, 3, t_z, .true.)  ! dT/dz
  end if 
 
  if(turbulence_model .eq. LES_SMAGORINSKY) then
    do c = 1, grid % n_cells
      lf = grid % vol(c)**ONE_THIRD

      ! if(nearest_wall_cell(c) .ne. 0) is needed for parallel version
      ! since the subdomains which do not "touch" wall
      ! has nearest_wall_cell(c) = 0. 
      if(nearest_wall_cell(c) .ne. 0) then
        u_f = sqrt( viscosity  &
                    * sqrt(  u % n(nearest_wall_cell(c)) ** 2   & 
                           + v % n(nearest_wall_cell(c)) ** 2   &
                           + w % n(nearest_wall_cell(c)) ** 2)  &
                   / (grid % wall_dist(nearest_wall_cell(c))+TINY) )
        y_plus(c) = grid % wall_dist(c) * u_f / viscosity
        cs = c_smag * (1.0 - exp(-y_plus(c) / 25.0))
      else  
        cs = c_smag
      end if
      vis_t(c) = density &
              * (lf*lf)           &          ! delta^2 
              * (cs*cs)           &          ! cs^2   
              * shear(c) 
    end do

  else if(turbulence_model .eq. LES_DYNAMIC) then
    if(buoyancy) then  
      do c = 1, grid % n_cells
        lf = grid % vol(c)**ONE_THIRD  
        vis_t(c) = density            &
                * (lf*lf)             &          ! delta^2 
                * c_dyn(c)            &          ! c_dynamic   
                * shear(c) 
      end do
    else
      ! lf is not initialized here
      do c = 1, grid % n_cells
        vis_t(c) = density                &
                * (lf*lf)                 &          ! delta^2 
                * c_dyn(c)                &          ! c_dynamic   
                * sqrt(shear(c)*shear(c)  &
                + 2.5*(grav_x*t_x(c) + grav_y*t_y(c) + grav_z*t_z(c)))  
      end do
    end if

  else if(turbulence_model .eq. LES_WALE) then
    do c = 1, grid % n_cells
      lf = grid % vol(c)**ONE_THIRD
      vis_t(c) = density           &
              * (lf*lf)            &          ! delta^2 
              * (0.5*0.5)          &          ! cs^2   
              * wale_v(c) 
    end do
  end if

  if(buoyancy) then
    do c = 1, grid % n_cells
      nc2 = -(  grav_x * t_x(c)   &
              + grav_y * t_y(c)   &
              + grav_z * t_z(c))  &
          / t_ref
      nc2 = max(0.0, nc2) 
      vis_t(c) = vis_t(c) * sqrt(1.0 - min(2.5*nc2/(shear(c)**2), 1.0))
    end do
  end if

  !-------------------!
  !   Wall function   !
  !-------------------+--------------!
  !   Law of the wall:               !
  !                                  !
  !   u+ = yz+  for z+ < 11.81       !
  !                                  !
  !   and                            !
  !                                  !
  !   u+ = A(y+)^B   for y+ > 11.81  !
  !                                  !
  !   with: A = 8.3 and B = 1/7      !
  !                                  !
  !----------------------------------+----------!
  !   The procedure below should be activated   !
  !   only if wall function approach is used.   !
  !----------------.----------------------------! 
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2  < 0) then 

      nx = grid % sx(s) / grid % s(s)
      ny = grid % sy(s) / grid % s(s)
      nz = grid % sz(s) / grid % s(s)

      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then

        u_tot = sqrt(  u % n(c1) * u % n(c1)     &
                     + v % n(c1) * v % n(c1)     & 
                     + w % n(c1) * w % n(c1)  )

        u_nor = (   u % n(c1) * nx  &
                  + v % n(c1) * ny  &
                  + w % n(c1) * nz )   

        if( abs(u_tot) > abs(u_nor) ) then
          u_tan = sqrt(u_tot * u_tot - u_nor * u_nor)
        else
          u_tan = TINY 
        end if

        a_pow = 8.3
        b_pow = 1.0/7.0
        nu = viscosity/density
        dely = grid % wall_dist(c1)

        ! Calculate u_tau_l
        u_tau_l = ( u_tan/a_pow * (nu/dely)**b_pow )                     &
          ** (1.0/(1.0+b_pow))

        ! Calculate tau_wall 
        tau_wall(c1) = viscosity * u_tan / dely 
 
        ! Calculate y+
        y_plus(c1)  = dely*u_tau_l/nu
        if(y_plus(c1)  >=  11.81) then
          ! This one is effective viscosity
          vis_wall(c1) = density*u_tau_l*u_tau_l*dely/abs(u_tan) 
        else 
          vis_wall(c1) = viscosity + fw(s)*vis_t(c1)+(1.0-fw(s))*vis_t(c2)
        end if
      end if  ! Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL or WALLFL
    end if    ! c2 < 0
  end do

  call Comm_Mod_Exchange_Real(grid, vis_t)
  call Comm_Mod_Exchange_Real(grid, vis_wall)

  end subroutine
