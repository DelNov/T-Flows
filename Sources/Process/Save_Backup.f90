!==============================================================================!
  subroutine Save_Backup(grid, time_step, name_save)
!------------------------------------------------------------------------------!
!   Writes backup files. name.backup                                           !
!----------------------------------[Modules]-----------------------------------!
  use Name_Mod, only: problem_name
  use Const_Mod
  use Flow_Mod
  use Comm_Mod
  use Rans_Mod
  use Grid_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  integer          :: time_step
  character(len=*) :: name_save
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: name_out, store_name
  integer           :: fh, d
!==============================================================================!

  store_name = problem_name

  problem_name = name_save

  ! Name backup file 
  call Name_File(0, name_out, '.backup')

  ! Open backup file
  call Comm_Mod_Open_File_Write(fh, name_out)

  ! Create new types
  call Comm_Mod_Create_New_Types()

  ! Initialize displacement
  d = 0

  !-----------------------------------------------------------------------!
  !   Save cell-centre coordinates.  Could be useful for interpolations   !
  !-----------------------------------------------------------------------! 
  call Write_Backup_Cell_Bnd(fh, d, 'x_coordinates', grid % xc(-nb_s:nc_s))
  call Write_Backup_Cell_Bnd(fh, d, 'y_coordinates', grid % yc(-nb_s:nc_s))
  call Write_Backup_Cell_Bnd(fh, d, 'z_coordinates', grid % zc(-nb_s:nc_s))

  !---------------!
  !               !
  !   Save data   !
  !               !
  !---------------! 

  ! Time step
  call Write_Backup_Int(fh, d, 'time_step', time_step)

  ! Number of processors
  call Write_Backup_Int(fh, d, 'n_proc', n_proc)

  ! Bulk flows and pressure drops in each direction
  call Write_Backup_Real(fh, d, 'bulk_flux_x',   bulk(1) % flux_x)
  call Write_Backup_Real(fh, d, 'bulk_flux_y',   bulk(1) % flux_y)
  call Write_Backup_Real(fh, d, 'bulk_flux_z',   bulk(1) % flux_z)
  call Write_Backup_Real(fh, d, 'bulk_p_drop_x', bulk(1) % p_drop_x)
  call Write_Backup_Real(fh, d, 'bulk_p_drop_y', bulk(1) % p_drop_y)
  call Write_Backup_Real(fh, d, 'bulk_p_drop_z', bulk(1) % p_drop_z)

  !--------------!
  !   Velocity   !
  !--------------!
  call Write_Backup_Variable(fh, d, 'u_velocity', u)
  call Write_Backup_Variable(fh, d, 'v_velocity', v)
  call Write_Backup_Variable(fh, d, 'w_velocity', w)

  !--------------------------------------!
  !   Pressure and pressure correction   !
  !--------------------------------------!
  call Write_Backup_Cell_Bnd(fh, d, 'press',       p  % n(-nb_s:nc_s))
  call Write_Backup_Cell_Bnd(fh, d, 'press_corr',  pp % n(-nb_s:nc_s))

  !----------------------------------------------------!
  !   Mass flow rates (ask Egor if name is correct?)   !
  !----------------------------------------------------!
  call Write_Backup_Face(fh, d, 'mass_flow_rate', flux(1:nf_s+nbf_s))

  !--------------!
  !   Etnhalpy   !
  !--------------!
  if(heat_transfer .eq. YES) then
    call Write_Backup_Variable(fh, d, 'temp', t)
  end if

  !-----------------------!
  !   Turbulence models   !
  !-----------------------!
  if(turbulence_model .eq. K_EPS    .or.  &
     turbulence_model .eq. K_EPS_ZETA_F     .or.  &
     turbulence_model .eq. HYBRID_K_EPS_ZETA_F) then

    ! K and epsilon
    call Write_Backup_Variable(fh, d, 'kin', kin)
    call Write_Backup_Variable(fh, d, 'eps', eps)

    ! Other turbulent quantities
    call Write_Backup_Cell_Bnd(fh, d, 'p_kin',    p_kin   (-nb_s:nc_s))
    call Write_Backup_Cell_Bnd(fh, d, 'u_tau',    u_tau   (-nb_s:nc_s))
    call Write_Backup_Cell_Bnd(fh, d, 'y_plus',   y_plus  (-nb_s:nc_s))
    call Write_Backup_Cell_Bnd(fh, d, 'vis_wall', vis_wall(-nb_s:nc_s))
    call Write_Backup_Cell    (fh, d, 'tau_wall', tau_wall  (1:nc_s))
  end if

  if(turbulence_model .eq. K_EPS_ZETA_F     .or.  &
     turbulence_model .eq. HYBRID_K_EPS_ZETA_F) then

    ! Zeta and f22
    call Write_Backup_Variable(fh, d, 'zeta', zeta)
    call Write_Backup_Variable(fh, d, 'f22',  f22)

    ! Other turbulent quantities
    call Write_Backup_Cell_Bnd(fh, d, 't_scale',  t_scale(-nb_s:nc_s))
    call Write_Backup_Cell_Bnd(fh, d, 'l_scale',  l_scale(-nb_s:nc_s))
  end if 

  if(turbulence_model .eq. REYNOLDS_STRESS .or.  &                          
     turbulence_model .eq. HANJALIC_JAKIRLIC) then                          

    ! Reynolds stresses
    call Write_Backup_Variable(fh, d, 'uu',  uu)
    call Write_Backup_Variable(fh, d, 'vv',  vv)
    call Write_Backup_Variable(fh, d, 'ww',  ww)
    call Write_Backup_Variable(fh, d, 'uv',  uv)
    call Write_Backup_Variable(fh, d, 'uw',  uw)
    call Write_Backup_Variable(fh, d, 'vw',  vw)

    ! Epsilon
    call Write_Backup_Variable(fh, d, 'eps', eps)

    ! F22
    if(turbulence_model .eq. HANJALIC_JAKIRLIC) then
      call Write_Backup_Variable(fh, d, 'f22',  f22)
    end if

    ! Other turbulent quantities ?
  end if 

  ! Close backup file
  call Comm_Mod_Close_File(fh)

  problem_name = store_name

  end subroutine

!TEST  ! Test face buffers
!TEST  do s = 1, grid % n_faces
!TEST    flux(s) = 1000000.0 * grid % xf(s) +  &
!TEST                 1000.0 * grid % yf(s) +  &
!TEST                    1.0 * grid % zf(s)
!TEST    print '(a6,i4.4,a4,f18.3)', ' flux(', s, ') = ', flux(s)
!TEST  end do
