!==============================================================================!
  subroutine Backup_Mod_Save(Fld, Swr, Tur, Vof, time, time_step, domain)
!------------------------------------------------------------------------------!
!   Saves backup files name.backup                                             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),  target :: Fld
  type(Swarm_Type),  target :: Swr
  type(Turb_Type),   target :: Tur
  type(Vof_Type),    target :: Vof
  real                      :: time            ! time of simulation
  integer                   :: time_step       ! current time step
  integer,         optional :: domain
!-----------------------------------[Locals]-----------------------------------!
  type(Comm_Type), pointer :: Comm
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: phi
  character(SL)            :: name_out, name_mean, a_name
  integer                  :: vc, sc, ua
  integer(DP)              :: d
!==============================================================================!

  call Profiler % Start('Backup_Mod_Save')

  ! Take aliases
  Grid => Fld % pnt_grid
  bulk => Fld % bulk
  Comm => Grid % Comm

  ! Name backup file
  call File % Set_Name(name_out, time_step=time_step,  &
                       extension='.backup', domain=domain)

  ! Open backup file
  call Comm % Open_File_Write(fh, name_out)

  ! Create new types
  call Comm % Create_New_Types()

  ! Initialize displacement
  d = 0

  ! Intialize number of stored variables
  vc = 0

  !-----------------------------------------------------------------------!
  !   Save cell-centre coordinates.  Could be useful for interpolations   !
  !-----------------------------------------------------------------------!
  call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'x_coords', Grid % xc)
  call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'y_coords', Grid % yc)
  call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'z_coords', Grid % zc)

  !---------------!
  !               !
  !   Save data   !
  !               !
  !---------------!

  ! Time step
  call Backup_Mod_Write_Int(Comm, d, vc, 'time_step', time_step)

  ! Simulation time
  call Backup_Mod_Write_Real(Comm, d, vc, 'time', time)

  ! Number of processors
  call Backup_Mod_Write_Int(Comm, d, vc, 'n_proc', n_proc)

  ! Bulk flows and pressure drops in each direction
  call Backup_Mod_Write_Real(Comm, d, vc, 'bulk_flux_x',   bulk % flux_x)
  call Backup_Mod_Write_Real(Comm, d, vc, 'bulk_flux_y',   bulk % flux_y)
  call Backup_Mod_Write_Real(Comm, d, vc, 'bulk_flux_z',   bulk % flux_z)
  call Backup_Mod_Write_Real(Comm, d, vc, 'bulk_flux_x_o', bulk % flux_x_o)
  call Backup_Mod_Write_Real(Comm, d, vc, 'bulk_flux_y_o', bulk % flux_y_o)
  call Backup_Mod_Write_Real(Comm, d, vc, 'bulk_flux_z_o', bulk % flux_z_o)
  call Backup_Mod_Write_Real(Comm, d, vc, 'bulk_p_drop_x', bulk % p_drop_x)
  call Backup_Mod_Write_Real(Comm, d, vc, 'bulk_p_drop_y', bulk % p_drop_y)
  call Backup_Mod_Write_Real(Comm, d, vc, 'bulk_p_drop_z', bulk % p_drop_z)

  !----------------------------!
  !                            !
  !   Navier-Stokes equation   !
  !                            !
  !----------------------------!

  !--------------!
  !   Velocity   !
  !--------------!
  call Backup_Mod_Write_Variable(d, vc, 'u_velocity', Fld % u)
  call Backup_Mod_Write_Variable(d, vc, 'v_velocity', Fld % v)
  call Backup_Mod_Write_Variable(d, vc, 'w_velocity', Fld % w)

  !------------------------------------------------------!
  !   Pressure, its gradients, and pressure correction   !
  !------------------------------------------------------!
  call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'press',      Fld %  p % n)
  call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'press_x',    Fld %  p % x)
  call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'press_y',    Fld %  p % y)
  call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'press_z',    Fld %  p % z)
  call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'press_corr', Fld % pp % n)

  !-------------------!
  !   Volume fluxes   ! -> don't use for the time being, too much trouble
  !-------------------!

  !-------------------------!
  !   Does it have outlet   !
  !-------------------------!
  ! Update on July 17, 2021: I have some reservations about this part, since
  ! there was another bug fix when computing fluxes in the meanwhile (check:
  ! 90f77a1c8bd4ca05330a4435ed6321782ef00199).  This balancing also caused a
  ! bug when loading backup file (also check "Compute_Pressure" as well as
  ! "Initialize_Variables" procedures)
  !
  ! Update on February 27, 2022: I have also added "has_outflow_boundary"
  ! to be able to tell PETSc if matrix for pressure is singular
  !
  ! Update on June 2, 2022: Unified all outlet boundaries into one
  ! to be able to tell PETSc if matrix for pressure is singular
  call Backup_Mod_Write_Log(Comm, d, vc, 'has_pressure', Fld % has_pressure)

  !--------------!
  !              !
  !   Enthalpy   !
  !              !
  !--------------!
  if(Fld % heat_transfer) then
    call Backup_Mod_Write_Variable(d, vc, 'temp', Fld % t)
  end if

  !--------------!
  !              !
  !  Multiphase  !
  !              !
  !--------------!
  if(Fld % with_interface) then
    call Backup_Mod_Write_Variable(d, vc, 'vof_fun', Vof % fun)
  end if

  !-----------------------!
  !                       !
  !   Turbulence models   !
  !                       !
  !-----------------------!

  !-----------------!
  !   K-eps model   !
  !-----------------!
  if(Tur % model .eq. K_EPS) then

    ! K and epsilon
    call Backup_Mod_Write_Variable(d, vc, 'kin', Tur % kin)
    call Backup_Mod_Write_Variable(d, vc, 'eps', Tur % eps)

    ! Other turbulent quantities
    call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'p_kin',  Tur % p_kin )
    call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'y_plus', Tur % y_plus)
    call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'vis_t',  Tur % vis_t )
    call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'vis_w',  Tur % vis_w )

    ! Turbulence quantities connected with heat transfer
    if(Fld % heat_transfer) then
      call Backup_Mod_Write_Variable(d, vc, 't2',    Tur % t2)
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'p_t2',  Tur % p_t2 )
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'con_w', Tur % con_w)
    end if

  end if

  !------------------------!
  !   K-eps-zeta-f model   !
  !------------------------!
  if(Tur % model .eq. K_EPS_ZETA_F .or.  &
     Tur % model .eq. HYBRID_LES_RANS) then

    ! K, eps, zeta and f22
    call Backup_Mod_Write_Variable(d, vc, 'kin',  Tur % kin)
    call Backup_Mod_Write_Variable(d, vc, 'eps',  Tur % eps)
    call Backup_Mod_Write_Variable(d, vc, 'zeta', Tur % zeta)
    call Backup_Mod_Write_Variable(d, vc, 'f22',  Tur % f22)

    ! Other turbulent quantities
    call Backup_Mod_Write_Cell_Real(Grid, d, vc,'p_kin',   Tur % p_kin  )
    call Backup_Mod_Write_Cell_Real(Grid, d, vc,'y_plus',  Tur % y_plus )
    call Backup_Mod_Write_Cell_Real(Grid, d, vc,'vis_t',   Tur % vis_t  )
    call Backup_Mod_Write_Cell_Real(Grid, d, vc,'vis_w',   Tur % vis_w  )
    call Backup_Mod_Write_Cell_Real(Grid, d, vc,'t_scale', Tur % t_scale)
    call Backup_Mod_Write_Cell_Real(Grid, d, vc,'l_scale', Tur % l_scale)

    if(Fld % heat_transfer) then
      call Backup_Mod_Write_Variable(d, vc, 't2',    Tur % t2)
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'p_t2',  Tur % p_t2 )
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'con_w', Tur % con_w)
    end if

  end if

  !----------------------------!
  !   Reynolds stress models   !
  !----------------------------!
  if(Tur % model .eq. RSM_MANCEAU_HANJALIC .or.  &
     Tur % model .eq. RSM_HANJALIC_JAKIRLIC) then

    ! Reynolds stresses
    call Backup_Mod_Write_Variable(d, vc, 'uu',  Tur % uu)
    call Backup_Mod_Write_Variable(d, vc, 'vv',  Tur % vv)
    call Backup_Mod_Write_Variable(d, vc, 'ww',  Tur % ww)
    call Backup_Mod_Write_Variable(d, vc, 'uv',  Tur % uv)
    call Backup_Mod_Write_Variable(d, vc, 'uw',  Tur % uw)
    call Backup_Mod_Write_Variable(d, vc, 'vw',  Tur % vw)

    ! Epsilon
    call Backup_Mod_Write_Variable(d, vc, 'eps', Tur % eps)

    ! F22
    if(Tur % model .eq. RSM_MANCEAU_HANJALIC) then
      call Backup_Mod_Write_Variable(d, vc, 'f22',  Tur % f22)
    end if

    ! Other turbulent quantities
    call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'vis_t', Tur % vis_t)
    call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'vis_w', Tur % vis_w)

    ! Turbulence quantities connected with heat transfer
    if(Fld % heat_transfer) then
      call Backup_Mod_Write_Cell_Real(Grid, d, vc,'con_w', Tur % con_w)
    end if
  end if

  !--------------!
  !   Roughness  !
  !--------------!
  if(Tur % rough_walls) then
    call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'z_o', Tur % z_o)
  end if

  !------------------!
  !   Save scalars   !
  !------------------!
  do sc = 1, Fld % n_scalars
    phi => Fld % scalar(sc)
    call Backup_Mod_Write_Variable(d, vc, phi % name, phi)
  end do

  !-----------------------------------------!
  !                                         !
  !   Turbulent statistics for all models   !
  !                                         !
  !-----------------------------------------!
  if(Tur % statistics) then

    call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'u_mean', Tur % u_mean)
    call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'v_mean', Tur % v_mean)
    call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'w_mean', Tur % w_mean)
    call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'p_mean', Tur % p_mean)
    if(Fld % heat_transfer) then
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 't_mean', Tur % t_mean)
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'q_mean', Tur % q_mean)
    end if

    ! K and epsilon
    if(Tur % model .eq. K_EPS) then
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'kin_mean',  &
                                      Tur % kin_mean)
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'eps_mean',  &
                                      Tur % eps_mean)
      if(Fld % heat_transfer) then
        call Backup_Mod_Write_Cell_Real(Grid, d, vc, 't2_mean',  &
                                        Tur % t2_mean)
        call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'ut_mean',  &
                                        Tur % ut_mean)
        call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'vt_mean',  &
                                        Tur % vt_mean)
        call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'wt_mean',  &
                                        Tur % wt_mean)
      end if
    end if

    ! K-eps-zeta-f and the hybrid model
    if(Tur % model .eq. K_EPS_ZETA_F .or.  &
       Tur % model .eq. HYBRID_LES_RANS) then
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'kin_mean',  &
                                      Tur % kin_mean )
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'eps_mean',  &
                                      Tur % eps_mean )
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'zeta_mean',  &
                                      Tur % zeta_mean)
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'f22_mean',  &
                                      Tur % f22_mean )
      if(Fld % heat_transfer) then
        call Backup_Mod_Write_Cell_Real(Grid, d, vc, 't2_mean',  &
                                        Tur % t2_mean)
        call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'ut_mean',  &
                                        Tur % ut_mean)
        call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'vt_mean',  &
                                        Tur % vt_mean)
        call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'wt_mean',  &
                                        Tur % wt_mean)
      end if
    end if

    ! Reynolds stress models
    if(Tur % model .eq. RSM_MANCEAU_HANJALIC .or.  &
       Tur % model .eq. RSM_HANJALIC_JAKIRLIC) then
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'uu_mean',  &
                                      Tur % uu_mean)
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'vv_mean',  &
                                      Tur % vv_mean)
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'ww_mean',  &
                                      Tur % ww_mean)
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'uv_mean',  &
                                      Tur % uv_mean)
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'uw_mean',  &
                                      Tur % uw_mean)
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'vw_mean',  &
                                      Tur % vw_mean)
    end if

    call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'uu_res', Tur % uu_res)
    call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'vv_res', Tur % vv_res)
    call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'ww_res', Tur % ww_res)
    call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'uv_res', Tur % uv_res)
    call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'uw_res', Tur % uw_res)
    call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'vw_res', Tur % vw_res)

    if(Fld % heat_transfer) then
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 't2_res', Tur % t2_res)
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'ut_res', Tur % ut_res)
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'vt_res', Tur % vt_res)
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'wt_res', Tur % wt_res)
    end if

    ! Scalars
    do sc = 1, Fld % n_scalars
      phi => Fld % scalar(sc)
      name_mean = phi % name
      name_mean(5:9) = '_mean'
      call Backup_Mod_Write_Cell_Real(Grid, d, vc, name_mean,  &
                                      Tur % scalar_mean(sc, :))
    end do

  end if

  !--------------------------!
  !                          !
  !   Swarm (of particles)   !
  !                          !
  !--------------------------!
  if(Fld % with_particles) then
    call Backup_Mod_Write_Swarm(d, vc, Swr)
    call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'n_deposited',      &
                                    Swr % n_deposited)
    call Backup_Mod_Write_Cell_Real(Grid, d, vc, 'n_reflected',      &
                                    Swr % n_reflected)
  end if

  !-----------------!
  !                 !
  !   User arrays   !
  !                 !
  !-----------------!

  do ua = 1, Grid % n_user_arrays
    a_name = 'A_??'
    write(a_name(3:4),'(I2.2)') ua
    call Backup_Mod_Write_Cell_Real(Grid, d, vc, a_name,  &
                                    Grid % user_array(ua, :))
  end do

  ! Variable count (store +1 to count its own self)
  call Backup_Mod_Write_Int(Comm, d, vc, 'variable_count', vc + 1)

  if(this_proc < 2) then
    print *, '# Wrote ', vc, ' variables!'
  end if

  ! Close backup file
  call Comm % Close_File(fh)

  call Profiler % Stop('Backup_Mod_Save')

  end subroutine
