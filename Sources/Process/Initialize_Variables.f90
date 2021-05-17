!==============================================================================!
  subroutine Initialize_Variables(Flow, turb, Vof, swarm, Sol)
!------------------------------------------------------------------------------!
!   Initialize dependent variables.  (It is a bit of a mess still)             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use File_Mod
  use Field_Mod,   only: Field_Type
  use Comm_Mod
  use Turb_Mod
  use Swarm_Mod
  use Grid_Mod
  use Bulk_Mod
  use User_Mod
  use Control_Mod
  use Vof_Mod
  use Numerics_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: Flow
  type(Turb_Type),  target :: turb
  type(Vof_Type)           :: Vof
  type(Swarm_Type)         :: swarm
  type(Solver_Type)        :: Sol
!----------------------------------[Calling]-----------------------------------!
  integer :: Key_Ind
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Bulk_Type),  pointer :: bulk
  type(Var_Type),   pointer :: u, v, w, t, phi
  type(Var_Type),   pointer :: kin, eps, f22, zeta, vis, t2
  type(Var_Type),   pointer :: uu, vv, ww, uv, uw, vw
  type(Face_Type),  pointer :: v_flux
  real, contiguous, pointer :: u_mean(:), v_mean(:), w_mean(:)
  integer                   :: i, c, c1, c2, m, s, nks, nvs, sc, fu
  integer                   :: n_wall, n_inflow, n_outflow, n_symmetry,  &
                               n_heated_wall, n_pressure, n_convect
  character(SL)             :: keys(128)
  character(SL)             :: keys_file(128)
  real                      :: vals(0:128)   ! note that they start from zero!
  real                      :: area

  integer                   :: n_points, k
  real, allocatable         :: prof(:,:), x(:), y(:), z(:), dist(:)
  logical                   :: found

  ! Default values for initial conditions
  real, parameter           :: u_def   = 0.0,  v_def    = 0.0,  w_def    = 0.0
  real, parameter           :: t_def   = 0.0,  t2_def   = 0.0,  phi_def  = 0.0
  real, parameter           :: kin_def = 0.0,  eps_def  = 0.0,  f22_def  = 0.0
  real, parameter           :: vis_def = 0.0,  zeta_def = 0.0
  real, parameter           :: uu_def  = 0.0,  vv_def   = 0.0,  ww_def   = 0.0
  real, parameter           :: uv_def  = 0.0,  uw_def   = 0.0,  vw_def   = 0.0
!==============================================================================!

  ! Take aliases
  grid     => Flow % pnt_grid
  bulk     => Flow % bulk
  v_flux   => Flow % v_flux
  vis      => turb % vis
  u_mean   => turb % u_mean
  v_mean   => turb % v_mean
  w_mean   => turb % w_mean
  call Flow % Alias_Momentum(u, v, w)
  call Flow % Alias_Energy  (t)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)
  call Turb_Mod_Alias_Stresses    (turb, uu, vv, ww, uv, uw, vw)
  call Turb_Mod_Alias_T2          (turb, t2)

  area  = 0.0
  if (this_proc < 2) print *, '# Grid name: ', grid % name

  ! Found the line where boundary condition definition is defined
  call Control_Mod_Position_At_One_Key('INITIAL_CONDITION', &
                                       found,               &
                                       .true.)

  !-----------------------------------------------!
  !                                               !
  !   Found the section with initial conditions   !
  !                                               !
  !-----------------------------------------------!
  if (found) then

    call Control_Mod_Read_Strings_On('VARIABLES', keys, nks, .true.)

    ! Input is valid, turn keys to upper case
    do i = 1, nks
      call To_Upper_Case(keys(i))
    end do

    ! Check if there is file specified
    call Control_Mod_Read_Strings_On('FILE', keys_file, nvs, .true.)

    !-----------------------------------------------!
    !   Initial conditions is prescribed in a file  !
    !-----------------------------------------------!
    if (nvs .eq. 1) then ! word 'file' was specified

      if (this_proc < 2) &
        print *, '# Values specified in the file: ', trim(keys_file(nvs))

      call File_Mod_Open_File_For_Reading(keys_file(1), fu, this_proc)

      ! number of points
      call File_Mod_Read_Line(fu)

      read(line % tokens(1), *) n_points

      if (this_proc < 2) print '(a,i0,2a)', " # Reading ", nks, &
        " columns in file " , trim(keys_file(1))

      allocate(prof(n_points, 0:nks)); prof = 0.
      allocate(x(n_points));           x    = 0.
      allocate(y(n_points));           y    = 0.
      allocate(z(n_points));           z    = 0.
      allocate(dist(n_points));        dist = 0.

      ! Read the entire profile file
      do m = 1, n_points
        call File_Mod_Read_Line(fu)
        do i = 1, nks
          read(line % tokens(i), *) prof(m, i)
        end do
      end do
      close(fu)

      ! A plane is defined
      if (keys(1) .eq. 'X' .and. keys(2) .eq. 'Y' .or.  &
          keys(1) .eq. 'X' .and. keys(2) .eq. 'Z' .or.  &
          keys(1) .eq. 'Y' .and. keys(2) .eq. 'Z') then

        ! Set the closest point
        do c = 1, grid % n_cells

          i=Key_Ind('X', keys, nks); x(:) = prof(:,i)
          i=Key_Ind('Y', keys, nks); y(:) = prof(:,i)
          i=Key_Ind('Z', keys, nks); z(:) = prof(:,i)

          ! do no waste time on sqrt((r-r0)^2) -> use (r-r0)^2
          if(keys(1) .eq. 'Y' .and. keys(2) .eq. 'Z') then
            dist(:) = (y(:)-grid % yc(c))**2 + (z(:)-grid % zc(c))**2
          else if(keys(1) .eq. 'X' .and. keys(2) .eq. 'Z') then
            dist(:) = (x(:)-grid % xc(c))**2 + (z(:)-grid % zc(c))**2
          else if(keys(1) .eq. 'X' .and. keys(2) .eq. 'Y') then
            dist(:) = (x(:)-grid % xc(c))**2 + (y(:)-grid % yc(c))**2
          end if

          ! Store closest point in k
          k = minloc(dist, dim = 1)

          i=Key_Ind('U',keys,nks);prof(k,0)=u_def;u%n(c)=prof(k,i)
          i=Key_Ind('V',keys,nks);prof(k,0)=v_def;v%n(c)=prof(k,i)
          i=Key_Ind('W',keys,nks);prof(k,0)=w_def;w%n(c)=prof(k,i)

          if(Flow % heat_transfer) then
            i=Key_Ind('T',keys,nks);prof(k,0)=t_def;t%n(c)=prof(k,i)
          end if

          if(turb % model .eq. K_EPS) then
            i=Key_Ind('KIN',keys,nks);prof(k,0)=kin_def; kin%n(c)=prof(k,i)
            i=Key_Ind('EPS',keys,nks);prof(k,0)=eps_def; eps%n(c)=prof(k,i)
            if(Flow % heat_transfer) then
              i=Key_Ind('T2', keys,nks);prof(k,0)=t2_def; t2 %n(c)=prof(k,i)
            end if
          end if

          if(turb % model .eq. K_EPS_ZETA_F .or.  &
             turb % model .eq. HYBRID_LES_RANS) then
            i=Key_Ind('KIN', keys,nks);prof(k,0)=kin_def; kin %n(c)=prof(k,i)
            i=Key_Ind('EPS', keys,nks);prof(k,0)=eps_def; eps %n(c)=prof(k,i)
            i=Key_Ind('ZETA',keys,nks);prof(k,0)=zeta_def;zeta%n(c)=prof(k,i)
            i=Key_Ind('F22', keys,nks);prof(k,0)=f22_def; f22 %n(c)=prof(k,i)
            if(Flow % heat_transfer) then
              i=Key_Ind('T2', keys,nks);prof(k,0)=t2_def; t2 %n(c)=prof(k,i)
            end if
           end if

          if(turb % model .eq. DES_SPALART) then
            i=Key_Ind('VIS',keys,nks); prof(k,0)=vis_def; vis%n(c)=prof(k,i)
          end if

          if(turb % model .eq. RSM_MANCEAU_HANJALIC .or. &
             turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
            i=Key_Ind('UU', keys,nks);prof(k,0)=uu_def; uu %n(c)=prof(k,i)
            i=Key_Ind('VV', keys,nks);prof(k,0)=vv_def; vv %n(c)=prof(k,i)
            i=Key_Ind('WW', keys,nks);prof(k,0)=ww_def; ww %n(c)=prof(k,i)
            i=Key_Ind('UV', keys,nks);prof(k,0)=uv_def; uv %n(c)=prof(k,i)
            i=Key_Ind('UW', keys,nks);prof(k,0)=uw_def; uw %n(c)=prof(k,i)
            i=Key_Ind('VW', keys,nks);prof(k,0)=vw_def; vw %n(c)=prof(k,i)
            i=Key_Ind('EPS',keys,nks);prof(k,0)=eps_def;eps%n(c)=prof(k,i)
            if (turb % model .eq. RSM_MANCEAU_HANJALIC) then
              i=Key_Ind('F22',keys,nks);prof(k,0)=f22_def;f22%n(c)=prof(k,i)
            end if
          end if

        end do ! c = 1, grid % n_cells

        call Comm_Mod_Wait
        deallocate(prof)
        deallocate(x)
        deallocate(y)
        deallocate(z)
        deallocate(dist)

      end if

    !---------------------------------------------------!
    !   Initial conditions is not prescribed in a file  !
    !---------------------------------------------------!
    else

      ! Go back to key and read again
      call Control_Mod_Position_At_One_Key('INITIAL_CONDITION', &
                                           found,               &
                                           .true.)

      call Control_Mod_Read_Strings_On('VARIABLES', keys, nks, .true.)

      ! Input is valid, turn keys to upper case
      do i = 1, nks
        call To_Upper_Case(keys(i))
      end do

      call Control_Mod_Read_Real_Array_On('VALUES', vals(1), nvs, .true.)

      ! Check validity of the input
      if(nks .eq. 0 .or. nvs .eq. 0 .and. this_proc < 2) then
        print '(2a)', '# Critical, for initial condition: ',        &
                      ' no values or variables have been provided'
        call Comm_Mod_End
        stop
      end if
      if(nks .ne. nvs .and. this_proc < 2) then
        print '(2a)', '# Critical for initial conditions, number of values ',  &
                      ' is not the same as number of provided variable names'
        call Comm_Mod_End
        stop
      end if

      ! Input is valid, turn keys to upper case
      do i = 1, nks
        call To_Upper_Case(keys(i))
      end do

      do c = 1, grid % n_cells

        if(turb % statistics) then
          u_mean(c) = 0.0
          v_mean(c) = 0.0
          w_mean(c) = 0.0
        end if

        vals(0) = u_def;  u % n(c) = vals(Key_Ind('U', keys, nks))
        vals(0) = v_def;  v % n(c) = vals(Key_Ind('V', keys, nks))
        vals(0) = w_def;  w % n(c) = vals(Key_Ind('W', keys, nks))

        u % o(c)  = u % n(c)
        u % oo(c) = u % n(c)
        v % o(c)  = v % n(c)
        v % oo(c) = v % n(c)
        w % o(c)  = w % n(c)
        w % oo(c) = w % n(c)

        if(Flow % heat_transfer) then
          vals(0) = t_def;  t % n(c) = vals(Key_Ind('T', keys, nks))
          t % o(c)  = t % n(c)
          t % oo(c) = t % n(c)
        end if

        ! Scalars
        do sc = 1, Flow % n_scalars
          phi => Flow % scalar(sc)
          vals(0) = phi_def
          phi % n(c)  = vals(Key_Ind(phi % name, keys, nks))
          phi % o(c)  = phi % n(c)
          phi % oo(c) = phi % n(c)
        end do

        if(turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
           turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
          vals(0) = uu_def;  uu  % n(c) = vals(Key_Ind('UU',  keys, nks))
          vals(0) = vv_def;  vv  % n(c) = vals(Key_Ind('VV',  keys, nks))
          vals(0) = ww_def;  ww  % n(c) = vals(Key_Ind('WW',  keys, nks))
          vals(0) = uv_def;  uv  % n(c) = vals(Key_Ind('UV',  keys, nks))
          vals(0) = uw_def;  uw  % n(c) = vals(Key_Ind('UW',  keys, nks))
          vals(0) = vw_def;  vw  % n(c) = vals(Key_Ind('VW',  keys, nks))
          vals(0) = eps_def; eps % n(c) = vals(Key_Ind('EPS', keys, nks))
          uu % o(c)  = uu % n(c)
          uu % oo(c) = uu % n(c)
          vv % o(c)  = vv % n(c)
          vv % oo(c) = vv % n(c)
          ww % o(c)  = ww % n(c)
          ww % oo(c) = ww % n(c)
          uv % o(c)  = uv % n(c)
          uv % oo(c) = uv % n(c)
          uw % o(c)  = uw % n(c)
          uw % oo(c) = uw % n(c)
          vw % o(c)  = vw % n(c)
          vw % oo(c) = vw % n(c)
          if(turb % model .eq. RSM_MANCEAU_HANJALIC) then
            vals(0) = f22_def; f22 % n(c) = vals(Key_Ind('F22', keys, nks))
            f22 % o(c)  = f22 % n(c)
            f22 % oo(c) = f22 % n(c)
          end if
        end if

        if(turb % model .eq. K_EPS) then
          vals(0) = kin_def; kin % n(c) = vals(Key_Ind('KIN', keys, nks))
          vals(0) = eps_def; eps % n(c) = vals(Key_Ind('EPS', keys, nks))
          kin % o(c)  = kin % n(c)
          kin % oo(c) = kin % n(c)
          eps % o(c)  = eps % n(c)
          eps % oo(c) = eps % n(c)
          turb % y_plus(c) = 0.001
          if(Flow % heat_transfer) then
            vals(0) = t2_def;  t2 % n(c) = vals(Key_Ind('T2',  keys, nks))
            t2 % o(c)  = t2 % n(c)
            t2 % oo(c) = t2 % n(c)
          end if
        end if

        if(turb % model .eq. K_EPS_ZETA_F .or.  &
           turb % model .eq. HYBRID_LES_RANS) then
          vals(0) = kin_def;  kin  % n(c) = vals(Key_Ind('KIN',  keys, nks))
          vals(0) = eps_def;  eps  % n(c) = vals(Key_Ind('EPS',  keys, nks))
          vals(0) = zeta_def; zeta % n(c) = vals(Key_Ind('ZETA', keys, nks))
          vals(0) = f22_def;  f22  % n(c) = vals(Key_Ind('F22',  keys, nks))
          kin  % o(c)  = kin  % n(c)
          kin  % oo(c) = kin  % n(c)
          eps  % o(c)  = eps  % n(c)
          eps  % oo(c) = eps  % n(c)
          zeta % o(c)  = zeta % n(c)
          zeta % oo(c) = zeta % n(c)
          f22  % o(c)  = f22  % n(c)
          f22  % oo(c) = f22  % n(c)
          turb % y_plus(c) = 0.001
          if(Flow % heat_transfer) then
            vals(0) = t2_def;  t2 % n(c) = vals(Key_Ind('T2',  keys, nks))
            t2 % o(c)  = t2 % n(c)
            t2 % oo(c) = t2 % n(c)
          end if
        end if

        if(turb % model .eq. SPALART_ALLMARAS .or.  &
           turb % model .eq. DES_SPALART) then
          vals(0) = vis_def; vis % n(c) = vals(Key_Ind('VIS', keys, nks))
          vis % o(c)  = vis % n(c)
          vis % oo(c) = vis % n(c)
        end if

      end do ! through cells

    end if

  end if

  call User_Mod_Initialize_Variables(Flow, turb, Vof, swarm, Sol)

  !--------------------------------!
  !      Calculate the inflow      !
  !   and initializes the v_flux   !
  !   at both inflow and outflow   !
  !--------------------------------!
  n_wall        = 0
  n_inflow      = 0
  n_outflow     = 0
  n_symmetry    = 0
  n_heated_wall = 0
  n_convect     = 0
  n_pressure    = 0

  bulk % vol_in = 0.0
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2  < 0) then
      v_flux % n(s) = u % n(c2) * grid % sx(s)  &
                    + v % n(c2) * grid % sy(s)  &
                    + w % n(c2) * grid % sz(s)

      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
        bulk % vol_in = bulk % vol_in - v_flux % n(s)
        area = area  + grid % s(s)
      end if
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL)      &
        n_wall        = n_wall        + 1
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW)    &
        n_inflow      = n_inflow      + 1
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. OUTFLOW)   &
        n_outflow     = n_outflow     + 1
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. SYMMETRY)  &
        n_symmetry    = n_symmetry    + 1
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL)    &
        n_heated_wall = n_heated_wall + 1
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. CONVECT)   &
        n_convect     = n_convect     + 1
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. PRESSURE)  &
        n_pressure    = n_pressure    + 1
    else
      v_flux % n(s) = 0.0
    end if
  end do

  call Comm_Mod_Global_Sum_Int(n_wall)
  call Comm_Mod_Global_Sum_Int(n_inflow)
  call Comm_Mod_Global_Sum_Int(n_outflow)
  call Comm_Mod_Global_Sum_Int(n_symmetry)
  call Comm_Mod_Global_Sum_Int(n_heated_wall)
  call Comm_Mod_Global_Sum_Int(n_convect)
  call Comm_Mod_Global_Sum_Int(n_pressure)
  call Comm_Mod_Global_Sum_Real(bulk % vol_in)
  call Comm_Mod_Global_Sum_Real(area)

  ! This parameter, has_pressure_outlet, is used in Compute_Pressure
  Flow % has_pressure_outlet = .false.
  if(n_pressure > 0) then
    Flow % has_pressure_outlet = .true.
  end if

  !----------------------!
  !   Initializes time   !
  !----------------------!
  if(this_proc  < 2) then
    if(n_inflow .gt. 0) then
      print '(a29,es12.5)', ' # Volume inflow           : ', bulk % vol_in
      if (Vof % model .eq. VOLUME_OF_FLUID) then
        ! Needs to be corrected
        print '(a29,es12.5)', ' # Average inflow velocity : ',  &
          bulk % vol_in / area
      else
        print '(a29,es12.5)', ' # Average inflow velocity : ',  &
          bulk % vol_in / area
      end if
    end if
    print *, '# Number of faces on the wall        : ', n_wall
    print *, '# Number of inflow faces             : ', n_inflow
    print *, '# Number of outflow faces            : ', n_outflow
    print *, '# Number of symetry faces            : ', n_symmetry
    print *, '# Number of faces on the heated wall : ', n_heated_wall
    print *, '# Number of convective outflow faces : ', n_convect
    print *, '# Number of pressure outflow faces   : ', n_pressure
    print *, '# Variables initialized !'
  end if

  end subroutine
