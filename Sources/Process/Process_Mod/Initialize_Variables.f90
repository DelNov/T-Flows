!==============================================================================!
  subroutine Initialize_Variables(Process, Flow, Turb, Vof, Swarm, Sol)
!------------------------------------------------------------------------------!
!   Initialize dependent variables.  (It is a bit of a mess still)             !
!                                                                              !
!   It is important to remember that this procedure is called only if backup   !
!   file wasn't read.  Hence, if you initialize something here which should    !
!   be kept after restart, it should be stored in backup file.
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
  class(Process_Type)        :: Process
  type(Field_Type),   target :: Flow
  type(Turb_Type),    target :: Turb
  type(Vof_Type)             :: Vof
  type(Swarm_Type)           :: Swarm
  type(Solver_Type)          :: Sol
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Bulk_Type),  pointer :: bulk
  type(Var_Type),   pointer :: u, v, w, t, phi
  type(Var_Type),   pointer :: kin, eps, f22, zeta, vis, t2
  type(Var_Type),   pointer :: uu, vv, ww, uv, uw, vw
  type(Face_Type),  pointer :: v_flux
  real,             pointer :: u_mean(:), v_mean(:), w_mean(:)
  integer                   :: i, c, c1, c2, m, s, nks, nvs, sc, fu
  integer                   :: n_wall, n_inflow, n_outflow, n_symmetry,  &
                               n_heated_wall, n_pressure, n_convect
  character(SL)             :: keys(128)
  character(SL)             :: keys_file(128)
  character(SL)             :: vals(0:128)   ! note that they start from zero!
  real                      :: area

  integer                   :: n_points, k
  real, allocatable         :: prof(:,:), x(:), y(:), z(:), dist(:)
  logical                   :: found

  ! Default values for initial conditions
  character(3) :: u_def   = '0.0',  v_def    = '0.0',  w_def   = '0.0'
  character(3) :: t_def   = '0.0',  t2_def   = '0.0',  phi_def = '0.0'
  character(3) :: kin_def = '0.0',  eps_def  = '0.0',  f22_def = '0.0'
  character(3) :: vis_def = '0.0',  zeta_def = '0.0'
  character(3) :: uu_def  = '0.0',  vv_def   = '0.0',  ww_def  = '0.0'
  character(3) :: uv_def  = '0.0',  uw_def   = '0.0',  vw_def  = '0.0'
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Process)
!==============================================================================!

  ! Take aliases
  Grid     => Flow % pnt_grid
  bulk     => Flow % bulk
  v_flux   => Flow % v_flux
  vis      => Turb % vis
  u_mean   => Turb % u_mean
  v_mean   => Turb % v_mean
  w_mean   => Turb % w_mean
  call Flow % Alias_Momentum    (u, v, w)
  call Flow % Alias_Energy      (t)
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)
  call Turb % Alias_Stresses    (uu, vv, ww, uv, uw, vw)
  call Turb % Alias_T2          (t2)

  area  = 0.0
  if (First_Proc()) print '(a,a)', ' # Grid name: ', trim(Grid % name)

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
      call String % To_Upper_Case(keys(i))
    end do

    ! Check if there is file specified
    call Control_Mod_Read_Strings_On('FILE', keys_file, nvs, .true.)

    !------------------------------------------------!
    !                                                !
    !   Initial conditions is prescribed in a file   !
    !                                                !
    !------------------------------------------------!
    if (nvs .eq. 1) then ! word 'file' was specified

      if (First_Proc()) &
        print *, '# Values specified in the file: ', trim(keys_file(nvs))

      call File % Open_For_Reading_Ascii(keys_file(1), fu, This_Proc())

      ! Number of points
      call File % Read_Line(fu)

      read(Line % tokens(1), *) n_points

      if (First_Proc()) print '(a,i0,2a)', " # Reading ", nks, &
        " columns in file " , trim(keys_file(1))

      allocate(prof(n_points, 0:nks)); prof = 0.
      allocate(x(n_points));           x    = 0.
      allocate(y(n_points));           y    = 0.
      allocate(z(n_points));           z    = 0.
      allocate(dist(n_points));        dist = 0.

      ! Read the entire profile file
      do m = 1, n_points
        call File % Read_Line(fu)
        do i = 1, nks
          read(Line % tokens(i), *) prof(m, i)
        end do
      end do
      close(fu)

      ! A plane is defined
      if (keys(1) .eq. 'X' .and. keys(2) .eq. 'Y' .or.  &
          keys(1) .eq. 'X' .and. keys(2) .eq. 'Z' .or.  &
          keys(1) .eq. 'Y' .and. keys(2) .eq. 'Z') then

        ! Set the closest point
        do c = 1, Grid % n_cells

          i=Key_Ind('X', keys, nks); x(:) = prof(:,i)
          i=Key_Ind('Y', keys, nks); y(:) = prof(:,i)
          i=Key_Ind('Z', keys, nks); z(:) = prof(:,i)

          ! do no waste time on sqrt((r-r0)^2) -> use (r-r0)^2
          if(keys(1) .eq. 'Y' .and. keys(2) .eq. 'Z') then
            dist(:) = (y(:)-Grid % yc(c))**2 + (z(:)-Grid % zc(c))**2
          else if(keys(1) .eq. 'X' .and. keys(2) .eq. 'Z') then
            dist(:) = (x(:)-Grid % xc(c))**2 + (z(:)-Grid % zc(c))**2
          else if(keys(1) .eq. 'X' .and. keys(2) .eq. 'Y') then
            dist(:) = (x(:)-Grid % xc(c))**2 + (y(:)-Grid % yc(c))**2
          end if

          ! Store closest point in k
          k = minloc(dist, dim = 1)

          !--------------!
          !   Momentum   !
          !--------------!
          i=Key_Ind('U',keys,nks); read(u_def, *)  prof(k,0); u % n(c)=prof(k,i)
          i=Key_Ind('V',keys,nks); read(v_def, *)  prof(k,0); v % n(c)=prof(k,i)
          i=Key_Ind('W',keys,nks); read(w_def, *)  prof(k,0); w % n(c)=prof(k,i)

          !-------------------!
          !   Heat transfer   !
          !-------------------!
          if(Flow % heat_transfer) then
            i=Key_Ind('T',keys,nks); read(t_def,*) prof(k,0); t % n(c)=prof(k,i)
          end if

          !------------------------------!
          !   Scalars are missing here   !
          !------------------------------!

          !-------------------------!
          !   Vof is missing here   !
          !-------------------------!

          !--------------------------!
          !   Turbulent quantities   !
          !--------------------------!
          if(Turb % model .eq. K_EPS) then
            i = Key_Ind('KIN',  keys, nks)
            read(kin_def,  *) prof(k, 0);  kin  % n(c) = prof(k, i)
            i = Key_Ind('EPS',  keys, nks)
            read(eps_def,  *) prof(k, 0);  eps  % n(c) = prof(k, i)
            if(Flow % heat_transfer) then
              i = Key_Ind('T2', keys, nks)
              read(t2_def, *) prof(k,0);   t2 % n(c) = prof(k, i)
            end if
          end if

          if(Turb % model .eq. K_EPS_ZETA_F .or.  &
             Turb % model .eq. HYBRID_LES_RANS) then
            i = Key_Ind('KIN',  keys, nks)
            read(kin_def,  *) prof(k, 0);  kin  % n(c) = prof(k, i)
            i = Key_Ind('EPS',  keys, nks)
            read(eps_def,  *) prof(k, 0);  eps  % n(c) = prof(k, i)
            i = Key_Ind('ZETA', keys, nks)
            read(zeta_def, *) prof(k, 0);  zeta % n(c) = prof(k, i)
            i = Key_Ind('F22',  keys, nks)
            read(f22_def,  *) prof(k, 0);  f22  % n(c) = prof(k, i)
            if(Flow % heat_transfer) then
              i = Key_Ind('T2', keys, nks)
              read(t2_def, *) prof(k,0);   t2 % n(c) = prof(k, i)
            end if
           end if

          if(Turb % model .eq. DES_SPALART) then
            i = Key_Ind('VIS', keys, nks)
            read(vis_def, *) prof(k, 0);  vis % n(c) = prof(k, i)
          end if

          if(Turb % model .eq. RSM_MANCEAU_HANJALIC .or. &
             Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
            i = Key_Ind('UU',  keys, nks)
            read(uu_def, *) prof(k, 0);   uu % n(c)  = prof(k, i)
            i = Key_Ind('VV',  keys, nks)
            read(vv_def, *) prof(k, 0);   vv % n(c)  = prof(k, i)
            i = Key_Ind('WW',  keys, nks)
            read(ww_def, *) prof(k, 0);   ww % n(c)  = prof(k, i)
            i = Key_Ind('UV',  keys, nks)
            read(uv_def, *) prof(k, 0);   uv % n(c)  = prof(k, i)
            i = Key_Ind('UW',  keys, nks)
            read(uw_def, *) prof(k, 0);   uw % n(c)  = prof(k, i)
            i = Key_Ind('VW',  keys, nks)
            read(vw_def, *) prof(k, 0);   vw % n(c)  = prof(k, i)
            i = Key_Ind('EPS', keys, nks)
            read(eps_def, *) prof(k, 0);  eps % n(c) = prof(k, i)
            if (Turb % model .eq. RSM_MANCEAU_HANJALIC) then
              i = Key_Ind('F22', keys, nks)
              read(f22_def, *) prof(k, 0)
              f22 % n(c) = prof(k, i)
            end if
          end if

        end do ! c = 1, Grid % n_cells

        call Comm_Mod_Wait
        deallocate(prof)
        deallocate(x)
        deallocate(y)
        deallocate(z)
        deallocate(dist)

      end if

    !---------------------------------------------------!
    !                                                   !
    !   Initial conditions is NOT prescribed in a file  !
    !                                                   !
    !---------------------------------------------------!
    else

      ! Go back to key and read again
      call Control_Mod_Position_At_One_Key('INITIAL_CONDITION', &
                                           found,               &
                                           .true.)

      call Control_Mod_Read_Strings_On('VARIABLES', keys, nks, .true.)

      ! Input is valid, turn keys to upper case
      do i = 1, nks
        call String % To_Upper_Case(keys(i))
      end do

      call Control_Mod_Read_Strings_On('VALUES', vals(1), nvs, .true.)

      ! Check validity of the input
      if(nks .eq. 0 .or. nvs .eq. 0 .and. First_Proc()) then
        print '(2a)', '# Critical, for initial condition: ',        &
                      ' no values or variables have been provided'
        call Comm_Mod_End
        stop
      end if
      if(nks .ne. nvs .and. First_Proc()) then
        print '(2a)', '# Critical for initial conditions, number of values ',  &
                      ' is not the same as number of provided variable names'
        call Comm_Mod_End
        stop
      end if

      ! Input is valid, turn keys to upper case
      do i = 1, nks
        call String % To_Upper_Case(keys(i))
      end do

      do c = 1, Grid % n_cells

        if(Turb % statistics) then
          u_mean(c) = 0.0
          v_mean(c) = 0.0
          w_mean(c) = 0.0
        end if

        !--------------!
        !   Momentum   !
        !--------------!
        vals(0) = u_def;  read(vals(Key_Ind('U', keys, nks)), *)  u % n(c)
        vals(0) = v_def;  read(vals(Key_Ind('V', keys, nks)), *)  v % n(c)
        vals(0) = w_def;  read(vals(Key_Ind('W', keys, nks)), *)  w % n(c)

        u % o(c)  = u % n(c)
        u % oo(c) = u % n(c)
        v % o(c)  = v % n(c)
        v % oo(c) = v % n(c)
        w % o(c)  = w % n(c)
        w % oo(c) = w % n(c)

        !-------------------!
        !   Heat transfer   !
        !-------------------!
        if(Flow % heat_transfer) then
          vals(0) = t_def;  read(vals(Key_Ind('T', keys, nks)), *)  t % n(c)
          t % o(c)  = t % n(c)
          t % oo(c) = t % n(c)
        end if

        !-------------!
        !   Scalars   !
        !-------------!
        do sc = 1, Flow % n_scalars
          phi => Flow % scalar(sc)
          vals(0) = phi_def
          read(vals(Key_Ind(phi % name, keys, nks)), *) phi % n(c)
          phi % o(c)  = phi % n(c)
          phi % oo(c) = phi % n(c)
        end do

        !--------------------------!
        !   Turbulent quantities   !
        !--------------------------!

        ! Reynols stress models
        if(Turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
           Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
          vals(0) = uu_def;  read(vals(Key_Ind('UU', keys,nks)), *)  uu  % n(c)
          vals(0) = vv_def;  read(vals(Key_Ind('VV', keys,nks)), *)  vv  % n(c)
          vals(0) = ww_def;  read(vals(Key_Ind('WW', keys,nks)), *)  ww  % n(c)
          vals(0) = uv_def;  read(vals(Key_Ind('UV', keys,nks)), *)  uv  % n(c)
          vals(0) = uw_def;  read(vals(Key_Ind('UW', keys,nks)), *)  uw  % n(c)
          vals(0) = vw_def;  read(vals(Key_Ind('VW', keys,nks)), *)  vw  % n(c)
          vals(0) = eps_def; read(vals(Key_Ind('EPS',keys,nks)), *)  eps % n(c)
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
          if(Turb % model .eq. RSM_MANCEAU_HANJALIC) then
            vals(0) = f22_def;  read(vals(Key_Ind('F22',keys,nks)),*) f22 % n(c)
            f22 % o(c)  = f22 % n(c)
            f22 % oo(c) = f22 % n(c)
          end if
        end if

        if(Turb % model .eq. K_EPS) then
          vals(0) = kin_def
          read(vals(Key_Ind('KIN', keys, nks)), *)  kin % n(c)
          kin % n(c) = max(0.01, kin % n(c))
          vals(0) = eps_def
          read(vals(Key_Ind('EPS', keys, nks)), *)  eps % n(c)
          eps % n(c) = max(0.001, eps % n(c))
          kin % o(c)  = kin % n(c)
          kin % oo(c) = kin % n(c)
          eps % o(c)  = eps % n(c)
          eps % oo(c) = eps % n(c)
          Turb % y_plus(c) = 0.001
          if(Flow % heat_transfer) then
            vals(0) = t2_def;  read(vals(Key_Ind('T2', keys, nks)), *) t2 % n(c)
            t2 % o(c)  = t2 % n(c)
            t2 % oo(c) = t2 % n(c)
          end if
        end if

        if(Turb % model .eq. K_EPS_ZETA_F .or.  &
           Turb % model .eq. HYBRID_LES_RANS) then
          vals(0) = kin_def;  read(vals(Key_Ind('KIN', keys,nks)),*) kin  % n(c)
          vals(0) = eps_def;  read(vals(Key_Ind('EPS', keys,nks)),*) eps  % n(c)
          vals(0) = zeta_def; read(vals(Key_Ind('ZETA',keys,nks)),*) zeta % n(c)
          vals(0) = f22_def;  read(vals(Key_Ind('F22', keys,nks)),*) f22  % n(c)
          kin  % o(c)  = kin  % n(c)
          kin  % oo(c) = kin  % n(c)
          eps  % o(c)  = eps  % n(c)
          eps  % oo(c) = eps  % n(c)
          zeta % o(c)  = zeta % n(c)
          zeta % oo(c) = zeta % n(c)
          f22  % o(c)  = f22  % n(c)
          f22  % oo(c) = f22  % n(c)
          Turb % y_plus(c) = 0.001
          if(Flow % heat_transfer) then
            vals(0) = t2_def;  read(vals(Key_Ind('T2',keys,nks)),*) t2 % n(c)
            t2 % o(c)  = t2 % n(c)
            t2 % oo(c) = t2 % n(c)
          end if
        end if

        if(Turb % model .eq. SPALART_ALLMARAS .or.  &
           Turb % model .eq. DES_SPALART) then
          vals(0) = vis_def;  read(vals(Key_Ind('VIS',keys,nks)),*) vis % n(c)
          vis % o(c)  = vis % n(c)
          vis % oo(c) = vis % n(c)
        end if

      end do ! through cells

    end if

    !---------------------------------!
    !   Vof initialization from STL   !
    !---------------------------------!
    if(Flow % with_interface) then
      read(vals(Key_Ind('VOF', keys, nks)), *)  Vof % name_stl
      call Vof % Initialize_From_Stl()
    end if

  end if

  call User_Mod_Initialize_Variables(Flow, Turb, Vof, Swarm, Sol)

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
  do s = 1, Grid % n_faces
    c1 = Grid % faces_c(1,s)
    c2 = Grid % faces_c(2,s)
    if(c2  < 0) then
      v_flux % n(s) = u % n(c2) * Grid % sx(s)  &
                    + v % n(c2) * Grid % sy(s)  &
                    + w % n(c2) * Grid % sz(s)

      if(Grid % Bnd_Cond_Type(c2) .eq. INFLOW) then
        bulk % vol_in = bulk % vol_in - v_flux % n(s)
        area = area  + Grid % s(s)
      end if
      if(Grid % Bnd_Cond_Type(c2) .eq. WALL)      &
        n_wall        = n_wall        + 1
      if(Grid % Bnd_Cond_Type(c2) .eq. INFLOW)    &
        n_inflow      = n_inflow      + 1
      if(Grid % Bnd_Cond_Type(c2) .eq. OUTFLOW)   &
        n_outflow     = n_outflow     + 1
      if(Grid % Bnd_Cond_Type(c2) .eq. SYMMETRY)  &
        n_symmetry    = n_symmetry    + 1
      if(Grid % Bnd_Cond_Type(c2) .eq. WALLFL)    &
        n_heated_wall = n_heated_wall + 1
      if(Grid % Bnd_Cond_Type(c2) .eq. CONVECT)   &
        n_convect     = n_convect     + 1
      if(Grid % Bnd_Cond_Type(c2) .eq. PRESSURE)  &
        n_pressure    = n_pressure    + 1
    else
      v_flux % n(s) = 0.0
    end if
  end do

  call Global % Sum_Int(n_wall)
  call Global % Sum_Int(n_inflow)
  call Global % Sum_Int(n_outflow)
  call Global % Sum_Int(n_symmetry)
  call Global % Sum_Int(n_heated_wall)
  call Global % Sum_Int(n_convect)
  call Global % Sum_Int(n_pressure)
  call Global % Sum_Real(bulk % vol_in)
  call Global % Sum_Real(area)

  !----------------------------------------------------------------------!
  !   This parameter, has_pressure_outlet, is used in Compute_Pressure   !
  !----------------------------------------------------------------------!
  ! Update on July 17, 2021: I have some reservations about this part, since
  ! there was another bug fix when computing fluxes in the meanwhile (check:
  ! 90f77a1c8bd4ca05330a4435ed6321782ef00199).  This balancing also caused a
  ! bug when loading backup file (also check "Compute_Pressure" as well as
  ! "Backup_Mod/Load and Backup_Mod/Save" procedures)
  !
  ! Update on February 27, 2022: I have also added "has_outflow_boundary"
  ! to be able to tell PETSc if matrix for pressure is singular.  Shall
  ! it also be included in this test?
  !
  ! Update on June 2, 2022: Unified all outlet boundaries into one
  ! to be able to tell PETSc if matrix for pressure is singular
  Flow % has_pressure = .false.
  if(n_pressure > 0) Flow % has_pressure = .true.

  !----------------------!
  !   Initializes time   !
  !----------------------!
  if(First_Proc()) then
    if(n_inflow .gt. 0) then
      print '(a29,es12.5)', ' # Volume inflow           : ', bulk % vol_in
      if(Flow % with_interface) then
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
