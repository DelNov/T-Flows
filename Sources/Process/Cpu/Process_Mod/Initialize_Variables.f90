!==============================================================================!
  subroutine Initialize_Variables(Process, Flow, Turb, Vof, Swarm, Sol)
!------------------------------------------------------------------------------!
!>  This subroutine initializes dependent variables in the simulation. It is   !
!>  crucial to note that it only executes if a backup file was not read.       !
!>  Therefore, any variable initialized here and needed post-restart should be !
!>  stored in the backup file.                                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Process_Type)        :: Process  !! parent class
  type(Field_Type),   target :: Flow     !! flow object
  type(Turb_Type),    target :: Turb     !! turbulence object
  type(Vof_Type)             :: Vof      !! VOF object
  type(Swarm_Type)           :: Swarm    !! swarm object
  type(Solver_Type)          :: Sol      !! solver object
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: Grid
  type(Bulk_Type),  pointer :: bulk
  type(Var_Type),   pointer :: u, v, w, t, phi
  type(Var_Type),   pointer :: kin, eps, f22, zeta, vis, t2, omega
  type(Var_Type),   pointer :: uu, vv, ww, uv, uw, vw
  type(Face_Type),  pointer :: v_flux
  real,             pointer :: u_mean(:), v_mean(:), w_mean(:)
  integer                   :: i, c, c1, c2, m, s, nks, nvs, sc, fu
  integer                   :: n_wall, n_inflow, n_outflow, n_symmetry,  &
                               n_heated_wall, n_pressure, n_ambient, n_convect
  character(SL)             :: keys(128)
  character(SL)             :: keys_file(128)
  character(SL)             :: vals(0:128)   ! note that they start from zero!
  real                      :: area

  integer                   :: n_points, k
  real, allocatable         :: prof(:,:), x(:), y(:), z(:), dist(:)
  logical                   :: found, file_exists

  ! Default values for initial conditions
  character(3) :: u_def   = '0.0',  v_def    = '0.0',  w_def   = '0.0'
  character(3) :: t_def   = '0.0',  t2_def   = '0.0',  phi_def = '0.0'
  character(3) :: vf_def  = '0.0'
  character(3) :: kin_def = '0.0',  eps_def  = '0.0',  f22_def = '0.0'
  character(3) :: vis_def = '0.0',  zeta_def = '0.0',  omg_def = '1.0'
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
  omega    => Turb % omega
  call Flow % Alias_Momentum    (u, v, w)
  call Flow % Alias_Energy      (t)
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)
  call Turb % Alias_Stresses    (uu, vv, ww, uv, uw, vw)
  call Turb % Alias_T2          (t2)

  area  = 0.0
  O_Print '(a,a)', ' # Grid name: ', trim(Grid % name)

  ! Found the line where boundary condition definition is defined
  call Control % Position_At_One_Key('INITIAL_CONDITION', &
                                     found,               &
                                     .true.)

  !-----------------------------------------------!
  !                                               !
  !   Found the section with initial conditions   !
  !                                               !
  !-----------------------------------------------!
  if (found) then

    call Control % Read_Strings_On('VARIABLES', keys, nks, .true.)

    ! Input is valid, turn keys to upper case
    do i = 1, nks
      call String % To_Upper_Case(keys(i))
    end do

    ! Check if there is file specified
    call Control % Read_Strings_On('FILE', keys_file, nvs, .true.)

    !------------------------------------------------!
    !                                                !
    !   Initial conditions is prescribed in a file   !
    !                                                !
    !------------------------------------------------!
    if (nvs .eq. 1) then ! word 'file' was specified

      if (First_Proc()) &
        print *, '# Values specified in the file: ', trim(keys_file(nvs))

      !call File % Open_For_Reading_Ascii(keys_file(1), fu, This_Proc())
      call File % Open_For_Reading_Ascii(keys_file(1), fu)

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
        do c = Cells_In_Domain_And_Buffers()

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

          if(Turb % model .eq. K_OMEGA_SST) then
            i = Key_Ind('KIN',  keys, nks)
            read(kin_def,  *) prof(k, 0);  kin   % n(c) = prof(k, i)
            i = Key_Ind('OMG',  keys, nks)
            read(eps_def,  *) prof(k, 0);  omega % n(c) = prof(k, i)
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

        end do ! c = 1, Grid % n_cells

        call Global % Wait
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
      call Control % Position_At_One_Key('INITIAL_CONDITION', &
                                         found,               &
                                         .true.)

      call Control % Read_Strings_On('VARIABLES', keys, nks, .true.)

      ! Input is valid, turn keys to upper case
      do i = 1, nks
        call String % To_Upper_Case(keys(i))
      end do

      call Control % Read_Strings_On('VALUES', vals(1), nvs, .true.)

      ! Check validity of the input
      if(nks .eq. 0 .or. nvs .eq. 0) then
        call Message % Error(72,                                     &
                      'Critical, for initial condition: '//          &
                      'no values or variables have been provided ',  &
                      file=__FILE__, line=__LINE__, one_proc=.true.)
      end if
      if(nks .ne. nvs) then
        call Message % Error(72,                                               &
                      'Critical, for initial condition: number of values '//   &
                      'is not the same as number of provided variable names.', &
                      file=__FILE__, line=__LINE__, one_proc=.true.)
      end if

      ! Input is valid, turn keys to upper case
      do i = 1, nks
        call String % To_Upper_Case(keys(i))
      end do

      do c = Cells_In_Domain_And_Buffers()

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

        !---------!
        !   Vof   !
        !---------!
        if(Flow % with_interface) then
          vals(0) = vf_def
          read(vals(Key_Ind('VOF', keys, nks)), *, err=999)  Vof % fun % n(c)
 999      continue  ! file name may be defined
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

        if(Turb % model .eq. K_OMEGA_SST) then
          vals(0) = kin_def
          read(vals(Key_Ind('KIN', keys, nks)), *)  kin % n(c)
          kin % n(c) = max(0.01, kin % n(c))
          vals(0) = omg_def
          read(vals(Key_Ind('OMG', keys, nks)), *)  omega % n(c)
          omega % n(c) = max(0.001, omega % n(c))
          kin % o(c)  = kin % n(c)
          kin % oo(c) = kin % n(c)
          omega % o(c)  = omega % n(c)
          omega % oo(c) = omega % n(c)
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
      inquire(file = trim(Vof % name_stl), exist = file_exists)

      ! File exists
      if(file_exists) then
        call Vof % Initialize_From_Stl()
      endif
    end if

  end if

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
  n_ambient     = 0

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
      if(Grid % Bnd_Cond_Type(c2) .eq. AMBIENT)  &
        n_ambient     = n_ambient     + 1
    else
      v_flux % n(s) = 0.0
    end if
  end do

  call Global % Sum_Ints(n_wall,         &
                         n_inflow,       &
                         n_outflow,      &
                         n_symmetry,     &
                         n_heated_wall,  &
                         n_convect,      &
                         n_pressure,     &
                         n_ambient)
  call Global % Sum_Reals(bulk % vol_in, area)

  !---------------------------------------------!
  !   Parameters has_pressure and has_ambient   !
  !---------------------------------------------!
  Flow % has_pressure = .false.
  if(n_pressure .gt. 0) Flow % has_pressure = .true.

  Flow % has_ambient = .false.
  if(n_ambient .gt. 0) Flow % has_ambient = .true.
  Flow % reached_ambient_pressure = .false.

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
    print *, '# Number of ambient faces            : ', n_ambient
    print *, '# Variables initialized !'
  end if

  end subroutine
