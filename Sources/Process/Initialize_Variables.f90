!==============================================================================!
  subroutine Initialize_Variables(flow, bulk)
!------------------------------------------------------------------------------!
!   Initialize dependent variables.  (It is a bit of a mess still)             !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Field_Mod, only: Field_Type, heat_transfer, density
  use Les_Mod
  use Comm_Mod
  use Rans_Mod
  use Grid_Mod
  use Bulk_Mod
  use User_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  type(Bulk_Type)          :: bulk
!----------------------------------[Calling]-----------------------------------!
  integer :: Key_Ind
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, v, w, t
  real,            pointer :: flux(:)
  integer                  :: i, c, c1, c2, m, s, nks, nvs
  integer                  :: n_wall, n_inflow, n_outflow, n_symmetry,  &
                              n_heated_wall, n_convect
  character(len=80)        :: keys(128)
  character(len=80)        :: keys_file(128)
  real                     :: vals(0:128)   ! note that they start from zero!
  real                     :: area

  integer                  :: n_points, k
  real, allocatable        :: prof(:,:), x(:), y(:), z(:), dist(:)
  logical                  :: found

  ! Default values for initial conditions
  real, parameter          :: u_def   = 0.0,  v_def    = 0.0,  w_def    = 0.0
  real, parameter          :: t_def   = 0.0
  real, parameter          :: kin_def = 0.0,  eps_def  = 0.0,  f22_def  = 0.0
  real, parameter          :: vis_def = 0.0,  zeta_def = 0.0
  real, parameter          :: uu_def  = 0.0,  vv_def   = 0.0,  ww_def   = 0.0
  real, parameter          :: uv_def  = 0.0,  uw_def   = 0.0,  vw_def   = 0.0
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  flux => flow % flux
  u    => flow % u
  v    => flow % v
  w    => flow % w
  t    => flow % t

  area  = 0.0
  if (this_proc < 2) print *, '# Grid material: ', grid % material % name

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

      open(9, file = keys_file(1))
      if (this_proc < 2) print *, '# Reading the file: ', trim(keys_file(1))

      ! number of points
      call Tokenizer_Mod_Read_Line(9)

      read (line % tokens(1),*) n_points

      if (this_proc < 2) print '(A,I0,2A)', " # Reading ", nks, &
        " columns in file " , trim(keys_file(1))

      allocate(prof(n_points, 0:nks)); prof = 0.
      allocate(x(n_points));           x    = 0.
      allocate(y(n_points));           y    = 0.
      allocate(z(n_points));           z    = 0.
      allocate(dist(n_points));        dist = 0.

      ! Read the entire profile file
      do m = 1, n_points
        call Tokenizer_Mod_Read_Line(9)
        do i = 1, nks
          read(line % tokens(i), *) prof(m, i)
        end do
      end do
      close(9)

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

          if(heat_transfer) then
            i=Key_Ind('T',keys,nks);prof(k,0)=t_def;t%n(c)=prof(k,i)
          end if

          if(turbulence_model .eq. K_EPS) then
            i=Key_Ind('KIN',keys,nks);prof(k,0)=kin_def; kin%n(c)=prof(k,i)
            i=Key_Ind('EPS',keys,nks);prof(k,0)=eps_def; eps%n(c)=prof(k,i)
          end if

          if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
             turbulence_model .eq. HYBRID_LES_RANS) then
            i=Key_Ind('KIN', keys,nks);prof(k,0)=kin_def; kin %n(c)=prof(k,i)
            i=Key_Ind('EPS', keys,nks);prof(k,0)=eps_def; eps %n(c)=prof(k,i)
            i=Key_Ind('ZETA',keys,nks);prof(k,0)=zeta_def;zeta%n(c)=prof(k,i)
            i=Key_Ind('F22', keys,nks);prof(k,0)=f22_def; f22 %n(c)=prof(k,i)
          end if

          if(turbulence_model .eq. DES_SPALART) then
            i=Key_Ind('VIS',keys,nks); prof(k,0)=vis_def; vis%n(c)=prof(k,i)
          end if

          if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or. &
             turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
            i=Key_Ind('UU', keys,nks);prof(k,0)=uu_def; uu %n(c)=prof(k,i)
            i=Key_Ind('VV', keys,nks);prof(k,0)=vv_def; vv %n(c)=prof(k,i)
            i=Key_Ind('WW', keys,nks);prof(k,0)=ww_def; ww %n(c)=prof(k,i)
            i=Key_Ind('UV', keys,nks);prof(k,0)=uv_def; uv %n(c)=prof(k,i)
            i=Key_Ind('UW', keys,nks);prof(k,0)=uw_def; uw %n(c)=prof(k,i)
            i=Key_Ind('VW', keys,nks);prof(k,0)=vw_def; vw %n(c)=prof(k,i)
            i=Key_Ind('EPS',keys,nks);prof(k,0)=eps_def;eps%n(c)=prof(k,i)
            if (turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
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
      end if
      if(nks .ne. nvs .and. this_proc < 2) then
        print '(2a)', '# Critical for initial conditions, number of values ',  &
                      ' is not the same as number of provided variable names'
        call Comm_Mod_End
      end if

      ! Input is valid, turn keys to upper case
      do i = 1, nks
        call To_Upper_Case(keys(i))
      end do

      do c = 1, grid % n_cells

        if(turbulence_statistics) then
          u % mean(c) = 0.0
          v % mean(c) = 0.0
          w % mean(c) = 0.0
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

        if(heat_transfer) then
          vals(0) = t_def;  t % n(c) = vals(Key_Ind('T', keys, nks))
          t % o(c)  = t % n(c)
          t % oo(c) = t % n(c)
        end if

        if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
           turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
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
          if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
            vals(0) = f22_def; f22 % n(c) = vals(Key_Ind('F22', keys, nks))
            f22 % o(c)  = f22 % n(c)
            f22 % oo(c) = f22 % n(c)
          end if
        end if

        if(turbulence_model .eq. K_EPS) then
          vals(0) = kin_def; kin % n(c) = vals(Key_Ind('KIN', keys, nks))
          vals(0) = eps_def; eps % n(c) = vals(Key_Ind('EPS', keys, nks))
          kin % o(c)  = kin % n(c)
          kin % oo(c) = kin % n(c)
          eps % o(c)  = eps % n(c)
          eps % oo(c) = eps % n(c)
          u_tau(c)  = 0.047
          y_plus(c) = 0.001
        end if

        if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
           turbulence_model .eq. HYBRID_LES_RANS) then
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
          u_tau(c)  = 0.047
          y_plus(c) = 0.001
        end if

        if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
           turbulence_model .eq. DES_SPALART) then
          vals(0) = vis_def; vis % n(c) = vals(Key_Ind('VIS', keys, nks))
          vis % o(c)  = vis % n(c)
          vis % oo(c) = vis % n(c)
        end if

      end do ! through cells

    end if

  end if

  call User_Mod_Initialize(grid)

  !---------------------------------!
  !      Calculate the inflow       !
  !   and initializes the flux(s)   !
  !   at both inflow and outflow    !
  !---------------------------------!
  n_wall        = 0
  n_inflow      = 0
  n_outflow     = 0
  n_symmetry    = 0
  n_heated_wall = 0
  n_convect     = 0

  bulk % mass_in = 0.0
  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)
    if(c2  < 0) then
      flux(s) = density*( u % n(c2) * grid % sx(s) + &
                          v % n(c2) * grid % sy(s) + &
                          w % n(c2) * grid % sz(s) )

      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. INFLOW) then
        bulk % mass_in = bulk % mass_in - flux(s)
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
    else
      flux(s) = 0.0
    end if
  end do
  call Comm_Mod_Global_Sum_Int(n_wall)
  call Comm_Mod_Global_Sum_Int(n_inflow)
  call Comm_Mod_Global_Sum_Int(n_outflow)
  call Comm_Mod_Global_Sum_Int(n_symmetry)
  call Comm_Mod_Global_Sum_Int(n_heated_wall)
  call Comm_Mod_Global_Sum_Int(n_convect)
  call Comm_Mod_Global_Sum_Real(bulk % mass_in)
  call Comm_Mod_Global_Sum_Real(area)

  !----------------------!
  !   Initializes time   !
  !----------------------!
  if(this_proc  < 2) then
    if(n_inflow .gt. 0) then
      print '(a29,es12.5)', ' # Mass inflow             : ', bulk % mass_in
      print '(a29,es12.5)', ' # Average inflow velocity : ', bulk % mass_in &
                                                             / area
    end if
    print *, '# Number of faces on the wall        : ', n_wall
    print *, '# Number of inflow faces             : ', n_inflow
    print *, '# Number of outflow faces            : ', n_outflow
    print *, '# Number of symetry faces            : ', n_symmetry
    print *, '# Number of faces on the heated wall : ', n_heated_wall
    print *, '# Number of convective outflow faces : ', n_convect
    print *, '# Variables initialized !'
  end if

  end subroutine
