!==============================================================================!
  subroutine Load_Boundary_Conditions(grid, restart)
!------------------------------------------------------------------------------!
!   Reads boundary condition from control file                                 !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use Flow_Mod
  use Rans_Mod
  use Comm_Mod, only: this_proc
  use Tokenizer_Mod
  use Grid_Mod
  use User_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
  logical         :: restart
!----------------------------------[Calling]-----------------------------------!
  real    :: Distance
  integer :: Key_Ind
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, m, l, k, i, n, n_points, nks, nvs, us
  character(len=80) :: name_prof(128), answer, name_in
  real              :: wi, dist_min, x, y, z, xp, dist
  real, allocatable :: prof(:,:)
  logical           :: here
  character(len=80) :: bc_type_name, try_str
  integer           :: bc_type_tag
  character(len=80) :: keys(128)
  real              :: vals(0:128)           ! note that they start from zero!
  integer           :: types_per_color(128)  ! how many types in each color
  character(len=80) :: types_names(128)      ! name of each type
  logical           :: types_file(128)       ! type specified in a file?
  integer           :: c_types               ! counter types
  logical           :: found
!==============================================================================!

  !-----------------------------------------!
  ! Full name is specified in control file  !
  !-----------------------------------------!
  call Control_Mod_Load_Backup_Name(name_in)

  !-----------------------------------!
  ! Check if restart is present       !
  !-----------------------------------!
  answer = name_in
  call To_Upper_Case(answer)

  restart = .true.
  if(answer .eq. 'SKIP') then
    restart = .false.
  end if

  !----------------------------------------------------------------!
  !   Count number of types per boundary condition, total number   !
  !        of types specified, and also extract their names        !
  !----------------------------------------------------------------!
  types_per_color(:) = 0
  types_file(:)      = .false.
  c_types            = 0

  do n = 1, grid % n_bnd_cond
    call Control_Mod_Position_At_Two_Keys('BOUNDARY_CONDITION',       &
                                          grid % bnd_cond % name(n),  &
                                          found,                      &
                                          .false.)
    if(found) then
1     continue

      ! Try to read next 'TYPE' in the control file
      call Control_Mod_Read_Char_Item_On('TYPE', 'VOID', bc_type_name, .false.)

      ! Get out of the loop if you fail
      if(bc_type_name .eq. 'VOID') goto 2

      ! Skip following two lines
      call Control_Mod_Read_Char_Item_On('VARIABLES', 'VOID', try_str, .false.)
      call Control_Mod_Read_Char_Item_On('VALUES',    'VOID', try_str, .false.)

      types_per_color(n) = types_per_color(n) + 1
      c_types = c_types + 1
      types_names(c_types) = bc_type_name

      ! If try_str is 'VOID', it didn't find 'VALUES'
      ! meaning that the keyword 'FILE' was specified
      if(try_str .eq. 'VOID') then
        types_file(c_types) = .true.
      end if

      goto 1
    else
      print *, '# Boundary conditions are not specified in control file!'
      print *, '# Exiting the program.'
      stop
    end if

2 continue

  end do

  !------------------------------------------------!
  !                                                !
  !                                                !
  !   Read boundary conditions from control file   !
  !                                                !
  !                                                !
  !------------------------------------------------!
  c_types = 0

  do n = 1, grid % n_bnd_cond

    ! Position yourself well
    call Control_Mod_Position_At_Two_Keys('BOUNDARY_CONDITION',       &
                                          grid % bnd_cond % name(n),  &
                                          found,                      &
                                          .false.)
    do l = 1, types_per_color(n)

      ! Update the counter
      c_types = c_types + 1

      !---------------------------------------------!
      !                                             !
      !   Read first line which is common for all   !
      !                                             !
      !---------------------------------------------!
      call Control_Mod_Read_Char_Item_On('TYPE', 'WALL', bc_type_name, .false.)
      call To_Upper_Case(bc_type_name)

      ! Copy boundary conditions which were given for the grid
      if( bc_type_name .eq. 'INFLOW') then
        bc_type_tag = INFLOW
        grid % bnd_cond % type(n) = INFLOW
      else if( bc_type_name .eq. 'WALL') then
        bc_type_tag = WALL
        grid % bnd_cond % type(n) = WALL
      else if( bc_type_name .eq. 'OUTFLOW') then
        bc_type_tag = OUTFLOW
        grid % bnd_cond % type(n) = OUTFLOW
      else if( bc_type_name .eq. 'SYMMETRY') then
        bc_type_tag = SYMMETRY
        grid % bnd_cond % type(n) = SYMMETRY
      else if( bc_type_name .eq. 'WALL_FLUX') then
        bc_type_tag = WALLFL
        grid % bnd_cond % type(n) = WALLFL
      else if( bc_type_name .eq. 'CONVECTIVE') then
        bc_type_tag = CONVECT
        grid % bnd_cond % type(n) = CONVECT
      else if( bc_type_name .eq. 'PRESSURE') then
        bc_type_tag = PRESSURE
        grid % bnd_cond % type(n) = PRESSURE
      else
        if(this_proc < 2)  &
          print *, '# Load_Boundary_Conditions: '//        &
                   '# Unknown boundary condition type: ',  &
                   bc_type_name
        stop
      end if

      !----------------------------------------------!
      !                                              !
      !   Read second line which is common for all   !
      !                                              !
      !----------------------------------------------!
      call Control_Mod_Read_Strings_On('VARIABLES', keys, nks, .false.)
      do i = 1, nks
        call To_Upper_Case(keys(i))
      end do

      !-----------------------------------------------------------------!
      !                                                                 !
      !   Boundary values are specified in a list (and not in a file)   !
      !                                                                 !
      !-----------------------------------------------------------------!
      if( .not. types_file(c_types) ) then
        call Control_Mod_Read_Real_Array_On('VALUES', vals(1), nvs, .false.)

        !--------------------------------------------------!
        !   Distribute boundary values to boundary cells   !
        !--------------------------------------------------!

        ! Distribute b.c. tags only.
        ! They are distributed if restart or not
        do c = -1, -grid % n_bnd_cells, -1
          if(grid % bnd_cond % color(c) .eq. n) then

            ! Temperature
            if(heat_transfer) then
              i = Key_Ind('T', keys, nks)
              if(i > 0) t % bnd_cell_type(c) = bc_type_tag
              i = Key_Ind('Q', keys, nks)
              if(i > 0) t % bnd_cell_type(c) = bc_type_tag
            end if

            ! For user scalars
            do us = 1, n_user_scalars
              i = Key_Ind(user_scalar(us) % name, keys, nks)
              if(i > 0) user_scalar(us) % bnd_cell_type(c) = bc_type_tag
              i = Key_Ind(user_scalar(us) % flux_name, keys, nks)
              if(i > 0) user_scalar(us) % bnd_cell_type(c) = bc_type_tag
            end do

          end if  ! bnd_color .eq. n

        end do

        ! Distribute b.c. values
        ! They are distributed only if not restart
        do c = -1, -grid % n_bnd_cells, -1
          if(grid % bnd_cond % color(c) .eq. n .and. .not. restart) then

            ! For velocity and pressure
            i = Key_Ind('U', keys, nks); if(i > 0) u % n(c) = vals(i)
            i = Key_Ind('V', keys, nks); if(i > 0) v % n(c) = vals(i)
            i = Key_Ind('W', keys, nks); if(i > 0) w % n(c) = vals(i)
            i = Key_Ind('P', keys, nks); if(i > 0) p % n(c) = vals(i)

            ! Temperature
            if(heat_transfer) then
              i = Key_Ind('T', keys, nks)
              if(i > 0) t % n(c) = vals(i)
              i = Key_Ind('Q', keys, nks)
              if(i > 0) t % q(c) = vals(i)
            end if

            ! For user scalars
            do us = 1, n_user_scalars
              i = Key_Ind(user_scalar(us) % name, keys, nks)
              if(i > 0) user_scalar(us) % n(c) = vals(i)
              i = Key_Ind(user_scalar(us) % flux_name, keys, nks)
              if(i > 0) user_scalar(us) % q(c) = vals(i)
            end do

            ! For turbulence models
            if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
               turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
              i = Key_Ind('UU',  keys, nks); if(i > 0) uu  % n(c) = vals(i)
              i = Key_Ind('VV',  keys, nks); if(i > 0) vv  % n(c) = vals(i)
              i = Key_Ind('WW',  keys, nks); if(i > 0) ww  % n(c) = vals(i)
              i = Key_Ind('UV',  keys, nks); if(i > 0) uv  % n(c) = vals(i)
              i = Key_Ind('UW',  keys, nks); if(i > 0) uw  % n(c) = vals(i)
              i = Key_Ind('VW',  keys, nks); if(i > 0) vw  % n(c) = vals(i)
              i = Key_Ind('EPS', keys, nks); if(i > 0) eps % n(c) = vals(i)

              if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
                i = Key_Ind('F22', keys, nks); if(i > 0) f22 % n(c) = vals(i)
              end if
            end if

            if(turbulence_model .eq. K_EPS) then
              i = Key_Ind('KIN', keys, nks); if(i > 0) kin % n(c) = vals(i)
              i = Key_Ind('EPS', keys, nks); if(i > 0) eps % n(c) = vals(i)
              u_tau(c)  = 0.047
              y_plus(c) = 1.1
            end if

            if(turbulence_model .eq. K_EPS_ZETA_F) then
              i = Key_Ind('KIN',  keys, nks); if(i > 0) kin  % n(c) = vals(i)
              i = Key_Ind('EPS',  keys, nks); if(i > 0) eps  % n(c) = vals(i)
              i = Key_Ind('ZETA', keys, nks); if(i > 0) zeta % n(c) = vals(i)
              i = Key_Ind('F22',  keys, nks); if(i > 0) f22  % n(c) = vals(i)
            end if

            if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
               turbulence_model .eq. DES_SPALART) then
              i = Key_Ind('VIS',  keys, nks); if(i > 0) vis % n(c) = vals(i)
            end if
          end if ! restart
        end do

      !---------------------------------------------!
      !                                             !
      !   Boundary values are specified in a file   !
      !                                             !
      !---------------------------------------------!
      else

        call Control_Mod_Read_Strings_On('FILE', name_prof, nvs, .false.)

        open(9, file=name_prof(1))
        if(this_proc < 2) print *, '# Reading the file: ', trim(name_prof(1))
        call Tokenizer_Mod_Read_Line(9)
        read(line % tokens(1),*) n_points  ! number of points

        allocate(prof(n_points, 0:nks))

        !----------------------------------!
        !   Read the entire profile file   !
        !----------------------------------!
        do m = 1, n_points
          call Tokenizer_Mod_Read_Line(9)
          do i = 1, nks
            read(line % tokens(i), *) prof(m, i)
          end do
        end do
        close(9)

        !------------------------!
        !   A plane is defined   !
        !------------------------!
        if(keys(1) .eq. 'X' .and. keys(2) .eq. 'Y' .or.  &
           keys(1) .eq. 'X' .and. keys(2) .eq. 'Z' .or.  &
           keys(1) .eq. 'Y' .and. keys(2) .eq. 'Z') then

          ! Set the closest point
          do c = -1, -grid % n_bnd_cells, -1

            ! Distribute b.c. types
            ! They are distributed only if not restart
            if(grid % bnd_cond % color(c) .eq. n) then

              ! For temperature
              if(heat_transfer) then
                i = Key_Ind('T', keys, nks)
                if(i > 0) t % bnd_cell_type(c) = bc_type_tag
                i = Key_Ind('Q', keys, nks)
                if(i > 0) t % bnd_cell_type(c) = bc_type_tag
              end if

              ! For user scalars
              do us = 1, n_user_scalars
                i = Key_Ind(user_scalar(us) % name, keys, nks)
                if(i > 0) user_scalar(us) % bnd_cell_type(c) = bc_type_tag
                i = Key_Ind(user_scalar(us) % flux_name, keys, nks)
                if(i > 0) user_scalar(us) % bnd_cell_type(c) = bc_type_tag
              end do

            end if

            ! Distribute b.c. values
            ! They are distributed only if not restart
            if(grid % bnd_cond % color(c) .eq. n .and. .not. restart) then

              dist_min = HUGE
              do m = 1, n_points

                i = Key_Ind('X', keys, nks); prof(m,0) = 0.0;  x = prof(m,i)
                i = Key_Ind('Y', keys, nks); prof(m,0) = 0.0;  y = prof(m,i)
                i = Key_Ind('Z', keys, nks); prof(m,0) = 0.0;  z = prof(m,i)

                if(keys(1) .eq. 'Y' .and. keys(2) .eq. 'Z') then
                  dist = Distance(y,            z,            0.0,  &
                                  grid % yc(c), grid % zc(c), 0.0)

                else if(keys(1) .eq. 'X' .and. keys(2) .eq. 'Z') then
                  dist = Distance(x,            z,            0.0,  &
                                  grid % xc(c), grid % zc(c), 0.0)

                else if(keys(1) .eq. 'X' .and. keys(2) .eq. 'Y') then
                  dist = Distance(x,            y,            0.0,  &
                                  grid % xc(c), grid % yc(c), 0.0)

                end if

                ! Store closest point in k
                if(dist < dist_min) then
                  dist_min = dist
                  k = m
                end if

              end do

              ! For velocity and pressure
              i = Key_Ind('U', keys, nks); if(i > 0) u % n(c) = prof(k,i)
              i = Key_Ind('V', keys, nks); if(i > 0) v % n(c) = prof(k,i)
              i = Key_Ind('W', keys, nks); if(i > 0) w % n(c) = prof(k,i)
              i = Key_Ind('P', keys, nks); if(i > 0) p % n(c) = prof(k,i)

              ! For temperature
              if(heat_transfer) then
                i = Key_Ind('T', keys, nks)
                if(i > 0) t % n(c) = prof(k,i)
                i = Key_Ind('Q', keys, nks)
                if(i > 0) t % q(c) = prof(k,i)
              end if

              ! For user scalars
              do us = 1, n_user_scalars
                i = Key_Ind(user_scalar(us) % name, keys, nks)
                if(i > 0) user_scalar(us) % n(c) = prof(k,i)
                i = Key_Ind(user_scalar(us) % flux_name, keys, nks)
                if(i > 0) user_scalar(us) % q(c) = prof(k,i)
              end do

              ! For turbulence models
              if(turbulence_model .eq. K_EPS) then
                i = Key_Ind('KIN', keys, nks); if(i > 0) kin % n(c) = prof(k,i)
                i = Key_Ind('EPS', keys, nks); if(i > 0) eps % n(c) = prof(k,i)
              end if

              if(turbulence_model .eq. K_EPS_ZETA_F) then
                i = Key_Ind('KIN',  keys, nks); if(i>0) kin  % n(c) = prof(k,i)
                i = Key_Ind('EPS',  keys, nks); if(i>0) eps  % n(c) = prof(k,i)
                i = Key_Ind('ZETA', keys, nks); if(i>0) zeta % n(c) = prof(k,i)
                i = Key_Ind('F22',  keys, nks); if(i>0) f22  % n(c) = prof(k,i)
              end if

              if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
                 turbulence_model .eq. DES_SPALART) then
                i = Key_Ind('VIS', keys, nks); if(i > 0) vis % n(c) = prof(k,i)
              end if

              if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
                 turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
                i = Key_Ind('UU', keys, nks); if(i > 0) uu  % n(c) = prof(k,i)
                i = Key_Ind('VV', keys, nks); if(i > 0) vv  % n(c) = prof(k,i)
                i = Key_Ind('WW', keys, nks); if(i > 0) ww  % n(c) = prof(k,i)
                i = Key_Ind('UV', keys, nks); if(i > 0) uv  % n(c) = prof(k,i)
                i = Key_Ind('UW', keys, nks); if(i > 0) uw  % n(c) = prof(k,i)
                i = Key_Ind('VW', keys, nks); if(i > 0) vw  % n(c) = prof(k,i)
                i = Key_Ind('EPS',keys, nks); if(i > 0) eps % n(c) = prof(k,i)

                if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
                  i = Key_Ind('F22', keys, nks); if(i>0) f22 % n(c) = prof(k,i)
                end if
              end if
            end if      !end if(grid % bnd_cond % color(c) .eq. n .and. restart)
          end do        !end do c = -1, -grid % n_bnd_cells, -1

        !----------------------------!
        !   A plane is not defined   !
        !----------------------------!
        else  ! dir .eq. "XPL" ...

          do c = -1, -grid % n_bnd_cells, -1

            ! If restart is set to true, set boundary values,
            ! otherwise, just the TypeBC remains set.
            if(grid % bnd_cond % color(c) .eq. n) then

              do m = 1, n_points-1
                here = .false.

                i = Key_Ind(keys(1), keys, nks)
                prof(m,   0) = 0.0;
                prof(m+1, 0) = 0.0;
                x  = prof(m,i)
                xp = prof(m+1,i)

                ! Compute the weight factors
                if( keys(1) .eq. 'X' .and.  &
                    grid % xc(c) >= x .and. grid % xc(c) <= xp ) then
                  wi = ( xp - grid % xc(c) ) / (xp - x)
                  here = .true.
                else if( keys(1) .eq. 'Y' .and.  &
                         grid % yc(c) >= x .and. grid % yc(c) <= xp ) then
                  wi = ( xp - grid % yc(c) ) / (xp - x)
                  here = .true.
                else if( keys(1) .eq. 'Z' .and.  &
                         grid % zc(c) >= x .and. grid % zc(c) <= xp ) then
                  wi = ( xp - grid % zc(c) ) / (xp - x)
                  here = .true.

                ! Beware; for cylindrical coordinates you have "inversion"
                else if( (keys(1) .eq. 'RX' .and.  &
                     sqrt(grid % yc(c)**2 + grid % zc(c)**2) >= xp .and.       &
                     sqrt(grid % yc(c)**2 + grid % zc(c)**2) <= x) ) then
                  wi = ( xp - sqrt(grid % yc(c)**2 + grid % zc(c)**2) ) / (xp-x)
                  here = .true.
                else if( (keys(1) .eq. 'RY' .and.  &
                     sqrt(grid % xc(c)**2 + grid % zc(c)**2) >= xp .and.       &
                     sqrt(grid % xc(c)**2 + grid % zc(c)**2) <= x) ) then
                  wi = ( xp - sqrt(grid % xc(c)**2 + grid % zc(c)**2) ) / (xp-x)
                  here = .true.
                else if( (keys(1) .eq. 'RZ' .and.  &
                     sqrt(grid % xc(c)**2 + grid % yc(c)**2) >= xp .and.       &
                     sqrt(grid % xc(c)**2 + grid % yc(c)**2) <= x) ) then
                  wi = ( xp - sqrt(grid % xc(c)**2 + grid % yc(c)**2) ) / (xp-x)
                  here = .true.
                end if

                if(here) then

                  ! For temperature
                  if(heat_transfer) then
                    i = Key_Ind('T',keys,nks)
                    if(i > 0) t % bnd_cell_type(c) = bc_type_tag
                    i = Key_Ind('Q',keys,nks)
                    if(i > 0) t % bnd_cell_type(c) = bc_type_tag
                  end if

                  ! For user scalars
                  do us = 1, n_user_scalars
                    i = Key_Ind(user_scalar(us) % name, keys, nks)
                    if(i > 0) user_scalar(us) % bnd_cell_type(c) = bc_type_tag
                    i = Key_Ind(user_scalar(us) % flux_name, keys, nks)
                    if(i > 0) user_scalar(us) % bnd_cell_type(c) = bc_type_tag
                  end do

                end if  ! here
              end do    ! m, points
            end if      ! bnd_color .eq. n

            ! If restart is set to true, set boundary values,
            ! otherwise, just the TypeBC remains set.
            if(grid % bnd_cond % color(c) .eq. n .and. .not. restart) then

              do m = 1, n_points-1
                here = .false.

                i = Key_Ind(keys(1), keys, nks)
                prof(m,   0) = 0.0;
                prof(m+1, 0) = 0.0;
                x  = prof(m,i)
                xp = prof(m+1,i)

                ! Compute the weight factors
                if( keys(1) .eq. 'X' .and.  &
                    grid % xc(c) >= x .and. grid % xc(c) <= xp ) then
                  wi = ( xp - grid % xc(c) ) / (xp - x)
                  here = .true.
                else if( keys(1) .eq. 'Y' .and.  &
                         grid % yc(c) >= x .and. grid % yc(c) <= xp ) then
                  wi = ( xp - grid % yc(c) ) / (xp - x)
                  here = .true.
                else if( keys(1) .eq. 'Z' .and.  &
                         grid % zc(c) >= x .and. grid % zc(c) <= xp ) then
                  wi = ( xp - grid % zc(c) ) / (xp - x)
                  here = .true.

                ! Beware; for cylindrical coordinates you have "inversion"
                else if( (keys(1) .eq. 'RX' .and.  &
                     sqrt(grid % yc(c)**2 + grid % zc(c)**2) >= xp .and.       &
                     sqrt(grid % yc(c)**2 + grid % zc(c)**2) <= x) ) then
                  wi = ( xp - sqrt(grid % yc(c)**2 + grid % zc(c)**2) ) / (xp-x)
                  here = .true.
                else if( (keys(1) .eq. 'RY' .and.  &
                     sqrt(grid % xc(c)**2 + grid % zc(c)**2) >= xp .and.       &
                     sqrt(grid % xc(c)**2 + grid % zc(c)**2) <= x) ) then
                  wi = ( xp - sqrt(grid % xc(c)**2 + grid % zc(c)**2) ) / (xp-x)
                  here = .true.
                else if( (keys(1) .eq. 'RZ' .and.  &
                     sqrt(grid % xc(c)**2 + grid % yc(c)**2) >= xp .and.       &
                     sqrt(grid % xc(c)**2 + grid % yc(c)**2) <= x) ) then
                  wi = ( xp - sqrt(grid % xc(c)**2 + grid % yc(c)**2) ) / (xp-x)
                  here = .true.
                end if

                ! Interpolate the profiles
                if(here) then

                  ! For velocity and pressure
                  i = Key_Ind('U',keys,nks)
                  if(i > 0) u % n(c) = wi*prof(m, i) + (1.-wi)*prof(m+1, i)
                  i = Key_Ind('V',keys,nks)
                  if(i > 0) v % n(c) = wi*prof(m, i) + (1.-wi)*prof(m+1, i)
                  i = Key_Ind('W',keys,nks)
                  if(i > 0) w % n(c) = wi*prof(m, i) + (1.-wi)*prof(m+1, i)
                  i = Key_Ind('P',keys,nks)
                  if(i > 0) p % n(c) = wi*prof(m, i) + (1.-wi)*prof(m+1, i)

                  ! For temperature
                  if(heat_transfer) then
                    i = Key_Ind('T',keys,nks)
                    if(i > 0) t % n(c) = wi*prof(m, i) + (1.-wi)*prof(m+1, i)
                    if(i > 0) t % bnd_cell_type(c) = bc_type_tag
                    i = Key_Ind('Q',keys,nks)
                    if(i > 0) t % q(c) = wi*prof(m, i) + (1.-wi)*prof(m+1, i)
                    if(i > 0) t % bnd_cell_type(c) = bc_type_tag
                  end if

                  ! For user scalars
                  do us = 1, n_user_scalars
                    i = Key_Ind(user_scalar(us) % name, keys, nks)
                    if(i > 0) &
                      user_scalar(us) % n(c)=wi*prof(m,i)+(1.-wi)*prof(m+1,i)
                    i = Key_Ind(user_scalar(us) % flux_name, keys, nks)
                    if(i > 0) &
                      user_scalar(us) % q(c)=wi*prof(m,i)+(1.-wi)*prof(m+1,i)
                  end do

                  ! For turbulence models
                  if(turbulence_model .eq. K_EPS) then

                    i = Key_Ind('KIN',keys,nks)
                    if(i > 0) kin % n(c) = wi*prof(m, i) + (1.-wi)*prof(m+1, i)

                    i = Key_Ind('EPS',keys,nks)
                    if(i > 0) eps % n(c) = wi*prof(m, i) + (1.-wi)*prof(m+1, i)
                  end if

                  if(turbulence_model .eq. K_EPS_ZETA_F) then

                    i = Key_Ind('KIN',keys,nks)
                    if(i > 0) kin % n(c) = wi*prof(m, i) + (1.-wi)*prof(m+1, i)

                    i = Key_Ind('EPS',keys,nks)
                    if(i > 0) eps % n(c) = wi*prof(m, i) + (1.-wi)*prof(m+1, i)

                    i = Key_Ind('ZETA',keys,nks)
                    if(i > 0) zeta % n(c) = wi*prof(m, i) + (1.-wi)*prof(m+1, i)

                    i = Key_Ind('F22',keys,nks)
                    if(i > 0) f22 % n(c) = wi*prof(m, i) + (1.-wi)*prof(m+1, i)
                  end if

                  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
                     turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then

                    i = Key_Ind('UU', keys, nks)
                    if(i > 0) uu % n(c) = wi*prof(m, i) + (1.-wi)*prof(m+1, i)

                    i = Key_Ind('VV', keys, nks)
                    if(i > 0) vv % n(c) = wi*prof(m, i) + (1.-wi)*prof(m+1, i)

                    i = Key_Ind('WW', keys, nks)
                    if(i > 0) ww % n(c) = wi*prof(m, i) + (1.-wi)*prof(m+1, i)

                    i = Key_Ind('UV', keys, nks)
                    if(i > 0) uv % n(c) = wi*prof(m, i) + (1.-wi)*prof(m+1, i)

                    i = Key_Ind('UW', keys, nks)
                    if(i > 0) uw % n(c) = wi*prof(m, i) + (1.-wi)*prof(m+1, i)

                    i = Key_Ind('VW', keys, nks)
                    if(i > 0) vw % n(c) = wi*prof(m, i) + (1.-wi)*prof(m+1, i)

                    i = Key_Ind('EPS', keys, nks)
                    if(i > 0) eps % n(c) = wi*prof(m, i) + (1.-wi)*prof(m+1, i)

                    if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
                      i = Key_Ind('F22', keys, nks)
                      if(i > 0)f22 % n(c) = wi*prof(m, i) + (1.-wi)*prof(m+1, i)
                    end if
                  end if

                  if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
                     turbulence_model .eq. DES_SPALART) then
                    i = Key_Ind('VIS',keys,nks)
                    if(i > 0) vis % n(c) = wi*prof(m, i) + (1.-wi)*prof(m+1, i)
                  end if

                end if  ! (here)
              end do  ! m = 1, n_points-1
            end if  ! if(restart)
          end do  ! c = -1, -grid % n_bnd_cells, -1
        end if  ! plane is defined?
        close(9)
      end if  ! boundary defined in a file
    end do

  end do

  end subroutine
