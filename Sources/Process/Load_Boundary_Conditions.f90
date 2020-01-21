!==============================================================================!
  subroutine Load_Boundary_Conditions(flow, turb, mult, turb_planes)
!------------------------------------------------------------------------------!
!   Reads boundary condition from control file                                 !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Const_Mod
  use File_Mod
  use Comm_Mod,       only: this_proc, Comm_Mod_End
  use Field_Mod,      only: Field_Type, heat_transfer
  use Turb_Mod
  use Multiphase_Mod, only: Multiphase_Type, multiphase_model, VOLUME_OF_FLUID
  use Grid_Mod,       only: Grid_Type
  use Eddies_Mod
  use User_Mod
  use Control_Mod
  use Var_Mod,        only: Var_Type
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Turb_Plane_Type)         :: turb_planes
!----------------------------------[Calling]-----------------------------------!
  integer :: Key_Ind
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: grid
  type(Var_Type),  pointer :: u, v, w, t, p, vof
  type(Var_Type),  pointer :: kin, eps, f22, zeta, vis, t2
  type(Var_Type),  pointer :: uu, vv, ww, uv, uw, vw
  type(Var_Type),  pointer :: scalar(:)
  integer                  :: c,m,l,k,i,bc,n_points,nks,nvs,sc,c1,c2,s,fu
  character(len=80)        :: name_prof(128)
  real                     :: wi, dist_min, x, y, z, xp, dist
  real, allocatable        :: prof(:,:)
  logical                  :: here
  character(len=80)        :: bc_type_name, try_str
  integer                  :: bc_type_tag
  character(len=80)        :: keys(128)
  real                     :: vals(0:128)           ! they start from zero!
  integer                  :: types_per_color(128)  ! how many types in a color
  character(len=80)        :: types_names(128)      ! name of each type
  logical                  :: types_file(128)       ! type specified in a file?
  integer                  :: c_types               ! counter types
  integer                  :: edd_n
  real                     :: edd_r
  real                     :: edd_i
  logical                  :: found
!==============================================================================!

  ! Take aliases
  grid   => flow % pnt_grid
  t      => flow % t
  p      => flow % p
  scalar => flow % scalar
  vis    => turb % vis
  t2     => turb % t2
  vof    => mult % vof

  call Field_Mod_Alias_Momentum   (flow, u, v, w)
  call Turb_Mod_Alias_K_Eps_Zeta_F(turb, kin, eps, zeta, f22)
  call Turb_Mod_Alias_Stresses    (turb, uu, vv, ww, uv, uw, vw)

  !-------------------------!
  !   Read wall roughness   !
  !-------------------------!
  call Control_Mod_Roughness_Coefficient(turb % z_o)

  !----------------------------------------------------------------!
  !   Count number of types per boundary condition, total number   !
  !        of types specified, and also extract their names        !
  !----------------------------------------------------------------!
  types_per_color(:) = 0
  types_file(:)      = .false.
  c_types            = 0

  do bc = 1, grid % n_bnd_cond
    call Control_Mod_Position_At_Two_Keys('BOUNDARY_CONDITION',        &
                                          grid % bnd_cond % name(bc),  &
                                          found,                       &
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

      types_per_color(bc) = types_per_color(bc) + 1
      c_types = c_types + 1
      types_names(c_types) = bc_type_name

      ! If try_str is 'VOID', it didn't find 'VALUES'
      ! meaning that the keyword 'FILE' was specified
      if(try_str .eq. 'VOID') then
        types_file(c_types) = .true.
      end if

      goto 1
    else
      if(this_proc < 2) then
        print *, '# ERROR!  Boundary conditions for ',  &
                 trim(grid % bnd_cond % name(bc)),      &
                 ' not specified in the control file!'
        print *, '# Exiting the program.'
      end if
      call Comm_Mod_End
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

  do bc = 1, grid % n_bnd_cond

    ! Position yourself well
    call Control_Mod_Position_At_Two_Keys('BOUNDARY_CONDITION',        &
                                          grid % bnd_cond % name(bc),  &
                                          found,                       &
                                          .false.)
    do l = 1, types_per_color(bc)

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
        grid % bnd_cond % type(bc) = INFLOW
      else if( bc_type_name .eq. 'WALL') then
        bc_type_tag = WALL
        grid % bnd_cond % type(bc) = WALL
      else if( bc_type_name .eq. 'OUTFLOW') then
        bc_type_tag = OUTFLOW
        grid % bnd_cond % type(bc) = OUTFLOW
      else if( bc_type_name .eq. 'SYMMETRY') then
        bc_type_tag = SYMMETRY
        grid % bnd_cond % type(bc) = SYMMETRY
      else if( bc_type_name .eq. 'WALL_FLUX') then
        bc_type_tag = WALLFL
        grid % bnd_cond % type(bc) = WALLFL
      else if( bc_type_name .eq. 'CONVECTIVE') then
        bc_type_tag = CONVECT
        grid % bnd_cond % type(bc) = CONVECT
      else if( bc_type_name .eq. 'PRESSURE') then
        bc_type_tag = PRESSURE
        grid % bnd_cond % type(bc) = PRESSURE
      else
        if(this_proc < 2)  &
          print *, '# ERROR!  Load_Boundary_Conditions: '//        &
                   '# Unknown boundary condition type: ',  &
                   bc_type_name
        call Comm_Mod_End
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
        do c = -1, -grid % n_bnd_cells, -1
          if(grid % bnd_cond % color(c) .eq. bc) then

            ! Temperature
            if(heat_transfer) then
              i = Key_Ind('T', keys, nks)
              if(i > 0) t % bnd_cond_type(c) = bc_type_tag
              i = Key_Ind('Q', keys, nks)
              if(i > 0) t % bnd_cond_type(c) = bc_type_tag
            end if

            ! Multiphase flow
            if (multiphase_model .eq. VOLUME_OF_FLUID) then
              i = Key_Ind('VOF', keys, nks)
              if(i > 0) vof % bnd_cond_type(c) = bc_type_tag
              i = Key_Ind('VOF_C_ANG', keys, nks)
              if(i > 0) vof % bnd_cond_type(c) = bc_type_tag
            end if

            ! For scalars
            do sc = 1, flow % n_scalars
              i = Key_Ind(scalar(sc) % name, keys, nks)
              if(i > 0) scalar(sc) % bnd_cond_type(c) = bc_type_tag
              i = Key_Ind(scalar(sc) % flux_name, keys, nks)
              if(i > 0) scalar(sc) % bnd_cond_type(c) = bc_type_tag
            end do

          end if  ! bnd_color .eq. bc

        end do

        ! Distribute b.c. values
        do c = -1, -grid % n_bnd_cells, -1
          if(grid % bnd_cond % color(c) .eq. bc) then

            ! For velocity and pressure
            i = Key_Ind('U', keys, nks); if(i > 0) u % b(c) = vals(i)
            i = Key_Ind('V', keys, nks); if(i > 0) v % b(c) = vals(i)
            i = Key_Ind('W', keys, nks); if(i > 0) w % b(c) = vals(i)
            i = Key_Ind('P', keys, nks); if(i > 0) p % b(c) = vals(i)

            ! Temperature
            if(heat_transfer) then
              i = Key_Ind('T', keys, nks)
              if(i > 0) t % b(c) = vals(i)
              i = Key_Ind('Q', keys, nks)
              if(i > 0) t % q(c) = vals(i)
            end if

            ! Multiphase flow
            if (multiphase_model .eq. VOLUME_OF_FLUID) then
              i = Key_Ind('VOF', keys, nks)
              if(i > 0) vof % b(c) = vals(i)
              i = Key_Ind('VOF_C_ANG', keys, nks)
              if(i > 0) vof % q(c) = vals(i)
            end if

            ! For scalars
            do sc = 1, flow % n_scalars
              i = Key_Ind(scalar(sc) % name, keys, nks)
              if(i > 0) scalar(sc) % b(c) = vals(i)
              i = Key_Ind(scalar(sc) % flux_name, keys, nks)
              if(i > 0) scalar(sc) % q(c) = vals(i)
            end do

            ! For turbulence models
            if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
               turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
              i = Key_Ind('UU',  keys, nks); if(i > 0) uu  % b(c) = vals(i)
              i = Key_Ind('VV',  keys, nks); if(i > 0) vv  % b(c) = vals(i)
              i = Key_Ind('WW',  keys, nks); if(i > 0) ww  % b(c) = vals(i)
              i = Key_Ind('UV',  keys, nks); if(i > 0) uv  % b(c) = vals(i)
              i = Key_Ind('UW',  keys, nks); if(i > 0) uw  % b(c) = vals(i)
              i = Key_Ind('VW',  keys, nks); if(i > 0) vw  % b(c) = vals(i)
              i = Key_Ind('EPS', keys, nks); if(i > 0) eps % b(c) = vals(i)

              if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
                i = Key_Ind('F22', keys, nks); if(i > 0) f22 % b(c) = vals(i)
              end if
            end if

            if(turbulence_model .eq. K_EPS) then
              i = Key_Ind('KIN', keys, nks); if(i > 0) kin % b(c) = vals(i)
              i = Key_Ind('EPS', keys, nks); if(i > 0) eps % b(c) = vals(i)
              turb % y_plus(c) = 1.1
              if(heat_transfer) then
                i = Key_Ind('T2',  keys, nks); if(i > 0) t2 % b(c) = vals(i)
              end if
            end if

            if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
               turbulence_model .eq. HYBRID_LES_RANS) then
              i = Key_Ind('KIN',  keys, nks); if(i > 0) kin  % b(c) = vals(i)
              i = Key_Ind('EPS',  keys, nks); if(i > 0) eps  % b(c) = vals(i)
              i = Key_Ind('ZETA', keys, nks); if(i > 0) zeta % b(c) = vals(i)
              i = Key_Ind('F22',  keys, nks); if(i > 0) f22  % b(c) = vals(i)
              if(heat_transfer) then
                i = Key_Ind('T2',  keys, nks); if(i > 0) t2  % b(c) = vals(i)
              end if
            end if

            if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
               turbulence_model .eq. DES_SPALART) then
              i = Key_Ind('VIS',  keys, nks); if(i > 0) vis % b(c) = vals(i)
            end if
          end if

        end do

      !---------------------------------------------!
      !                                             !
      !   Boundary values are specified in a file   !
      !                                             !
      !---------------------------------------------!
      else

        call Control_Mod_Read_Strings_On('FILE', name_prof, nvs, .false.)

        call File_Mod_Open_File_For_Reading(name_prof(1), fu)
        call File_Mod_Read_Line(fu)
        read(line % tokens(1),*) n_points  ! number of points

        allocate(prof(n_points, 0:nks))

        !----------------------------------!
        !   Read the entire profile file   !
        !----------------------------------!
        do m = 1, n_points
          call File_Mod_Read_Line(fu)
          do i = 1, nks
            read(line % tokens(i), *) prof(m,i)
          end do
        end do
        close(fu)

        !------------------------!
        !   A plane is defined   !
        !------------------------!
        if(keys(1) .eq. 'X' .and. keys(2) .eq. 'Y' .or.  &
           keys(1) .eq. 'X' .and. keys(2) .eq. 'Z' .or.  &
           keys(1) .eq. 'Y' .and. keys(2) .eq. 'Z') then

          ! Set the closest point
          do c = -1, -grid % n_bnd_cells, -1

            ! Distribute b.c. types
            if(grid % bnd_cond % color(c) .eq. bc) then

              ! For temperature
              if(heat_transfer) then
                i = Key_Ind('T', keys, nks)
                if(i > 0) t % bnd_cond_type(c) = bc_type_tag
                i = Key_Ind('Q', keys, nks)
                if(i > 0) t % bnd_cond_type(c) = bc_type_tag
              end if

              ! For scalars
              do sc = 1, flow % n_scalars
                i = Key_Ind(scalar(sc) % name, keys, nks)
                if(i > 0) scalar(sc) % bnd_cond_type(c) = bc_type_tag
                i = Key_Ind(scalar(sc) % flux_name, keys, nks)
                if(i > 0) scalar(sc) % bnd_cond_type(c) = bc_type_tag
              end do

            end if

            ! Distribute b.c. values
            if(grid % bnd_cond % color(c) .eq. bc) then

              dist_min = HUGE
              do m = 1, n_points

                i = Key_Ind('X', keys, nks); prof(m,0) = 0.0;  x = prof(m,i)
                i = Key_Ind('Y', keys, nks); prof(m,0) = 0.0;  y = prof(m,i)
                i = Key_Ind('Z', keys, nks); prof(m,0) = 0.0;  z = prof(m,i)

                if(keys(1) .eq. 'Y' .and. keys(2) .eq. 'Z') then
                  dist = Math_Mod_Distance(                         &
                                  y,            z,            0.0,  &
                                  grid % yc(c), grid % zc(c), 0.0)

                else if(keys(1) .eq. 'X' .and. keys(2) .eq. 'Z') then
                  dist = Math_Mod_Distance(                         &
                                  x,            z,            0.0,  &
                                  grid % xc(c), grid % zc(c), 0.0)

                else if(keys(1) .eq. 'X' .and. keys(2) .eq. 'Y') then
                  dist = Math_Mod_Distance(                         &
                                  x,            y,            0.0,  &
                                  grid % xc(c), grid % yc(c), 0.0)

                end if

                ! Store closest point in k
                if(dist < dist_min) then
                  dist_min = dist
                  k = m
                end if

              end do

              ! For velocity and pressure
              i = Key_Ind('U', keys, nks); if(i > 0) u % b(c) = prof(k,i)
              i = Key_Ind('V', keys, nks); if(i > 0) v % b(c) = prof(k,i)
              i = Key_Ind('W', keys, nks); if(i > 0) w % b(c) = prof(k,i)
              i = Key_Ind('P', keys, nks); if(i > 0) p % b(c) = prof(k,i)

              ! For temperature
              if(heat_transfer) then
                i = Key_Ind('T', keys, nks)
                if(i > 0) t % b(c) = prof(k,i)
                i = Key_Ind('Q', keys, nks)
                if(i > 0) t % q(c) = prof(k,i)
              end if

              ! For scalars
              do sc = 1, flow % n_scalars
                i = Key_Ind(scalar(sc) % name, keys, nks)
                if(i > 0) scalar(sc) % b(c) = prof(k,i)
                i = Key_Ind(scalar(sc) % flux_name, keys, nks)
                if(i > 0) scalar(sc) % q(c) = prof(k,i)
              end do

              ! For turbulence models
              if(turbulence_model .eq. K_EPS) then
                i = Key_Ind('KIN', keys, nks); if(i > 0) kin % b(c) = prof(k,i)
                i = Key_Ind('EPS', keys, nks); if(i > 0) eps % b(c) = prof(k,i)
                if(heat_transfer) then
                  i = Key_Ind('T2',  keys, nks); if(i>0) t2  % b(c) = prof(k,i)
                end if
              end if

              if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
                 turbulence_model .eq. HYBRID_LES_RANS) then
                i = Key_Ind('KIN',  keys, nks); if(i>0) kin  % b(c) = prof(k,i)
                i = Key_Ind('EPS',  keys, nks); if(i>0) eps  % b(c) = prof(k,i)
                i = Key_Ind('ZETA', keys, nks); if(i>0) zeta % b(c) = prof(k,i)
                i = Key_Ind('F22',  keys, nks); if(i>0) f22  % b(c) = prof(k,i)
                if(heat_transfer) then
                  i = Key_Ind('T2',  keys, nks); if(i>0) t2  % b(c) = prof(k,i)
                end if
              end if

              if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
                 turbulence_model .eq. DES_SPALART) then
                i = Key_Ind('VIS', keys, nks); if(i > 0) vis % b(c) = prof(k,i)
              end if

              if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
                 turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
                i = Key_Ind('UU', keys, nks); if(i > 0) uu  % b(c) = prof(k,i)
                i = Key_Ind('VV', keys, nks); if(i > 0) vv  % b(c) = prof(k,i)
                i = Key_Ind('WW', keys, nks); if(i > 0) ww  % b(c) = prof(k,i)
                i = Key_Ind('UV', keys, nks); if(i > 0) uv  % b(c) = prof(k,i)
                i = Key_Ind('UW', keys, nks); if(i > 0) uw  % b(c) = prof(k,i)
                i = Key_Ind('VW', keys, nks); if(i > 0) vw  % b(c) = prof(k,i)
                i = Key_Ind('EPS',keys, nks); if(i > 0) eps % b(c) = prof(k,i)

                if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
                  i = Key_Ind('F22', keys, nks); if(i>0) f22 % b(c) = prof(k,i)
                end if
              end if
            end if      !end if(grid % bnd_cond % color(c) .eq. n)
          end do        !end do c = -1, -grid % n_bnd_cells, -1

        !----------------------------!
        !   A plane is not defined   !
        !----------------------------!
        else  ! dir .eq. "XPL" ...

          do c = -1, -grid % n_bnd_cells, -1

            if(grid % bnd_cond % color(c) .eq. bc) then

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

                ! Wall distance too
                else if( (keys(1) .eq. 'WD'           .and.  &
                     grid % wall_dist(c) >= min(x,xp) .and.  &
                     grid % wall_dist(c) <= max(x,xp)) ) then
                  wi = ( max(x,xp) - grid % wall_dist(c) )   &
                     / ( max(x,xp) - min(x,xp) )
                  here = .true.
                end if

                if(here) then

                  ! For temperature
                  if(heat_transfer) then
                    i = Key_Ind('T',keys,nks)
                    if(i > 0) t % bnd_cond_type(c) = bc_type_tag
                    i = Key_Ind('Q',keys,nks)
                    if(i > 0) t % bnd_cond_type(c) = bc_type_tag
                  end if

                  ! For scalars
                  do sc = 1, flow % n_scalars
                    i = Key_Ind(scalar(sc) % name, keys, nks)
                    if(i > 0) scalar(sc) % bnd_cond_type(c) = bc_type_tag
                    i = Key_Ind(scalar(sc) % flux_name, keys, nks)
                    if(i > 0) scalar(sc) % bnd_cond_type(c) = bc_type_tag
                  end do

                end if  ! here
              end do    ! m, points
            end if      ! bnd_color .eq. bc

            if(grid % bnd_cond % color(c) .eq. bc) then

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

                ! Wall distance too
                else if( (keys(1) .eq. 'WD'           .and.  &
                     grid % wall_dist(c) >= min(x,xp) .and.  &
                     grid % wall_dist(c) <= max(x,xp)) ) then
                  wi = ( max(x,xp) - grid % wall_dist(c) )   &
                     / ( max(x,xp) - min(x,xp) )
                  here = .true.
                end if

                ! Interpolate the profiles
                if(here) then

                  ! For velocity and pressure
                  i = Key_Ind('U',keys,nks)
                  if(i > 0) u % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)
                  i = Key_Ind('V',keys,nks)
                  if(i > 0) v % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)
                  i = Key_Ind('W',keys,nks)
                  if(i > 0) w % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)
                  i = Key_Ind('P',keys,nks)
                  if(i > 0) p % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)

                  ! For temperature
                  if(heat_transfer) then
                    i = Key_Ind('T',keys,nks)
                    if(i > 0) t % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)
                    if(i > 0) t % bnd_cond_type(c) = bc_type_tag
                    i = Key_Ind('Q',keys,nks)
                    if(i > 0) t % q(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)
                    if(i > 0) t % bnd_cond_type(c) = bc_type_tag
                  end if

                  ! For scalars
                  do sc = 1, flow % n_scalars
                    i = Key_Ind(scalar(sc) % name, keys, nks)
                    if(i > 0) &
                      scalar(sc) % b(c)=wi*prof(m,i)+(1.-wi)*prof(m+1,i)
                    i = Key_Ind(scalar(sc) % flux_name, keys, nks)
                    if(i > 0) &
                      scalar(sc) % q(c)=wi*prof(m,i)+(1.-wi)*prof(m+1,i)
                  end do

                  ! For turbulence models
                  if(turbulence_model .eq. K_EPS) then

                    i = Key_Ind('KIN',keys,nks)
                    if(i > 0) kin % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)

                    i = Key_Ind('EPS',keys,nks)
                    if(i > 0) eps % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)

                    if(heat_transfer) then
                      i = Key_Ind('T2',keys,nks)
                      if(i > 0) t2 % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)
                    end if

                  end if

                  if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
                     turbulence_model .eq. HYBRID_LES_RANS) then

                    i = Key_Ind('KIN',keys,nks)
                    if(i > 0) kin % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)

                    i = Key_Ind('EPS',keys,nks)
                    if(i > 0) eps % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)

                    i = Key_Ind('ZETA',keys,nks)
                    if(i > 0) zeta % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)

                    i = Key_Ind('F22',keys,nks)
                    if(i > 0) f22 % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)

                    if(heat_transfer) then
                      i = Key_Ind('T2',keys,nks)
                      if(i > 0) t2 % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)
                    end if

                  end if

                  if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
                     turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then

                    i = Key_Ind('UU', keys, nks)
                    if(i > 0) uu % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)

                    i = Key_Ind('VV', keys, nks)
                    if(i > 0) vv % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)

                    i = Key_Ind('WW', keys, nks)
                    if(i > 0) ww % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)

                    i = Key_Ind('UV', keys, nks)
                    if(i > 0) uv % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)

                    i = Key_Ind('UW', keys, nks)
                    if(i > 0) uw % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)

                    i = Key_Ind('VW', keys, nks)
                    if(i > 0) vw % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)

                    i = Key_Ind('EPS', keys, nks)
                    if(i > 0) eps % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)

                    if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
                      i = Key_Ind('F22', keys, nks)
                      if(i > 0)f22 % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)
                    end if
                  end if

                  if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
                     turbulence_model .eq. DES_SPALART) then
                    i = Key_Ind('VIS',keys,nks)
                    if(i > 0) vis % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)
                  end if

                end if  ! (here)
              end do  ! m = 1, n_points-1
            end if
          end do  ! c = -1, -grid % n_bnd_cells, -1
        end if  ! plane is defined?
        close(fu)
      end if  ! boundary defined in a file
    end do

  end do

  !-----------------------------------!
  !                                   !
  !   Read data on synthetic eddies   !
  !                                   !
  !-----------------------------------!
  turb_planes % n_planes = 0
  do bc = 1, grid % n_bnd_cond  ! imagine there are as many eddies as bcs
    call Control_Mod_Position_At_Two_Keys('SYNTHETIC_EDDIES',          &
                                          grid % bnd_cond % name(bc),  &
                                          found,                       &
                                          .false.)
    if(found) then
      turb_planes % n_planes = turb_planes % n_planes + 1
      call Control_Mod_Read_Int_Item_On ('NUMBER_OF_EDDIES', 24, edd_n, .false.)
      call Control_Mod_Read_Real_Item_On('MAX_EDDY_RADIUS',  .2, edd_r, .false.)
      call Control_Mod_Read_Real_Item_On('EDDY_INTENSITY',   .1, edd_i, .false.)
      call Eddies_Mod_Allocate(turb_planes % plane(turb_planes % n_planes),  &
                               edd_n,                                        &
                               edd_r,                                        &
                               edd_i,                                        &
                               flow,                                         &
                               grid % bnd_cond % name(bc))
    end if
  end do
  if(turb_planes % n_planes > 0 .and. this_proc < 2) then
    print *, '# Found ', turb_planes % n_planes, ' turbulent planes'
  end if

  !---------------------------------------!
  !                                       !
  !                                       !
  !   Copy all "b" values to "n" values   !
  !                                       !
  !                                       !
  !---------------------------------------!
! if(.not. backup) then
    do c = -1, -grid % n_bnd_cells, -1

      u % n(c) = u % b(c)
      v % n(c) = v % b(c)
      w % n(c) = w % b(c)
      p % n(c) = p % b(c)

      if(heat_transfer) then
        t % n(c) = t % b(c)
      end if

      if (multiphase_model .eq. VOLUME_OF_FLUID) then
        vof % n(c) = vof % b(c)
      end if

      do sc = 1, flow % n_scalars
        scalar(sc) % n(c) = scalar(sc) % b(c)
      end do

      if(turbulence_model .eq. RSM_MANCEAU_HANJALIC .or.  &
         turbulence_model .eq. RSM_HANJALIC_JAKIRLIC) then
        uu  % n(c) = uu  % b(c)
        vv  % n(c) = vv  % b(c)
        ww  % n(c) = ww  % b(c)
        uv  % n(c) = uv  % b(c)
        uw  % n(c) = uw  % b(c)
        vw  % n(c) = vw  % b(c)
        eps % n(c) = eps % b(c)

        if(turbulence_model .eq. RSM_MANCEAU_HANJALIC) then
          f22 % n(c) = f22 % b(c)
        end if
      end if

      if(turbulence_model .eq. K_EPS) then
        kin % n(c) = kin % b(c)
        eps % n(c) = eps % b(c)
        if(heat_transfer) then
          t2 % n(c) = t2 % b(c)
        end if
      end if

      if(turbulence_model .eq. K_EPS_ZETA_F .or.  &
         turbulence_model .eq. HYBRID_LES_RANS) then
        kin  % n(c) = kin  % b(c)
        eps  % n(c) = eps  % b(c)
        zeta % n(c) = zeta % b(c)
        f22  % n(c) = eps  % b(c)
        if(heat_transfer) then
          t2 % n(c) = t2 % b(c)
        end if
      end if

      if(turbulence_model .eq. SPALART_ALLMARAS .or.  &
         turbulence_model .eq. DES_SPALART) then
        vis % n(c) = vis % b(c)
      end if

    end do  ! through boundary cells
! end if ! backup

  !------------------------------!
  !   Find the near-wall cells   !
  !------------------------------!
  grid % cell_near_wall = .false.

  do s = 1, grid % n_faces
    c1 = grid % faces_c(1,s)
    c2 = grid % faces_c(2,s)

    if(c2 < 0) then
      if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALL .or.  &
         Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
        grid % cell_near_wall(c1) = .true.
      end if
    end if

  end do  ! faces

  end subroutine
