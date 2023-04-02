!==============================================================================!
  subroutine Boundary_Conditions(Rc, Flow, Turb, Vof, turb_planes)
!------------------------------------------------------------------------------!
!   Reads boundary condition from control file                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Read_Controls_Type)         :: Rc
  type(Field_Type),          target :: Flow
  type(Turb_Type),           target :: Turb
  type(Vof_Type),            target :: Vof
  type(Turb_Plane_Type)             :: turb_planes
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Var_Type),  pointer :: u, v, w, t, p, fun
  type(Var_Type),  pointer :: kin, eps, f22, zeta, vis, t2
  type(Var_Type),  pointer :: uu, vv, ww, uv, uw, vw
  type(Var_Type),  pointer :: scalar(:)
  real,            pointer :: z_o(:)
  integer                  :: c,m,l,k,i,bc,n_points,nks,nvs,sc,c1,c2,s,fu
  character(SL)            :: name_prof(128)
  real                     :: wi, dist_min, x, y, z, xp, dist
  real, allocatable        :: prof(:,:)
  logical                  :: here
  character(SL)            :: bc_type_name, try_str
  integer                  :: bc_type_tag
  character(SL)            :: keys(128)
  real                     :: vals(0:128)         ! they start from zero!
  integer                  :: types_per_reg(128)  ! how many types in a region
  character(SL)            :: types_names(128)    ! name of each type
  logical                  :: types_file(128)     ! type specified in a file?
  integer                  :: c_types             ! counter types
  integer                  :: edd_n
  real                     :: edd_r
  real                     :: edd_i
  logical                  :: found
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Rc)
!==============================================================================!

  ! Take aliases
  Grid   => Flow % pnt_grid
  t      => Flow % t
  p      => Flow % p
  scalar => Flow % scalar
  vis    => Turb % vis
  t2     => Turb % t2
  z_o    => Turb % z_o
  fun    => Vof % fun

  call Flow % Alias_Momentum    (u, v, w)
  call Turb % Alias_K_Eps_Zeta_F(kin, eps, zeta, f22)
  call Turb % Alias_Stresses    (uu, vv, ww, uv, uw, vw)

  !-------------------------------------------------------------------!
  !   Read wall roughness, if specified as a constant for all walls   !
  !   (If it is not specified in the control file, it will be zero)   !
  !-------------------------------------------------------------------!
  call Control_Mod_Roughness_Coefficient(Turb % z_o)

  !----------------------------------------------------------------!
  !   Count number of types per boundary condition, total number   !
  !        of types specified, and also extract their names        !
  !----------------------------------------------------------------!
  types_per_reg(:) = 0
  types_file(:)      = .false.
  c_types            = 0

  do bc = Boundary_Regions()
    call Control % Position_At_Two_Keys('BOUNDARY_CONDITION',      &
                                        Grid % region % name(bc),  &
                                        found,                     &
                                        .false.)
    if(found) then
1     continue

      ! Try to read next 'TYPE' in the control file
      call Control % Read_Char_Item_On('TYPE', 'VOID', bc_type_name, .false.)

      ! Get out of the loop if you fail
      if(bc_type_name .eq. 'VOID') goto 2

      ! Skip following two lines
      call Control % Read_Char_Item_On('VARIABLES', 'VOID', try_str, .false.)
      call Control % Read_Char_Item_On('VALUES',    'VOID', try_str, .false.)

      types_per_reg(bc) = types_per_reg(bc) + 1
      c_types = c_types + 1
      types_names(c_types) = bc_type_name

      ! If try_str is 'VOID', it didn't find 'VALUES'
      ! meaning that the keyword 'FILE' was specified
      if(try_str .eq. 'VOID') then
        types_file(c_types) = .true.
      end if

      goto 1
    else
      call Message % Error(72,                                             &
             'Boundary conditions for '//trim(Grid % region % name(bc))//  &
             ' not specified in the control file!  \n \n  Exiting!',       &
             file=__FILE__, line=__LINE__)
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

  do bc = Boundary_Regions()

    ! Position yourself well
    call Control % Position_At_Two_Keys('BOUNDARY_CONDITION',      &
                                        Grid % region % name(bc),  &
                                        found,                     &
                                        .false.)
    do l = 1, types_per_reg(bc)

      ! Update the counter
      c_types = c_types + 1

      !---------------------------------------------!
      !                                             !
      !   Read first line which is common for all   !
      !                                             !
      !---------------------------------------------!
      call Control % Read_Char_Item_On('TYPE', 'WALL', bc_type_name, .false.)
      call String % To_Upper_Case(bc_type_name)

      ! Copy boundary conditions which were given for the Grid
      if( bc_type_name .eq. 'INFLOW') then
        bc_type_tag = INFLOW
        Grid % region % type(bc) = INFLOW
      else if( bc_type_name .eq. 'WALL') then
        bc_type_tag = WALL
        Grid % region % type(bc) = WALL
      else if( bc_type_name .eq. 'OUTFLOW') then
        bc_type_tag = OUTFLOW
        Grid % region % type(bc) = OUTFLOW
      else if( bc_type_name .eq. 'SYMMETRY') then
        bc_type_tag = SYMMETRY
        Grid % region % type(bc) = SYMMETRY
      else if( bc_type_name .eq. 'WALL_FLUX') then
        bc_type_tag = WALLFL
        Grid % region % type(bc) = WALLFL
      else if( bc_type_name .eq. 'CONVECTIVE') then
        bc_type_tag = CONVECT
        Grid % region % type(bc) = CONVECT
      else if( bc_type_name .eq. 'PRESSURE') then
        bc_type_tag = PRESSURE
        Grid % region % type(bc) = PRESSURE
      else
        call Message % Error(72,                                            &
                 'Unknown boundary condition type: '//trim(bc_type_name)//  &
                 '. \n \n This error is critical.  Exiting!',               &
                 file=__FILE__, line=__LINE__)
      end if

      !----------------------------------------------!
      !                                              !
      !   Read second line which is common for all   !
      !                                              !
      !----------------------------------------------!
      call Control % Read_Strings_On('VARIABLES', keys, nks, .false.)
      do i = 1, nks
        call String % To_Upper_Case(keys(i))
      end do

      !-----------------------------------------------------------------!
      !                                                                 !
      !   Boundary values are specified in a list (and not in a file)   !
      !                                                                 !
      !-----------------------------------------------------------------!
      if( .not. types_file(c_types) ) then
        call Control % Read_Real_Vector_On('VALUES', vals(1), nvs, .false.)

        !--------------------------------------------------!
        !   Distribute boundary values to boundary cells   !
        !--------------------------------------------------!

        ! Distribute b.c. tags only.
        do c = Cells_In_Region(bc)

          ! Temperature
          if(Flow % heat_transfer) then
            i = Key_Ind('T', keys, nks)
            if(i > 0) then
              t % bnd_cond_type(c) = bc_type_tag
              if(bc_type_tag .eq. WALLFL) t    % bnd_cond_type(c)    = WALL
              if(bc_type_tag .eq. WALLFL) Grid % region % type(bc) = WALL
            end if
            i = Key_Ind('Q', keys, nks)
            if(i > 0) then
              t % bnd_cond_type(c) = bc_type_tag
              if(bc_type_tag .eq. WALL) t    % bnd_cond_type(c)    = WALLFL
              if(bc_type_tag .eq. WALL) Grid % region % type(bc) = WALLFL
            end if
          end if

          ! Volume of fluid -> still to be worked around
          if (Flow % with_interface) then
            i = Key_Ind('VOF', keys, nks)
            if(i > 0) fun % bnd_cond_type(c) = bc_type_tag
            i = Key_Ind('VOF_C_ANG', keys, nks)
            if(i > 0) fun % bnd_cond_type(c) = bc_type_tag
          end if

          ! For scalars
          do sc = 1, Flow % n_scalars
            i = Key_Ind(scalar(sc) % name, keys, nks)
            if(i > 0) then
              scalar(sc) % bnd_cond_type(c) = bc_type_tag
            end if
            i = Key_Ind(scalar(sc) % flux_name, keys, nks)
            if(i > 0) then
              scalar(sc) % bnd_cond_type(c) = WALLFL
            end if
          end do

        end do    ! Cells_In_Region

        ! Distribute b.c. values
        do c = Cells_In_Region(bc)

          ! For velocity, pressure and wall roughness
          i = Key_Ind('U',   keys, nks); if(i > 0) u % b(c) = vals(i)
          i = Key_Ind('V',   keys, nks); if(i > 0) v % b(c) = vals(i)
          i = Key_Ind('W',   keys, nks); if(i > 0) w % b(c) = vals(i)
          i = Key_Ind('P',   keys, nks); if(i > 0) p % b(c) = vals(i)
          i = Key_Ind('Z_O', keys, nks); if(i > 0) z_o  (c) = vals(i)

          ! Temperature
          if(Flow % heat_transfer) then
            i = Key_Ind('T', keys, nks)
            if(i > 0) t % b(c) = vals(i)
            i = Key_Ind('Q', keys, nks)
            if(i > 0) t % q(c) = vals(i)
          end if

          ! Multiphase Flow
          if (Flow % with_interface) then
            i = Key_Ind('VOF', keys, nks)
            if(i > 0) fun % b(c) = vals(i)
            i = Key_Ind('VOF_C_ANG', keys, nks)
            if(i > 0) fun % q(c) = vals(i)
          end if

          ! For scalars
          do sc = 1, Flow % n_scalars
            i = Key_Ind(scalar(sc) % name, keys, nks)
            if(i > 0) scalar(sc) % b(c) = vals(i)
            i = Key_Ind(scalar(sc) % flux_name, keys, nks)
            if(i > 0) scalar(sc) % q(c) = vals(i)
          end do

          ! For turbulence models
          if(Turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
             Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
            i = Key_Ind('UU',  keys, nks); if(i > 0) uu  % b(c) = vals(i)
            i = Key_Ind('VV',  keys, nks); if(i > 0) vv  % b(c) = vals(i)
            i = Key_Ind('WW',  keys, nks); if(i > 0) ww  % b(c) = vals(i)
            i = Key_Ind('UV',  keys, nks); if(i > 0) uv  % b(c) = vals(i)
            i = Key_Ind('UW',  keys, nks); if(i > 0) uw  % b(c) = vals(i)
            i = Key_Ind('VW',  keys, nks); if(i > 0) vw  % b(c) = vals(i)
            i = Key_Ind('EPS', keys, nks); if(i > 0) eps % b(c) = vals(i)

            if(Turb % model .eq. RSM_MANCEAU_HANJALIC) then
              i = Key_Ind('F22', keys, nks); if(i > 0) f22 % b(c) = vals(i)
            end if
          end if

          if(Turb % model .eq. K_EPS) then
            i = Key_Ind('KIN', keys, nks); if(i > 0) kin % b(c) = vals(i)
            i = Key_Ind('EPS', keys, nks); if(i > 0) eps % b(c) = vals(i)
            Turb % y_plus(c) = 1.1
            if(Flow % heat_transfer) then
              i = Key_Ind('T2',  keys, nks); if(i > 0) t2 % b(c) = vals(i)
            end if
          end if

          if(Turb % model .eq. K_EPS_ZETA_F .or.  &
             Turb % model .eq. HYBRID_LES_RANS) then
            i = Key_Ind('KIN',  keys, nks); if(i > 0) kin  % b(c) = vals(i)
            i = Key_Ind('EPS',  keys, nks); if(i > 0) eps  % b(c) = vals(i)
            i = Key_Ind('ZETA', keys, nks); if(i > 0) zeta % b(c) = vals(i)
            i = Key_Ind('F22',  keys, nks); if(i > 0) f22  % b(c) = vals(i)
            if(Flow % heat_transfer) then
              i = Key_Ind('T2',  keys, nks); if(i > 0) t2  % b(c) = vals(i)
            end if
          end if

          if(Turb % model .eq. SPALART_ALLMARAS .or.  &
             Turb % model .eq. DES_SPALART) then
            i = Key_Ind('VIS',  keys, nks); if(i > 0) vis % b(c) = vals(i)
          end if

        end do    ! Cells_In_Region

      !---------------------------------------------!
      !                                             !
      !   Boundary values are specified in a file   !
      !                                             !
      !---------------------------------------------!
      else  !  types_file(c_types) == .true.

        call Control % Read_Strings_On('FILE', name_prof, nvs, .false.)

        call File % Open_For_Reading_Ascii(name_prof(1), fu)
        call File % Read_Line(fu)
        read(Line % tokens(1),*) n_points  ! number of points

        !---------------------------------!
        !   Allocate memory for profile   !
        !---------------------------------!
        allocate(prof(n_points, 0:nks))

        !----------------------------------!
        !   Read the entire profile file   !
        !----------------------------------!
        do m = 1, n_points
          call File % Read_Line(fu)
          do i = 1, nks
            read(Line % tokens(i), *) prof(m,i)
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
          do c = Cells_In_Region(bc)

            ! Distribute b.c. types

            ! For temperature
            if(Flow % heat_transfer) then
              i = Key_Ind('T', keys, nks)
              if(i > 0) then
                t % bnd_cond_type(c) = bc_type_tag
                if(bc_type_tag .eq. WALLFL) t    % bnd_cond_type(c)    = WALL
                if(bc_type_tag .eq. WALLFL) Grid % region % type(bc) = WALL
              end if
              i = Key_Ind('Q', keys, nks)
              if(i > 0) then
                t % bnd_cond_type(c) = bc_type_tag
                if(bc_type_tag .eq. WALL) t    % bnd_cond_type(c)    = WALLFL
                if(bc_type_tag .eq. WALL) Grid % region % type(bc) = WALLFL
              end if
            end if

            ! For scalars
            do sc = 1, Flow % n_scalars
              i = Key_Ind(scalar(sc) % name, keys, nks)
              if(i > 0) then
                scalar(sc) % bnd_cond_type(c) = bc_type_tag
              end if
              i = Key_Ind(scalar(sc) % flux_name, keys, nks)
              if(i > 0) then
                scalar(sc) % bnd_cond_type(c) = WALLFL
              end if
            end do

            ! Distribute b.c. values
            dist_min = HUGE
            do m = 1, n_points

              i = Key_Ind('X', keys, nks); prof(m,0) = 0.0;  x = prof(m,i)
              i = Key_Ind('Y', keys, nks); prof(m,0) = 0.0;  y = prof(m,i)
              i = Key_Ind('Z', keys, nks); prof(m,0) = 0.0;  z = prof(m,i)

              if(keys(1) .eq. 'Y' .and. keys(2) .eq. 'Z') then
                dist = Math % Distance(                           &
                                y,            z,            0.0,  &
                                Grid % yc(c), Grid % zc(c), 0.0)

              else if(keys(1) .eq. 'X' .and. keys(2) .eq. 'Z') then
                dist = Math % Distance(                           &
                                x,            z,            0.0,  &
                                Grid % xc(c), Grid % zc(c), 0.0)

              else if(keys(1) .eq. 'X' .and. keys(2) .eq. 'Y') then
                dist = Math % Distance(                           &
                                x,            y,            0.0,  &
                                Grid % xc(c), Grid % yc(c), 0.0)

              end if

              ! Store closest point in k
              if(dist < dist_min) then
                dist_min = dist
                k = m
              end if

            end do

            ! For velocity, pressure and wall roughness
            i = Key_Ind('U',   keys, nks); if(i > 0) u % b(c) = prof(k,i)
            i = Key_Ind('V',   keys, nks); if(i > 0) v % b(c) = prof(k,i)
            i = Key_Ind('W',   keys, nks); if(i > 0) w % b(c) = prof(k,i)
            i = Key_Ind('P',   keys, nks); if(i > 0) p % b(c) = prof(k,i)
            i = Key_Ind('Z_O', keys, nks); if(i > 0) z_o  (c) = prof(k,i)

            ! For temperature
            if(Flow % heat_transfer) then
              i = Key_Ind('T', keys, nks)
              if(i > 0) t % b(c) = prof(k,i)
              i = Key_Ind('Q', keys, nks)
              if(i > 0) t % q(c) = prof(k,i)
            end if

            ! For scalars
            do sc = 1, Flow % n_scalars
              i = Key_Ind(scalar(sc) % name, keys, nks)
              if(i > 0) scalar(sc) % b(c) = prof(k,i)
              i = Key_Ind(scalar(sc) % flux_name, keys, nks)
              if(i > 0) scalar(sc) % q(c) = prof(k,i)
            end do

            ! For turbulence models
            if(Turb % model .eq. K_EPS) then
              i = Key_Ind('KIN', keys, nks); if(i > 0) kin % b(c) = prof(k,i)
              i = Key_Ind('EPS', keys, nks); if(i > 0) eps % b(c) = prof(k,i)
              if(Flow % heat_transfer) then
                i = Key_Ind('T2',  keys, nks); if(i>0) t2  % b(c) = prof(k,i)
              end if
            end if

            if(Turb % model .eq. K_EPS_ZETA_F .or.  &
               Turb % model .eq. HYBRID_LES_RANS) then
              i = Key_Ind('KIN',  keys, nks); if(i>0) kin  % b(c) = prof(k,i)
              i = Key_Ind('EPS',  keys, nks); if(i>0) eps  % b(c) = prof(k,i)
              i = Key_Ind('ZETA', keys, nks); if(i>0) zeta % b(c) = prof(k,i)
              i = Key_Ind('F22',  keys, nks); if(i>0) f22  % b(c) = prof(k,i)
              if(Flow % heat_transfer) then
                i = Key_Ind('T2',  keys, nks); if(i>0) t2  % b(c) = prof(k,i)
              end if
            end if

            if(Turb % model .eq. SPALART_ALLMARAS .or.  &
               Turb % model .eq. DES_SPALART) then
              i = Key_Ind('VIS', keys, nks); if(i > 0) vis % b(c) = prof(k,i)
            end if

            if(Turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
               Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
              i = Key_Ind('UU', keys, nks); if(i > 0) uu  % b(c) = prof(k,i)
              i = Key_Ind('VV', keys, nks); if(i > 0) vv  % b(c) = prof(k,i)
              i = Key_Ind('WW', keys, nks); if(i > 0) ww  % b(c) = prof(k,i)
              i = Key_Ind('UV', keys, nks); if(i > 0) uv  % b(c) = prof(k,i)
              i = Key_Ind('UW', keys, nks); if(i > 0) uw  % b(c) = prof(k,i)
              i = Key_Ind('VW', keys, nks); if(i > 0) vw  % b(c) = prof(k,i)
              i = Key_Ind('EPS',keys, nks); if(i > 0) eps % b(c) = prof(k,i)

              if(Turb % model .eq. RSM_MANCEAU_HANJALIC) then
                i = Key_Ind('F22', keys, nks); if(i>0) f22 % b(c) = prof(k,i)
              end if
            end if

          end do        ! Cells_In_Region

        !----------------------------!
        !   A plane is not defined   !
        !----------------------------!
        else  ! dir .eq. "XPL" ...

          do c = Cells_In_Region(bc)

            do m = 1, n_points-1
              here = .false.

              i = Key_Ind(keys(1), keys, nks)
              prof(m,   0) = 0.0;
              prof(m+1, 0) = 0.0;
              x  = prof(m,i)
              xp = prof(m+1,i)

              ! Compute the weight factors
              if( keys(1) .eq. 'X' .and.  &
                  Grid % xc(c) >= x .and. Grid % xc(c) <= xp ) then
                wi = ( xp - Grid % xc(c) ) / (xp - x)
                here = .true.
              else if( keys(1) .eq. 'Y' .and.  &
                       Grid % yc(c) >= x .and. Grid % yc(c) <= xp ) then
                wi = ( xp - Grid % yc(c) ) / (xp - x)
                here = .true.
              else if( keys(1) .eq. 'Z' .and.  &
                       Grid % zc(c) >= x .and. Grid % zc(c) <= xp ) then
                wi = ( xp - Grid % zc(c) ) / (xp - x)
                here = .true.

              ! Beware; for cylindrical coordinates you have "inversion"
              else if( (keys(1) .eq. 'RX' .and.  &
                   sqrt(Grid % yc(c)**2 + Grid % zc(c)**2) >= xp .and.       &
                   sqrt(Grid % yc(c)**2 + Grid % zc(c)**2) <= x) ) then
                wi = ( xp - sqrt(Grid % yc(c)**2 + Grid % zc(c)**2) ) / (xp-x)
                here = .true.
              else if( (keys(1) .eq. 'RY' .and.  &
                   sqrt(Grid % xc(c)**2 + Grid % zc(c)**2) >= xp .and.       &
                   sqrt(Grid % xc(c)**2 + Grid % zc(c)**2) <= x) ) then
                wi = ( xp - sqrt(Grid % xc(c)**2 + Grid % zc(c)**2) ) / (xp-x)
                here = .true.
              else if( (keys(1) .eq. 'RZ' .and.  &
                   sqrt(Grid % xc(c)**2 + Grid % yc(c)**2) >= xp .and.       &
                   sqrt(Grid % xc(c)**2 + Grid % yc(c)**2) <= x) ) then
                wi = ( xp - sqrt(Grid % xc(c)**2 + Grid % yc(c)**2) ) / (xp-x)
                here = .true.

              ! Wall distance too
              else if( (keys(1) .eq. 'WD'           .and.  &
                   Grid % wall_dist(c) >= min(x,xp) .and.  &
                   Grid % wall_dist(c) <= max(x,xp)) ) then
                wi = ( max(x,xp) - Grid % wall_dist(c) )   &
                   / ( max(x,xp) - min(x,xp) )
                here = .true.
              end if

              if(here) then

                ! For temperature
                if(Flow % heat_transfer) then
                  i = Key_Ind('T',keys,nks)
                  if(i > 0) then
                    t % bnd_cond_type(c) = bc_type_tag
                    if(bc_type_tag .eq. WALLFL) t % bnd_cond_type(c) = WALL
                    if(bc_type_tag .eq. WALLFL)  &
                      Grid % region % type(bc) = WALL
                  end if
                  i = Key_Ind('Q',keys,nks)
                  if(i > 0) then
                    t % bnd_cond_type(c) = bc_type_tag
                    if(bc_type_tag .eq. WALL) t % bnd_cond_type(c) = WALLFL
                    if(bc_type_tag .eq. WALL)  &
                      Grid % region % type(bc) = WALLFL
                  end if
                end if

                ! For scalars
                do sc = 1, Flow % n_scalars
                  i = Key_Ind(scalar(sc) % name, keys, nks)
                  if(i > 0) then
                    scalar(sc) % bnd_cond_type(c) = bc_type_tag
                  end if
                  i = Key_Ind(scalar(sc) % flux_name, keys, nks)
                  if(i > 0) then
                    scalar(sc) % bnd_cond_type(c) = WALLFL
                  end if
                end do

              end if  ! here
            end do    ! m, points

            do m = 1, n_points-1
              here = .false.

              i = Key_Ind(keys(1), keys, nks)
              prof(m,   0) = 0.0;
              prof(m+1, 0) = 0.0;
              x  = prof(m,i)
              xp = prof(m+1,i)

              ! Compute the weight factors
              if( keys(1) .eq. 'X' .and.  &
                  Grid % xc(c) >= x .and. Grid % xc(c) <= xp ) then
                wi = ( xp - Grid % xc(c) ) / (xp - x)
                here = .true.
              else if( keys(1) .eq. 'Y' .and.  &
                       Grid % yc(c) >= x .and. Grid % yc(c) <= xp ) then
                wi = ( xp - Grid % yc(c) ) / (xp - x)
                here = .true.
              else if( keys(1) .eq. 'Z' .and.  &
                       Grid % zc(c) >= x .and. Grid % zc(c) <= xp ) then
                wi = ( xp - Grid % zc(c) ) / (xp - x)
                here = .true.

              ! Beware; for cylindrical coordinates you have "inversion"
              else if( (keys(1) .eq. 'RX' .and.  &
                   sqrt(Grid % yc(c)**2 + Grid % zc(c)**2) >= xp .and.       &
                   sqrt(Grid % yc(c)**2 + Grid % zc(c)**2) <= x) ) then
                wi = ( xp - sqrt(Grid % yc(c)**2 + Grid % zc(c)**2) ) / (xp-x)
                here = .true.
              else if( (keys(1) .eq. 'RY' .and.  &
                   sqrt(Grid % xc(c)**2 + Grid % zc(c)**2) >= xp .and.       &
                   sqrt(Grid % xc(c)**2 + Grid % zc(c)**2) <= x) ) then
                wi = ( xp - sqrt(Grid % xc(c)**2 + Grid % zc(c)**2) ) / (xp-x)
                here = .true.
              else if( (keys(1) .eq. 'RZ' .and.  &
                   sqrt(Grid % xc(c)**2 + Grid % yc(c)**2) >= xp .and.       &
                   sqrt(Grid % xc(c)**2 + Grid % yc(c)**2) <= x) ) then
                wi = ( xp - sqrt(Grid % xc(c)**2 + Grid % yc(c)**2) ) / (xp-x)
                here = .true.

              ! Wall distance too
              else if( (keys(1) .eq. 'WD'           .and.  &
                   Grid % wall_dist(c) >= min(x,xp) .and.  &
                   Grid % wall_dist(c) <= max(x,xp)) ) then
                wi = ( max(x,xp) - Grid % wall_dist(c) )   &
                   / ( max(x,xp) - min(x,xp) )
                here = .true.
              end if

              ! Interpolate the profiles
              if(here) then

                ! For velocity, pressure and wall roughness
                i = Key_Ind('U',keys,nks)
                if(i > 0) u % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)
                i = Key_Ind('V',keys,nks)
                if(i > 0) v % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)
                i = Key_Ind('W',keys,nks)
                if(i > 0) w % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)
                i = Key_Ind('P',keys,nks)
                if(i > 0) p % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)
                i = Key_Ind('Z_O',keys,nks)
                if(i > 0) z_o  (c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)

                ! For temperature
                if(Flow % heat_transfer) then
                  i = Key_Ind('T',keys,nks)
                  if(i > 0) t % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)
                  if(i > 0) then
                    t % bnd_cond_type(c) = bc_type_tag
                    if(bc_type_tag .eq. WALLFL) t % bnd_cond_type(c) = WALL
                    if(bc_type_tag .eq. WALLFL)  &
                      Grid % region % type(bc) = WALL
                  end if
                  i = Key_Ind('Q',keys,nks)
                  if(i > 0) t % q(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)
                  if(i > 0) then
                    t % bnd_cond_type(c) = bc_type_tag
                    if(bc_type_tag .eq. WALL) t % bnd_cond_type(c) = WALLFL
                    if(bc_type_tag .eq. WALL)  &
                      Grid % region % type(bc) = WALLFL
                  end if
                end if

                ! For scalars
                do sc = 1, Flow % n_scalars
                  i = Key_Ind(scalar(sc) % name, keys, nks)
                  if(i > 0) then
                    scalar(sc) % b(c)=wi*prof(m,i)+(1.-wi)*prof(m+1,i)
                    scalar(sc) % bnd_cond_type(c) = bc_type_tag
                  end if
                  i = Key_Ind(scalar(sc) % flux_name, keys, nks)
                  if(i > 0) then
                    scalar(sc) % q(c)=wi*prof(m,i)+(1.-wi)*prof(m+1,i)
                    scalar(sc) % bnd_cond_type(c) = WALLFL
                  end if
                end do

                ! For turbulence models
                if(Turb % model .eq. K_EPS) then

                  i = Key_Ind('KIN',keys,nks)
                  if(i > 0) kin % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)

                  i = Key_Ind('EPS',keys,nks)
                  if(i > 0) eps % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)

                  if(Flow % heat_transfer) then
                    i = Key_Ind('T2',keys,nks)
                    if(i > 0) t2 % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)
                  end if

                end if

                if(Turb % model .eq. K_EPS_ZETA_F .or.  &
                   Turb % model .eq. HYBRID_LES_RANS) then

                  i = Key_Ind('KIN',keys,nks)
                  if(i > 0) kin % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)

                  i = Key_Ind('EPS',keys,nks)
                  if(i > 0) eps % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)

                  i = Key_Ind('ZETA',keys,nks)
                  if(i > 0) zeta % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)

                  i = Key_Ind('F22',keys,nks)
                  if(i > 0) f22 % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)

                  if(Flow % heat_transfer) then
                    i = Key_Ind('T2',keys,nks)
                    if(i > 0) t2 % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)
                  end if

                end if

                if(Turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
                   Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then

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

                  if(Turb % model .eq. RSM_MANCEAU_HANJALIC) then
                    i = Key_Ind('F22', keys, nks)
                    if(i > 0)f22 % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)
                  end if
                end if

                if(Turb % model .eq. SPALART_ALLMARAS .or.  &
                   Turb % model .eq. DES_SPALART) then
                  i = Key_Ind('VIS',keys,nks)
                  if(i > 0) vis % b(c) = wi*prof(m,i) + (1.-wi)*prof(m+1,i)
                end if

              end if  ! (here)
            end do    ! m = 1, n_points-1

          end do      ! Cells_In_Region
        end if        ! plane is defined?
        close(fu)

        !-----------------------------!
        !   Free memory for profile   !
        !-----------------------------!
        deallocate(prof)

      end if  ! boundary defined in a file
    end do

  end do

  !-----------------------------------!
  !                                   !
  !   Read data on synthetic eddies   !
  !                                   !
  !-----------------------------------!
  turb_planes % n_planes = 0
  do bc = Boundary_Regions()  ! imagine there are as many eddies as bcs
    call Control % Position_At_Two_Keys('SYNTHETIC_EDDIES',        &
                                        Grid % region % name(bc),  &
                                        found,                     &
                                        .false.)
    if(found) then
      turb_planes % n_planes = turb_planes % n_planes + 1
      call Control % Read_Int_Item_On ('NUMBER_OF_EDDIES', 24, edd_n, .false.)
      call Control % Read_Real_Item_On('MAX_EDDY_RADIUS',  .2, edd_r, .false.)
      call Control % Read_Real_Item_On('EDDY_INTENSITY',   .1, edd_i, .false.)
      call Eddies_Mod_Allocate(turb_planes % plane(turb_planes % n_planes),  &
                               edd_n,                                        &
                               edd_r,                                        &
                               edd_i,                                        &
                               Flow,                                         &
                               Grid % region % name(bc))
    end if
  end do
  if(turb_planes % n_planes > 0 .and. First_Proc()) then
    print *, '# Found ', turb_planes % n_planes, ' turbulent planes'
  end if

  !---------------------------------------!
  !                                       !
  !                                       !
  !   Copy all "b" values to "n" values   !
  !                                       !
  !                                       !
  !---------------------------------------!
  do bc = Boundary_Regions()
    do c = Cells_In_Region(bc)

      u % n(c) = u % b(c)
      v % n(c) = v % b(c)
      w % n(c) = w % b(c)
      p % n(c) = p % b(c)

      if(Flow % heat_transfer) then
        t % n(c) = t % b(c)
      end if

      if (Flow % with_interface) then
        fun % n(c) = fun % b(c)
      end if

      do sc = 1, Flow % n_scalars
        scalar(sc) % n(c) = scalar(sc) % b(c)
      end do

      if(Turb % model .eq. RSM_MANCEAU_HANJALIC .or.  &
         Turb % model .eq. RSM_HANJALIC_JAKIRLIC) then
        uu  % n(c) = uu  % b(c)
        vv  % n(c) = vv  % b(c)
        ww  % n(c) = ww  % b(c)
        uv  % n(c) = uv  % b(c)
        uw  % n(c) = uw  % b(c)
        vw  % n(c) = vw  % b(c)
        eps % n(c) = eps % b(c)

        if(Turb % model .eq. RSM_MANCEAU_HANJALIC) then
          f22 % n(c) = f22 % b(c)
        end if
      end if

      if(Turb % model .eq. K_EPS) then
        kin % n(c) = kin % b(c)
        eps % n(c) = eps % b(c)
        if(Flow % heat_transfer) then
          t2 % n(c) = t2 % b(c)
        end if
      end if

      if(Turb % model .eq. K_EPS_ZETA_F .or.  &
         Turb % model .eq. HYBRID_LES_RANS) then
        kin  % n(c) = kin  % b(c)
        eps  % n(c) = eps  % b(c)
        zeta % n(c) = zeta % b(c)
        f22  % n(c) = f22  % b(c)
        if(Flow % heat_transfer) then
          t2 % n(c) = t2 % b(c)
        end if
      end if

      if(Turb % model .eq. SPALART_ALLMARAS .or.  &
         Turb % model .eq. DES_SPALART) then
        vis % n(c) = vis % b(c)
      end if

    end do  ! boundary cells
  end do    ! boundary regions

  !------------------------------!
  !   Find the near-wall cells   !
  !------------------------------!
  Grid % cell_near_wall = .false.

  do bc = Boundary_Regions()
    if(Grid % region % type(bc) .eq. WALL   .or.  &
       Grid % region % type(bc) .eq. WALLFL) then
      do s = Faces_In_Region(bc)
        c1 = Grid % faces_c(1,s)
        c2 = Grid % faces_c(2,s)

        Grid % cell_near_wall(c1) = .true.
      end do
    end if

  end do  ! boundary regions

  call Grid % Exchange_Cells_Log(Grid % cell_near_wall)

  end subroutine
