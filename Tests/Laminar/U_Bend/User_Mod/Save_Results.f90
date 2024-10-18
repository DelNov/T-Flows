!==============================================================================!
  subroutine User_Mod_Save_Results(Flow, Turb, Vof, Swarm, domain)
!------------------------------------------------------------------------------!
!   This subroutine is called each RESULTS_SAVE_INTERVAL (set in control       !
!   file), at the end of a simulation and after 'save_now' command.            !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target   :: Flow
  type(Turb_Type),  target   :: Turb
  type(Vof_Type),   target   :: Vof
  type(Swarm_Type), target   :: Swarm
  integer,          optional :: domain
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  type(Bulk_Type), pointer :: bulk
  type(Var_Type),  pointer :: u, v, w, t
  integer                  :: c, c1, c2, i_nod, n, n_inner, n_outer, reg, s
  real                     :: x_upp, x_low, r_outer, r_inner, s_coord
  real                     :: tan_alpha, alpha_deg
  real,    allocatable     :: s_inner(:), s_outer(:), buffer(:)
  integer, allocatable     :: s_inner_c1(:), s_outer_c1(:)
  integer, allocatable     :: s_inner_c2(:), s_outer_c2(:)
  integer                  :: fu
  character(SL)            :: res_name
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid
  bulk => Flow % bulk
  u    => Flow % u
  v    => Flow % v
  w    => Flow % w
  t    => Flow % t

  !-----------------------------!
  !                             !
  !                             !
  !   We are in U_BEND domain   !
  !                             !
  !                             !
  !-----------------------------!
  if(Grid % name .eq. 'U_BEND') then

    !-------------------------------------------------------------!
    !                                                             !
    !   Find the smallest x, it will be used to compute s_coord   !
    !                                                             !
    !-------------------------------------------------------------!
    x_upp = HUGE
    x_low = HUGE

    do c = Cells_In_Domain()
      do i_nod = 1, Grid % cells_n_nodes(c)  ! through cell's nodes
        n = Grid % cells_n(i_nod, c)

        if(Grid % yn(n) > 0.0) x_upp = min(Grid % xn(n), x_upp)
        if(Grid % yn(n) < 0.0) x_low = min(Grid % xn(n), x_low)
      end do
    end do

    print *, "x upp = ", x_upp
    print *, "x low = ", x_low

    !--------------------------------------------------------!
    !                                                        !
    !   Count the number of cells in outer and inner walls   !
    !      and work out the outer and inner radius too.      !
    !                                                        !
    !--------------------------------------------------------!
    n_outer = 0;
    r_outer = 0;
    do reg = Boundary_Regions()
      if(Grid % region % name(reg) .eq. 'WALLS_OUTER') then
        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1, s)
          c2 = Grid % faces_c(2, s)
          Assert(c2 .lt. 0)
          if(Math % Approx_Real(0.0, Grid % zc(c2))) n_outer = n_outer + 1
          r_outer = max(r_outer, Grid % yc(c2))
        end do  ! faces in region
      end if    ! region is outer wall
    end do      ! through regions

    n_inner = 0;
    r_inner = 0;
    do reg = Boundary_Regions()
      if(Grid % region % name(reg) .eq. 'WALLS_INNER') then
        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1, s)
          c2 = Grid % faces_c(2, s)
          Assert(c2 .lt. 0)
          if(Math % Approx_Real(0.0, Grid % zc(c2))) n_inner = n_inner + 1
          r_inner = max(r_inner, Grid % yc(c2))
        end do  ! faces in region
      end if    ! region is outer wall
    end do      ! through regions


    print *, "r_outer = ", r_outer
    print *, "r_inner = ", r_inner
    Assert(n_outer .eq. n_inner)

    !--------------------------------------!
    !                                      !
    !   Allocate memory for local arrays   !
    !                                      !
    !--------------------------------------!
    allocate(s_outer   (n_outer));  s_outer   (:) = 0.0
    allocate(s_inner   (n_inner));  s_inner   (:) = 0.0
    allocate(s_outer_c1(n_outer));  s_outer_c1(:) = 0
    allocate(s_inner_c1(n_inner));  s_inner_c1(:) = 0
    allocate(s_outer_c2(n_outer));  s_outer_c2(:) = 0
    allocate(s_inner_c2(n_inner));  s_inner_c2(:) = 0

    !-------------------------------------------------!
    !                                                 !
    !   Browse through outer and inner wall regions   !
    !                                                 !
    !-------------------------------------------------!
    n_outer = 0
    n_inner = 0
    do reg = Boundary_Regions()

      !-----------------------------!
      !   Outer wall region found   !
      !-----------------------------!
      if(Grid % region % name(reg) .eq. 'WALLS_OUTER') then

        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1, s)
          c2 = Grid % faces_c(2, s)
          Assert(c2 .lt. 0)

          ! Work out s and n_inner if you are at the center line
          if(Math % Approx_Real(0.0, Grid % zc(c2))) then

            ! If you are in one of the legs
            if(Grid % xc(c2) < 0) then
              if(Grid % yc(c2) > 0) then  ! you are in the upper leg
                s_coord = abs(x_upp) + Grid % xc(c2)
              end if
              if(Grid % yc(c2) < 0) then  ! you are in the lower leg
                s_coord = abs(x_upp) + r_outer * PI + abs(Grid % xc(c2))
              end if

            ! You are in the turn
            else
              tan_alpha = Grid % xc(c2) / Grid % yc(c2)
              alpha_deg = atan(tan_alpha) * 180 / PI
              if(alpha_deg < 0)  alpha_deg = 180.0 + alpha_deg
              s_coord = abs(x_upp) + alpha_deg / 180 * r_outer * PI
            end if

            ! Store the next entry along the "s" coordinate
            n_outer = n_outer + 1
            s_outer   (n_outer) = s_coord
            s_outer_c1(n_outer) = c1
            s_outer_c2(n_outer) = c2

          end if  ! in the centerline

        end do    ! faces in the region

      !-----------------------------!
      !   Inner wall region found   !
      !-----------------------------!
      else if(Grid % region % name(reg) .eq. 'WALLS_INNER') then

        do s = Faces_In_Region(reg)
          c1 = Grid % faces_c(1, s)
          c2 = Grid % faces_c(2, s)
          Assert(c2 .lt. 0)

          ! If you are at the center line
          if(Math % Approx_Real(0.0, Grid % zc(c2))) then

            ! If you are in one of the legs
            if(Grid % xc(c2) < 0) then
              if(Grid % yc(c2) > 0) then  ! you are in the upper leg
                s_coord = abs(x_upp) + Grid % xc(c2)
              end if
              if(Grid % yc(c2) < 0) then  ! you are in the lower leg
                s_coord = abs(x_upp) + r_inner * PI + abs(Grid % xc(c2))
              end if

            ! You are in the turn
            else
              tan_alpha = Grid % xc(c2) / Grid % yc(c2)
              alpha_deg = atan(tan_alpha) * 180 / PI
              if(alpha_deg < 0)  alpha_deg = 180 + alpha_deg
              s_coord = abs(x_upp) + alpha_deg / 180 * r_inner * PI
            end if

            ! Store the next entry along the "s" coordinate
            n_inner = n_inner + 1
            s_inner   (n_inner) = s_coord
            s_inner_c1(n_inner) = c1
            s_inner_c2(n_inner) = c2

          end if  ! in the centerline

        end do    ! faces in the region

        print *, "n_outer = ", n_outer
        print *, "n_inner = ", n_inner

      end if  ! regions is outer or inner wall
    end do    ! through regions

    !-------------------------------------------------!
    !                                                 !
    !   Sort cell indices along the "s" coordinates   !
    !                                                 !
    !-------------------------------------------------!
    call Sort % Real_Carry_Two_Int(s_outer, s_outer_c1, s_outer_c2)
    call Sort % Real_Carry_Two_Int(s_inner, s_inner_c1, s_inner_c2)

    !------------------------------------------------!
    !                                                !
    !   Now save some values along the centerlines   !
    !                                                !
    !------------------------------------------------!
    allocate(buffer(n_outer));  buffer(:) = 0.0

    !-------------------------------------------------!
    !   Save x velocity component on the outer path   !
    !-------------------------------------------------!
    do n = 1, n_outer
      c1 = s_outer_c1(n)  ! if you need value from inside
      c2 = s_outer_c2(n)  ! if you need value from the wall

      buffer(n) = Flow % u % n(c1)  ! could have used alias u as well
    end do

    call File % Set_Name(res_name,                        &
                         domain    = domain,              &
                         time_step = Time % Curr_Dt(),    &
                         appendix  = '-u-along-outer-s',  &
                         extension = '.dat')
    call File % Open_For_Writing_Ascii(res_name, fu)

    do n = 1, n_outer
      write(fu, '(i6, es15.5, es15.5)')  n, s_outer(n), buffer(n)
    end do

    close(fu)

    !-------------------------------------------------!
    !   Save x velocity component on the inner path   !
    !-------------------------------------------------!
    do n = 1, n_inner
      c1 = s_inner_c1(n)  ! if you need value from inside
      c2 = s_inner_c2(n)  ! if you need value from the wall

      buffer(n) = Flow % u % n(c1)  ! could have used alias u as well
    end do

    call File % Set_Name(res_name,                        &
                         domain    = domain,              &
                         time_step = Time % Curr_Dt(),    &
                         appendix  = '-u-along-inner-s',  &
                         extension = '.dat')
    call File % Open_For_Writing_Ascii(res_name, fu)

    do n = 1, n_inner
      write(fu, '(i6, es15.5, es15.5)')  n, s_inner(n), buffer(n)
    end do

    close(fu)

    !-----------------------------------------------!
    !   Save pressure component on the outer path   !
    !-----------------------------------------------!
    do n = 1, n_outer
      c1 = s_outer_c1(n)  ! if you need value from inside
      c2 = s_outer_c2(n)  ! if you need value from the wall

      buffer(n) = Flow % p % n(c1)  ! could have used alias u as well
    end do

    call File % Set_Name(res_name,                        &
                         domain    = domain,              &
                         time_step = Time % Curr_Dt(),    &
                         appendix  = '-p-along-outer-s',  &
                         extension = '.dat')
    call File % Open_For_Writing_Ascii(res_name, fu)

    do n = 1, n_outer
      write(fu, '(i6, es15.5, es15.5)')  n, s_outer(n), buffer(n)
    end do

    close(fu)

    !-----------------------------------------------!
    !   Save pressure component on the inner path   !
    !-----------------------------------------------!
    do n = 1, n_inner
      c1 = s_inner_c1(n)  ! if you need value from inside
      c2 = s_inner_c2(n)  ! if you need value from the wall

      buffer(n) = Flow % p % n(c1)  ! could have used alias u as well
    end do

    call File % Set_Name(res_name,                        &
                         domain    = domain,              &
                         time_step = Time % Curr_Dt(),    &
                         appendix  = '-p-along-inner-s',  &
                         extension = '.dat')
    call File % Open_For_Writing_Ascii(res_name, fu)

    do n = 1, n_inner
      write(fu, '(i6, es15.5, es15.5)')  n, s_inner(n), buffer(n)
    end do

    close(fu)

  end if

  end subroutine
