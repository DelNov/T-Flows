!==============================================================================!
  subroutine Compress_Surf_Vertices(Surf, verbose)
!------------------------------------------------------------------------------!
!   Compresses vertices' list                                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Surf_Type), target :: Surf
  logical                   :: verbose
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: Vert(:)
  type(Elem_Type), pointer :: Elem(:)
  integer,         pointer :: nv, ne
  integer                  :: e, v, n_vert, i_ver, j_ver, nv_tot
  real,    allocatable     :: xv(:), yv(:), zv(:)
  integer, allocatable     :: ni(:), new_n(:)
!==============================================================================!

  ! Take aliases
  nv   => Surf % n_verts
  ne   => Surf % n_elems
  Vert => Surf % Vert
  Elem => Surf % Elem

  ! Check sanity of the elements so far
  do e = 1, ne
    do i_ver = 1, Elem(e) % nv-1
      do j_ver = i_ver+1, Elem(e) % nv
        if(Elem(e) % v(i_ver) .eq. Elem(e) % v(j_ver)) then
          call Message % Error(44,                                    &
                               "Error in the sanity check for  "  //  &
                               "elements. \n  Some element(s)  "  //  &
                               "have duplicate vertices.",            &
                               file=__FILE__, line=__LINE__)
        end if
      end do
    end do
  end do

  if(nv > 0) then
    allocate(xv(nv));     xv    = 0.0
    allocate(yv(nv));     yv    = 0.0
    allocate(zv(nv));     zv    = 0.0
    allocate(ni(nv));     ni    = 0
    allocate(new_n(nv));  new_n = 0

    do v = 1, nv
      xv(v) = Vert(v) % x_n
      yv(v) = Vert(v) % y_n
      zv(v) = Vert(v) % z_n
      ni(v) = v
    end do
    call Sort % Three_Real_Carry_Int(xv, yv, zv, ni)
  end if

  !-----------------------------------------!
  !   Count compressed number of vertices   !
  !-----------------------------------------!
  if(nv > 0) then
    n_vert = 1
    new_n(1) = n_vert
    do v = 2, nv
      if(.not. Math % Approx_Real(xv(v), xv(v-1), NANO)) then
        n_vert = n_vert + 1

      ! xi(v) .eq. xi(v-1)
      else
        if(.not. Math % Approx_Real(yv(v), yv(v-1), NANO)) then
          n_vert = n_vert + 1

        ! xi(v) .eq. xi(v-1) and yi(v) .eq. yi(v-1)
        else
          if(.not. Math % Approx_Real(zv(v), zv(v-1), NANO)) then
            n_vert = n_vert + 1
          end if
        end if
      end if
      new_n(v) = n_vert
    end do
  else
    n_vert = 0
  end if

  !----------------------------------------!
  !   Copy compressed vertex coordinates   !
  !----------------------------------------!
  do v = 1, nv
    Vert(new_n(v)) % x_n = xv(v)
    Vert(new_n(v)) % y_n = yv(v)
    Vert(new_n(v)) % z_n = zv(v)
  end do

  !--------------------------------!
  !   Correct elements' vertices   !
  !--------------------------------!
  if(nv > 0) then
    call Sort % Int_Carry_Int(ni, new_n)

    do e = 1, ne
      do i_ver = 1, Elem(e) % nv
        Elem(e) % v(i_ver) = new_n(Elem(e) % v(i_ver))
      end do
    end do
  end if

  ! Store compressed number of vertices
  nv     = n_vert
  nv_tot = n_vert
  if(Surf % mesh_divided) then
    call Comm_Mod_Global_Sum_Int(nv_tot)
  end if
  if(verbose .and. this_proc < 2) then
    print '(a40,i8)', ' # Compressed number of vertices:       ', nv_tot
  end if

  ! Check sanity of the elements in the end
  do e = 1, ne
    do i_ver = 1, Elem(e) % nv-1
      do j_ver = i_ver+1, Elem(e) % nv
        if(Elem(e) % v(i_ver) .eq. Elem(e) % v(j_ver)) then
          call Message % Error(44,                                          &
                               "Error in the final sanity check for  "  //  &
                               "elements. \n  Some element(s)        "  //  &
                               "have duplicate vertices.",                  &
                               file=__FILE__, line=__LINE__)
        end if
      end do
    end do
  end do

  end subroutine
