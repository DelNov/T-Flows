!==============================================================================!
  subroutine Front_Mod_Compress_Vertices(front, verbose)
!------------------------------------------------------------------------------!
!   Compresses vertices' list                                                  !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Front_Type), target :: front
  logical                  :: verbose
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ne
  integer                  :: e, v, n_vert
  real,    allocatable     :: xv(:), yv(:), zv(:)
  integer, allocatable     :: ni(:), new_n(:)
!==============================================================================!

  ! Take aliases
  nv   => front % n_verts
  ne   => front % n_elems
  vert => front % vert
  elem => front % elem

  ! Check sanity of the elements so far
  do e = 1, ne
    if( (elem(e) % i .eq. elem(e) % j) .or.  &
        (elem(e) % i .eq. elem(e) % k) .or.  &
        (elem(e) % j .eq. elem(e) % k) ) then
      print '(a)',      ' # ERROR in the beginning of Compress_Vertices'
      print '(a,i6,a)', ' # element ', e, 'has same vertices'
    end if
  end do

  allocate(xv(nv));     xv    = 0.0
  allocate(yv(nv));     yv    = 0.0
  allocate(zv(nv));     zv    = 0.0
  allocate(ni(nv));     ni    = 0
  allocate(new_n(nv));  new_n = 0

  do v = 1, nv
    xv(v) = vert(v) % x_n
    yv(v) = vert(v) % y_n
    zv(v) = vert(v) % z_n
    ni(v) = v
  end do
  call Sort_Mod_3_Real_Carry_Int(xv, yv, zv, ni)

  !-----------------------------------------!
  !   Count compressed number of vertices   !
  !-----------------------------------------!
  n_vert = 1
  new_n(1) = n_vert
  do v = 2, nv
    if(.not. Math_Mod_Approx_Real(xv(v), xv(v-1), MICRO)) then
      n_vert = n_vert + 1

    ! xi(v) .eq. xi(v-1)
    else
      if(.not. Math_Mod_Approx_Real(yv(v), yv(v-1), MICRO)) then
        n_vert = n_vert + 1

      ! xi(v) .eq. xi(v-1) and yi(v) .eq. yi(v-1)
      else
        if(.not. Math_Mod_Approx_Real(zv(v), zv(v-1), MICRO)) then
          n_vert = n_vert + 1
        end if
      end if
    end if
    new_n(v) = n_vert
  end do

  !----------------------------------------!
  !   Copy compressed vertex coordinates   !
  !----------------------------------------!
  do v = 1, nv
    vert(new_n(v)) % x_n = xv(v)
    vert(new_n(v)) % y_n = yv(v)
    vert(new_n(v)) % z_n = zv(v)
  end do

  !--------------------------------!
  !   Correct elements' vertices   !
  !--------------------------------!
  call Sort_Mod_Int_Carry_Int(ni, new_n)

  do e = 1, ne
    elem(e) % i = new_n(elem(e) % i)
    elem(e) % j = new_n(elem(e) % j)
    elem(e) % k = new_n(elem(e) % k)
  end do

  ! Store compressed number of vertices
  nv = n_vert
  if(verbose) print *, '# Compressed number of vertices: ', nv

  ! Check sanity of the elements in the end
  do e = 1, ne
    if( (elem(e) % i .eq. elem(e) % j) .or.  &
        (elem(e) % i .eq. elem(e) % k) .or.  &
        (elem(e) % j .eq. elem(e) % k) ) then
      print '(a)',      ' # ERROR in the end of Compress_Vertices'
      print '(a,i6,a)', ' # element', e, ' has same vertices'
      stop
    end if
  end do

  end subroutine
