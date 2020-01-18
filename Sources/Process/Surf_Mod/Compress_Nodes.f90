!==============================================================================!
  subroutine Surf_Mod_Compress_Nodes(surf, verbose)
!------------------------------------------------------------------------------!
!   Compresses nodes' list                                                     !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Surf_Type), target :: surf
  logical                 :: verbose
!-----------------------------------[Locals]-----------------------------------!
  type(Vert_Type), pointer :: vert(:)
  type(Elem_Type), pointer :: elem(:)
  integer,         pointer :: nv, ne
  integer                  :: e, v, n_vert
  integer, allocatable     :: xi(:), yi(:), zi(:)
  integer, allocatable     :: ni(:), new_n(:)
  integer, parameter       :: INT_SCL = 100000000  ! increased this for VOF
  integer, parameter       :: INT_TOL =       100
!==============================================================================!

  ! Take aliases
  nv   => surf % n_verts
  ne   => surf % n_elems
  vert => surf % vert
  elem => surf % elem

  ! Check sanity of the elements so far
  do e = 1, ne
    if( (elem(e) % i .eq. elem(e) % j) .or.  &
        (elem(e) % i .eq. elem(e) % k) .or.  &
        (elem(e) % j .eq. elem(e) % k) ) then
      print '(a)',      ' # ERROR in the beginning of Compress_Nodes'
      print '(a,i6,a)', ' # element ', e, 'has same nodes'
    end if
  end do

  allocate(xi(nv));     xi    = 0
  allocate(yi(nv));     yi    = 0
  allocate(zi(nv));     zi    = 0
  allocate(ni(nv));     ni    = 0
  allocate(new_n(nv));  new_n = 0

  do v = 1, nv
    xi(v) = nint(vert(v) % x_n * INT_SCL)
    yi(v) = nint(vert(v) % y_n * INT_SCL)
    zi(v) = nint(vert(v) % z_n * INT_SCL)
    ni(v) = v
  end do
  call Sort_Mod_3_Int_Carry_Int(xi, yi, zi, ni)
! do v = 1, nv
!   WRITE(100, '(3i6,3i16)') v, ni(v), new_n(v), xi(v), yi(v), zi(v)
! end do

  !-----------------------------------------!
  !   Count compressed number of vertices   !
  !-----------------------------------------!
  n_vert = 1
  new_n(1) = n_vert
  do v = 2, nv
    if(abs(xi(v) - xi(v-1)) > INT_TOL) then
      n_vert = n_vert + 1

    ! xi(v) .eq. xi(v-1)
    else
      if(abs(yi(v) - yi(v-1)) > INT_TOL) then
        n_vert = n_vert + 1

      ! xi(v) .eq. xi(v-1) and yi(v) .eq. yi(v-1)
      else
        if(abs(zi(v) - zi(v-1)) > INT_TOL) then
          n_vert = n_vert + 1
        end if
      end if
    end if
    new_n(v) = n_vert
  end do

! do v = 1, nv
!   WRITE(150, '(3i6,3f12.6)') v, ni(v), new_n(v),  &
!              real(xi(v))/INT_SCL, real(yi(v))/INT_SCL, real(zi(v))/INT_SCL
! end do

  !----------------------------------------!
  !   Copy compressed vertex coordinates   !
  !----------------------------------------!
  do v = 1, nv
    vert(new_n(v)) % x_n = real(xi(v))/INT_SCL
    vert(new_n(v)) % y_n = real(yi(v))/INT_SCL
    vert(new_n(v)) % z_n = real(zi(v))/INT_SCL
  end do

  !-----------------------------!
  !   Correct elements' nodes   !
  !-----------------------------!
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
      print '(a)',      ' # ERROR in the end of Compress_Nodes'
      print '(a,i6,a)', ' # element', e, ' has same nodes'
      stop
    end if
  end do

  end subroutine
