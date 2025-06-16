!==============================================================================!
  subroutine Store_Level(Amg, level, max_lev,  &
                         a, u, f, ia, ja,   &  ! defining system
                         iw, icg, ifg)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type) :: Amg
  integer         :: level, max_lev
  real            :: a(:), u(:), f(:)
  integer         :: ia(:), ja(:)
  integer         :: iw(:), icg(:), ifg(:)
!-----------------------------------[locals]-----------------------------------!
  integer :: i, iaux, j
  integer :: i_loc, j_loc, k_loc, k, n, nw, nnz, nnw
  integer :: if, test_level, nc
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  !----------------------------------------------!
  !                                              !
  !   Part I: ia(Amg % imax(level)+1) is fixed   !
  !                                              !
  !----------------------------------------------!

  ! See comment in source "Coarsening.f90" at line 180
  iaux = ia(Amg % imax(level)+1)
  ia(Amg % imax(level)+1) = iw(Amg % iminw(level))

  ! Calculate level dimensions
  n   = Amg % imax (level) - Amg % imin (level) + 1
  nw  = Amg % imaxw(level) - Amg % iminw(level) + 1
  nnz = ia(Amg % imax (level)+1) - ia(Amg % imin (level))
  nnw = iw(Amg % imaxw(level)+1) - iw(Amg % iminw(level))

  ! Size of the coarser level (for checking)
  nc = Amg % imax(level+1) - Amg % imin(level+1) + 1

# define TO_MAKE_SURE_I_GET_IT

# ifdef  TO_MAKE_SURE_I_GET_IT
  print *, __FILE__, __LINE__
  print *, ' level = ', level
  print *, ' unknowns range:         ', Amg % imin (level),  &
                        ' ---> ',       Amg % imax (level),  &
                        ' tot: ',       Amg % imax (level) - Amg % imin (level) + 1
  print *, ' matrix storage range:   ', ia(Amg % imin(level)),    &
                        ' ---> ',       ia(Amg % imax(level)+1),  &  ! +1, stores the end
                        ' tot: ',       ia(Amg % imax(level)+1) - ia(Amg % imin(level))
  print *, ' operator storage range: ', iw(Amg % iminw(level)),  &
                        ' ---> ',       iw(Amg % imaxw(level)),  &
                        ' tot: ',       iw(Amg % imaxw(level)) - iw(Amg % iminw(level)) + 1
# endif

  Amg % lev(level) % n   = n
  Amg % lev(level) % nw  = nw
  Amg % lev(level) % nnz = nnz  ! needed? yes!

  ! Allocate memory for the level
  allocate(Amg % lev(level) % a(nnz),     &
           Amg % lev(level) % u(n),       &
           Amg % lev(level) % f(n),       &
           Amg % lev(level) % ia(n+1),    &
           Amg % lev(level) % ja(nnz),    &
           Amg % lev(level) % icg(n),     &
           Amg % lev(level) % ifg_1(n),   &
           Amg % lev(level) % ifg_2(nw),  &
           Amg % lev(level) % w(nnw),     &
           Amg % lev(level) % jw(nnw),    &
           Amg % lev(level) % iw(nw+1))  ! to store the end

  ! Initialize values to zero
  Amg % lev(level) % a(:)     = 0.0
  Amg % lev(level) % u(:)     = 0.0
  Amg % lev(level) % f(:)     = 0.0
  Amg % lev(level) % ia(:)    = 0
  Amg % lev(level) % ja(:)    = 0
  Amg % lev(level) % icg(:)   = 0
  Amg % lev(level) % ifg_1(:) = 0
  Amg % lev(level) % ifg_2(:) = 0
  Amg % lev(level) % iw(:)    = 0

  !---------------------------------------------!
  !   Copy vectors u and f to level's storage   !
  !---------------------------------------------!
  do i = Amg % imin(level), Amg % imax(level)
    i_loc = i - Amg % imin(level) + 1
    Amg % lev(level) % u(i_loc) = u(i)
    Amg % lev(level) % f(i_loc) = f(i)
  end do

  !-----------------------!
  !   Copy and shift ia   !
  !-----------------------!
  do i = Amg % imin(level), Amg % imax(level) + 1
    i_loc = i - Amg % imin(level) + 1
    Amg % lev(level) % ia(i_loc) = ia(i) - ia(Amg % imin(level)) + 1
  end do

  if(Amg % lev(level) % ia(n+1) - 1 /= nnz) then
    print *, "IA Copy Error: Expected nnz=", nnz, " Got=", &
             Amg % lev(level) % ia(n+1)-1
    stop
  end if

  !--------------------------------------------------------------!
  !   Copy a and ja (only the parts which hold linear systems)   !
  !--------------------------------------------------------------!
  j_loc = 0
  do i = Amg % imin(level), Amg % imax(level)
    do j = ia(i), ia(i+1) - 1
      Assert(Amg % Level_Of_Cell(ja(j)) == level)  ! stays at its own level
      j_loc = j_loc + 1
      Amg % lev(level) %  a(j_loc) =  a(j)
      Amg % lev(level) % ja(j_loc) = ja(j) - Amg % imin(level) + 1
    end do
  end do

  !-------------------------------------------------------!
  !   Copy icg ... now this one is a bit more elaborate   !
  !-------------------------------------------------------!
  do i = Amg % imin(level), Amg % imax(level)

    ! Retreive local index, at this level
    i_loc = i - Amg % imin(level) + 1

    if(icg(i) .ne. -AMG_BIG_INTEGER) then

      ! Negative icg's point to current level
      if(icg(i) .lt. 0) then
        if(level.lt.max_lev)  Assert(Amg % Level_Of_Cell(-icg(i)).eq.level)
        j_loc = -icg(i) - Amg % imin(level) + 1
        if(level.lt.max_lev)  Assert(j_loc .gt. 0)  ! must be > 0
        if(level.lt.max_lev)  Assert(j_loc .le. n)  ! must fall in this level
        Amg % lev(level) % icg(i_loc) = -j_loc
      end if

      ! Positive icg's point to coarser level
      if(icg(i) .gt. 0) then
        if(level.lt.max_lev)  Assert(Amg % Level_Of_Cell(icg(i)).eq.level+1)
        j_loc = icg(i) - Amg % imin(level + 1) + 1
        if(level.lt.max_lev)  Assert(j_loc .gt.  0)  ! must be > 0
        if(level.lt.max_lev)  Assert(j_loc .le. nc)  ! must fall in the coarser
        Amg % lev(level) % icg(i_loc) = +j_loc
      end if

    ! Case when icg(i) is equal to -AMG_BIG_INTEGER, just copy
    else
      Amg % lev(level) % icg(i_loc) = -AMG_BIG_INTEGER
    end if
  end do

   !------------------------------------------------!
   !   Copy vectors ifg, fine grid level pointers   !
   !   (I strongly believe these are the pointers   !
   !   to the finer grid and Assert confirms it.)   !
   !------------------------------------------------!
   if(level .gt. 1) then
     do i = Amg % imin(level), Amg % imax(level)
       Assert(Amg % Level_Of_Cell(ifg(i)) == level-1)  ! it points to finer grid
       i_loc = i - Amg % imin(level) + 1
       Amg % lev(level) % ifg_1(i_loc) = ifg(i) - Amg % imin(level-1) + 1
                                                      ! index from finer grid
     end do
   end if

   ! This is very dubious, I am not sure if it makes any sense
   do k = Amg % iminw(level), Amg % imaxw(level)
     k_loc = k - Amg % iminw(level) + 1
     Amg % lev(level) % ifg_2(k_loc) = ifg(k) - Amg % imin(level) + 1
   end do                                        ! not sure about this :-( !

  !----------------------------------------------------!
  !   End of Part I: restore ia(Amg % imax(level)+1)   !
  !----------------------------------------------------!
  ia(Amg % imax(level)+1) = iaux

  !-----------------------------------------------!
  !                                               !
  !   Part II: iw(Amg % imax(level)+1) is fixed   !
  !                                               !
  !-----------------------------------------------!

  IF(LEVEL .eq. 7) PRINT *, __FILE__, 'iw(Amg % imaxw(level)+1) = ', iw(Amg % imaxw(level)+1)
  IF(LEVEL .eq. 7) PRINT *, __FILE__, 'ia(Amg % imin(level+1))  = ', ia(Amg % imin(level+1))
  if(level .lt. max_lev) then
    iaux = iw(Amg % imaxw(level)+1)
    iw(Amg % imaxw(level)+1) = ia(Amg % imin(level+1))
  end if

  !------------------------------------------!
  !   Test run from Interpolate_Correction   !
  !------------------------------------------!
  if(level .lt. max_lev) then
    do i = Amg % iminw(level), Amg % imaxw(level)
      if = ifg(i)              ! at max_lev, this holds some garbage
      Assert(Amg % Level_Of_Cell(if) == level)
      do j = iw(i), iw(i+1)-1  ! j browses through work space, THIS workspace
        Assert(Amg % Level_Of_Cell(ja(j)) == level+1)  ! coarser cell
      end do
    end do
  end if

  if(level .lt. max_lev) then
    j_loc = 0
    do i = Amg % iminw(level), Amg % imaxw(level)
      if = ifg(i)              ! at max_lev, this holds some garbage
      Assert(Amg % Level_Of_Cell(if) == level)
      k_loc = i - Amg % iminw(level) + 1
      Amg % lev(level) % iw(k_loc) = iw(i) - iw(Amg % iminw(level)) + 1
      do j = iw(i), iw(i+1)-1
        j_loc = j_loc + 1
        Assert(Amg % Level_Of_Cell(ja(j)) == level+1)
        Amg % lev(level) %  w(j_loc) =  a(j)  ! copy factor to "w"
        Amg % lev(level) % jw(j_loc) = ja(j) - Amg % imin(level+1) + 1  ! position at the coarser level
      end do
    end do
    Amg % lev(level) % iw(nw+1) = iw(Amg % imaxw(level)+1) - iw(Amg % iminw(level)) + 1

  end if

  !------------------------------------------------------!
  !   End of Part II: restore iw(Amg % imaxw(level)+1)   !
  !------------------------------------------------------!
  iw(Amg % imaxw(level)+1) = iaux

  !----------------------------------!
  !                                  !
  !   Start of each color handling   !
  !                                  !
  !----------------------------------!
  i = Amg % start_of_color(level)
  if(level .lt. max_lev)  Assert(Amg % Level_Of_Cell(i) == level)

  Amg % lev(level) % start_of_color = i - Amg % imin(level) + 1
  j = Amg % lev(level) % start_of_color
  if(level .lt. max_lev) then
    Assert(j .gt. 0)
    Assert(j .le. n)
  end if

  end subroutine
