!========================= =====================================================!
  subroutine Row_Sort(Amg, level,  &
                      a, ia, ja,   &
                      iw, ifg,     &  ! created from "kwork"
                      jtr)
!------------------------------------------------------------------------------!
!   Row-sort algorithm for rows of a(level). in detail:
!
!   - Sorts the elements of each row of a(level) such that the strong connections
!     just follow the diagonal element.  One of the strongest is always first
!     (even if partial sorting is performed!).  Non-negative connections are
!     always defined to be weak.  ja(ia(i)) is re-defined to point to the last
!     strong connection of point i (or to ia(i) if there is no strong
!     connection).
!
!   - The strong transpose connections are loaded logically into jtr, i.e.,
!     its pointers are stored in jtr.
!
!   Comments on work space used:
!
!   - ja(ia(i)) ..... i=imin(level),...,imax(level)
!
!   is defined to point to the last strong connection of point i.  Note that
!   the original contents of ja(ia(i)) is overwritten.  (It is put back in
!   subroutine wint.)
!
!   - jtr(j) ........ j=1,iw(imax(level)+iws+1)-1
!   - iw(i)  ........ i=imin(level)+iws,...,imax(level)+iws
!   - ifg(i) ........ i=imin(level),...,imax(level)
!
!   (iws=0, if level=1; iws=imaxw(level-1)+2-imin(level) otherwise)
!
!   jtr is initialized to have same form as ja. jtr(j) contains information
!   on strong transpose connections:  jtr(j) with iw(i+iws)<=j<=iw(i+iws+1)-1
!   points to the strong transpose connections of i.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type) :: Amg
  integer         :: level
  real            :: a(:)
  integer         :: ia(:), ja(:)
  integer         :: iw(:), ifg(:), jtr(:)
!-----------------------------------[locals]-----------------------------------!
  real    :: rs, amn, amx, ast, atmp
  integer :: i, ihi, ii, ilo, imx, itmp, iws, j, jhi, jlo, jmx
  integer :: ndjtr, mdjtr
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  ndjtr = size(jtr, 1)

  !----------------------------------!
  !                                  !
  !   Initialization of work space   !
  !                                  !
  !----------------------------------!
  ilo = Amg % imin(level)
  ihi = Amg % imax(level)
  if(level .ne. 1) then
    iws = Amg % imaxw(level-1) + 2 - ilo
  else
    iws = 0
  end if
  do i = ilo, ihi
    ifg(i) = 0
    ja(ia(i)) = ia(i)
  end do

  !---------------------------!
  !                           !
  !   Partial standard sort   !
  !                           !
  !- - - - - - - - - - - - - -+------------------------------------------!
  !   Note: Non-negative connections are always defined to be weak!      !
  !         The strongest connection is always following the diagonal.   !
  !----------------------------------------------------------------------!
  do i = ilo, ihi
    jlo = ia(i)+1
    jhi = ia(i+1)-1
    if(jhi .lt. jlo) exit

    !------------------------------------------------------------------!
    !   Find strongest connection and sum of off-diagonals| of row i   !
    !------------------------------------------------------------------!
    amx = a(jlo)
    amn = a(jlo)
    jmx = jlo
    if(Amg % ecg1 .ne. 0.0) then
      rs = 0.0
      do j = jlo+1, jhi
        rs = rs+abs(a(j))
        if (a(j).lt.amx) then
          amx = a(j)
          jmx = j
        elseif (a(j).gt.amn) then
          amn = a(j)
        end if
      end do

      !----------------------------------------------------------!
      !   Test for positive off-diagonals / diagonal dominance   !
      !----------------------------------------------------------!
      if(amx .ge. 0.0 .or. rs .le. Amg % ecg1 * a(ia(i))) exit
    else
      do j = jlo+1, jhi
        if (a(j).lt.amx) then
          amx = a(j)
          jmx = j
        elseif (a(j).gt.amn) then
          amn = a(j)
        end if
      end do

      !-------------------------------------!
      !   Test for positive off-diagonals   !
      !-------------------------------------!
      if (amx .ge. 0.0) exit
    end if

    !------------------------------------------------!
    !   Put strongest connection in first position   !
    !------------------------------------------------!
    ast = Amg % ecg2 * amx
    imx = ja(jmx)
    a(jmx)  = a(jlo)
    ja(jmx) = ja(jlo)
    a(jlo)  = amx
    ja(jlo) = imx
    if(amn .gt. ast) then
      jhi = jhi+1

      !------------------------------------------------------!
      !   Decrease jhi until a strong connection is found    !
      !   (if jlo >= jhi stop: all connections are sorted)   !
      !------------------------------------------------------!
      outer: do
        jhi = jhi-1
        if(jlo .ge. jhi) exit outer
        if(a(jhi) .le. ast) then

          !------------------------------------------------------!
          !   Increase jlo until a weak connection is found      !
          !   (if jlo >= jhi stop: all connections are sorted)   !
          !------------------------------------------------------!
          inner: do
            jlo = jlo+1
            if(jlo .ge. jhi) exit outer
            if(a(jlo) .gt. ast) exit inner
          end do inner

          !-----------------------------------!
          !   Interchange a(jhi) and a(jlo)   !
          !-----------------------------------!
          atmp = a(jhi)
          itmp = ja(jhi)
          a(jhi)  = a(jlo)
          ja(jhi) = ja(jlo)
          a(jlo)  = atmp
          ja(jlo) = itmp
        end if
      end do outer
    end if

    !------------------------------------------------------------!
    !   Row sorted --  set ja(ia(i)) to last strong connection   !
    !------------------------------------------------------------!
    ja(ia(i)) = jhi

    !-------------------------------------------------!
    !   Count strong transpose connections in row i   !
    !-------------------------------------------------!
    do j = ia(i)+1, jhi
      ifg(ja(j)) = ifg(ja(j))+1
    end do
  end do

  !-------------------------------------------------------------------!
  !   Initialization of work space for strong transpose connections   !
  !-------------------------------------------------------------------!
  iw(ilo+iws) = 1
  do i = ilo, ihi
    iw(i+iws+1) = iw(i+iws)+ifg(i)
    ifg(i) = iw(i+iws)
  end do
  mdjtr = iw(ihi+iws+1)-1
  if(mdjtr .gt. ndjtr) then
    write(6, '(a)')  &
      ' *** error in Row_Sort: nda too small ***'
    Amg % ierr = AMG_ERR_DIM_A_TOO_SMALL
    return
  end if

  !------------------------------------------------------------!
  !   Load pointers to strong transpose connections into jtr   !
  !------------------------------------------------------------!
  do i = ilo, ihi
    do j = ia(i) + 1, ja(ia(i))
      ii = ja(j)
      jtr(ifg(ii)) = i
      ifg(ii) = ifg(ii)+1
    end do
  end do

  end subroutine
