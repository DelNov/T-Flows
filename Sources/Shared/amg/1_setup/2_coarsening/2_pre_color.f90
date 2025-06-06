!==============================================================================!
  subroutine pre_color(amg, level,    &
                       ia, ja,        &  ! defining system
                       iw, icg, ifg,  &  ! created from "kwork"
                       ir,            &  ! created from "iwork"
                       jtr,           &
                       iias)
!------------------------------------------------------------------------------!
!     Pre-coloring algorithm for grid k. this is the version as
!     described in ruge/stueben (bristol). the goal is to obtain quickly
!     a tentative set of c-points with the following properties:
!
!       - the c-points are only weakly connected among each other;
!       - each f-point has (at least) one strong connection to a c-point
!         (except for some exceptional f-points, e.g., those which do
!         not have a strong connection at all: the "forced f-points").
!
!     on exit, icg(i) (imin(k)<=i<=imax(k)) contains all the information
!     on the coloring. in detail:
!
!        icg(i) = 1: i is c-point;
!        icg(i) = 0: i is forced f-point, i.e., it is an f-point which
!                    has no strong connection at all;
!        icg(i) =-1: i is "regular" f-point, i.e., i has at least one
!                    strong connection to a c-point;
!        icg(i) =-2: i is "exceptional" f-point, i.e., i has no strong
!                    connection to a c-point. furthermore, no regular
!                    f-point is strongly connected to i and all points
!                    with icg=-2 have no strong connection among each
!                    other. (note: these points don't contribute from
!                    any c-point. also, it does not make sense to make
!                    these points c-points as no f-point contributes
!                    from them.)
!
!     Comments on input
!
!     - ja(ia(i)) ----- i=imin(k),...,imax(k)
!
!     used as defined in row_sort.  not changed.
!
!     - jtr(j)--------- j=1,iw(imax(k)+iws+1)-1
!     - iw(i) --------- i=imin(k)+iws,...,imax(k)+iws
!
!     (iws=0, if k=1; iws=imaxw(k-1)+2-imin(k) otherwise)
!
!     used as defined in row_sort. not changed.
!
!     Comments on work space used
!
!     - icg(i) -------- i=imin(k),...,imax(k)
!
!     used to distinguish f-, c- and u- (undecided) points.
!     on entry, icg is defined to be 0 for forced f-points and -2
!     for u-pnts. on exit: see above.
!
!     - ifg(i) -------- i=imin(k),...,imax(k)
!
!     defined to be a measure for importance of making the u-pnt i
!     a c-point. this measure is a value between jval0 and jvalx with
!
!       jval0=imax(k)+npts+2,   npts=imax(k)-imin(k)+1;
!       jvalx=jval0+4*ntrav+1,  ntrav=avverage number of strong trans-
!                                     pose connections
!
!     the higher this value, the more important is it to make u-point i
!     a c-point. also, the value jv=ifg(i) is just the "origin" of
!     a list which connects all points having the same measure jv of
!     importance.
!
!     - icg(ii) -------- ii=imax(k)+2,...,jval0-1
!     - ifg(ii) -------- ii=imax(k)+2,...,jval0-1
!
!     left and right stack pointers in a doubly linked list. there
!     are several such lists, each of them linking u-points which
!     have the same measure of importance jv to be made c-points.
!     each of the values jv may be regarded as the "origin" of such
!     a list described in the following.
!
!     - icg(jv) -------- jv=jval0,...,jvalx
!     - ifg(jv) -------- jv=jval0,...,jvalx
!
!     used as pointers into the lists to indicate its beginning and its
!     end: ifg(jv) points to the first point in the list and icg(jv)
!     points to the last one. the rough picture of the list which
!     corresponds to a value jval0 <= jv <= jvalx looks as  follows:
!
!
!                                  ifg
!       ----------------------------------------------------------
!       |                                                        |
!       |            ifg           ifg           ifg             |
!       ---> ------ ----> ------- ----> ------- ----> ------- ----
!            |    |       |     |       |     |       |     |
!            | jv |       | ii1 |       | ii2 |       | ii3 |
!            |    |  icg  |     |  icg  |     |  icg  |     |
!       ---- ------ <---- ------- <---- ------- <---- ------- <---
!       |                                                        |
!       |                          icg                           |
!       ----------------------------------------------------------
!
!
!     here, ii1, ii2,... logically represent the "physical" points
!     i1, i2,... (which are values between imin(k) and imax(k)).
!     the interconnection between these two representations of the
!     same points is given by the relation
!
!                          ii := i+npts+1.
!
!     an empty list is characterized by icg(jv)=jv, ifg(jv)=jv.
!     obviously, it is quite easy to remove or add points to the list.
!     such re-arrangements are necessary as (in the coloring part of
!     the coarsening algorithm) undecided points become c- or f-points
!     (they have to be removed from their list) and the measure values
!     jv of some points change during the algorithm. such points have
!     to be moved from one list to another one. in order to keep track
!     of those points which have to be re-arranged, a reset-stack is
!     used which contains all these points i (i.e. their "logical"
!     numbers ii). the pointer ir is used to point from one point in
!     the stack to the previous one.
!
!     - ir(j) --------- j=1,...,npts
!
!     used below as stack-pointer for points to be reset in performing
!     the coloring algorithm. itop is the top-of-stack pointer.
!     the stack is empty if itop=-1.
!
!     the global picture of the pointers is sketched in the following
!
!
!                                          list end
!                                       <-------------|
!                                         list start  |
!                                       <------------||
!                                                    ||
!                  |-- ii=i+npts+1 ->|            ifg||icg
!                  |                 |               ||
!                  |                 |               ||
!       imin(k)          imax(k)           jval0             jvalx
!          |       i        |        ii      |       jv        |
!          |=======*========|========*=======|========*========|
!          |                |                |                 |
!               grid k          logical num      list origins
!                                              (measure values)
!                 ||                                  |
!                 ||         ifg                      |
!                 ||--------------------------------->|
!                 |
!                 |     |---> 1 (c)
!                 |     |
!                 | icg |---> 0 (ff)
!                 |-----|
!                       |--->-1 (f)
!                       |
!                       |--->-2 (u)
!
!
!          1                npts
!          |        j        |
!          |========*========|
!          |                 |
!           shifted grid k
!                   |
!                   | ir
!                   |-----> reset stack (itop = top-of-stack)
!
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(amg_type) :: amg
  integer         :: level
  integer         :: ia(:), ja(:)
  integer         :: iw(:), ifg(:), icg(:), ir(:), jtr(:)
  integer         :: iias
!-----------------------------------[locals]-----------------------------------!
  integer :: i,ii,iii,iiii, ilo,ilo1, ihi, ic, iip, in, itop, iws, iic
  integer :: j, jcnbhi, jj, jv, jval0, jvalx
  integer :: npts, npts1, nscn, nscnmx, ntrlim, ndicg, ndir
  logical :: picknext
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  ndicg = size(icg, 1)
  ndir  = size(ir,  1)

  !-----------------!
  !   Preparation   !
  !-----------------!
  call amg % timer_start()
  ilo = amg % imin(level)
  ihi = amg % imax(level)
  if (level.ne.1) then
    iws = amg % imaxw(level-1)+2-ilo
  else
    iws = 0
  endif
  ilo1 = ilo-1
  npts = amg % imax(level) - amg % imin(level)+1
  npts1 = npts+1
  ntrlim = 2*(iw(ihi+iws+1)-iw(ilo+iws))/npts
  nscnmx = 0
  jval0 = ihi+npts1+1
  jvalx = jval0+2*ntrlim+1
  if(jvalx .gt. ndicg) then
    write(6, '(a)')  &
      ' *** error in pre_color: ndw too small ***'
    amg % ierr = AMG_ERR_DIM_ICG_TOO_SMALL
    return
  end if
  if(npts .gt. ndir) then
    write(6, '(a)')  &
      ' *** error in pre_color: ndu too small ***'
    amg % ierr = AMG_ERR_DIM_U_TOO_SMALL
    return
  end if

  !---------------------------------------------------------!
  !   Put initial "measure" for each point i into ifg(i).   !
  !---------------------------------------------------------!
  do i = ilo, ihi
    nscn = iw(i+iws+1)-iw(i+iws)
    if (nscn.le.ntrlim) then
      ifg(i) = jval0+nscn
    else
      ifg(i) = jvalx
    endif
  end do

  !----------------------------------------------------------!
  !   Set circularly linked lists and reset-stack to empty   !
  !----------------------------------------------------------!
  itop = -1
  do j = jval0, jvalx
    icg(j) = j
    ifg(j) = j
  end do

  !-----------------------------------------------------------------------!
  !   Put all u-points of grid level into lists (i.e. no forced f-points) !
  !   add points always to the end of their corresponding list.           !
  !   in the following, jcnbhi denotes the actual highest measure value   !
  !   (among u-points).                                                   !
  !-----------------------------------------------------------------------!
  jcnbhi = 0
  do i = ilo, ihi
    if(ja(ia(i)).gt.ia(i)) then
      icg(i) = -2
      ii = i+npts1
      ir(i-ilo1) = 0
      jv = ifg(i)
      if (jv.gt.jcnbhi) jcnbhi = jv
      icg(ii) = icg(jv)
      ifg(ii) = jv
      icg(jv) = ii
      ifg(icg(ii)) = ii
    else
      icg(i) = 0
    end if
  end do

  !----------------------------------------------------------------------!
  !   Pre-coloring                                                       !
  !                                                                      !
  !   Pick a u-point ic with maximal measure as given by jcnbhi.         !
  !   (take the first from corresponding list: first-in/first-out)       !
  !   make that point a c-point and remove it from lists. then make      !
  !   all strong transpose u-connections f-points and remove them from   !
  !   their lists. update the measure of importance for u-points to      !
  !   become c-points and add these points to the reset-stack. finally,  !
  !   re-arrange the lists by sweeping through the reset stack. then     !
  !   pick another u-point ic.                                           !
  !                                                                      !
  !   if list corresponding to jcnbhi is empty, go to next lower value.  !
  !   pre-colouring is finished if jcnbhi<=jval0. all u-points left at   !
  !   that time, will be regarded as f-points later.                     !
  !----------------------------------------------------------------------!
  do
    picknext = .false.
    if(jcnbhi .le. jval0) then
      call amg % timer_stop(2)
      return
    end if
    iic = ifg(jcnbhi)
    if(iic .eq. jcnbhi) then
      jcnbhi = jcnbhi-1
      ! Setting picknext is just for fun
      picknext = .true.
      cycle
    end if

    !--------------------!
    !   Create c-point   !
    !--------------------!
    ic = iic-npts1
    icg(ic) = 1
    icg(ifg(iic)) = icg(iic)
    ifg(icg(iic)) = ifg(iic)

    !-----------------------------------------------------------------!
    !   For point ic with eccessive number of strong transpose        !
    !   connections: make it a c-point but let its connected points   !
    !   remain undecided.                                             !
    !-----------------------------------------------------------------!
    if(jcnbhi .ne. jvalx) then

      !------------------------------------------!
      !   Create f-points around above c-point   !
      !------------------------------------------!
      do j = iw(ic+iws), iw(ic+iws+1) - 1
        i = jtr(j)
        if(icg(i) .eq. -2) then
          icg(i) = -1
          ii = i+npts1
          icg(ifg(ii)) = icg(ii)
          ifg(icg(ii)) = ifg(ii)

          !-------------------------------------------------------------!
          !   Increment measure for all strong u-connections iii of i   !
          !   (if not yet measure=jvalx) and put them on reset stack    !
          !   (if not yet there)                                        !
          !-------------------------------------------------------------!
          do jj = ia(i) + 1, ja(ia(i))
            iii = ja(jj)
            if(icg(iii).eq.-2 .and. ifg(iii).lt.jvalx) then
              ifg(iii) = ifg(iii)+1
              iiii = iii-ilo1
              if(ir(iiii) .eq. 0) then
                ir(iiii) = itop
                itop = iiii
              end if
            end if
          end do
        end if
      end do

      !------------------------------------------------------------!
      !   Decrement measure for all strong u-connections i of ic   !
      !   and put them on reset-stack (if not yet there)           !
      !   (Why Ruge-Stueben sometimes put comments at the end of   !
      !   if-endif construct is beyond me.  But I keep it here.)   !
      !------------------------------------------------------------!
    end if

    do j = ia(ic) + 1, ja(ia(ic))
      i = ja(j)
      if(icg(i) .eq. -2) then
        ifg(i) = ifg(i)-1
        ii = i-ilo1
        if(ir(ii) .eq. 0) then
          ir(ii) = itop
          itop = ii
        end if
      end if
    end do

    !----------------------------------------------------------!
    !   Rearrange the lists by sweeping through reset-stack.   !
    !   then go back to pick another u-point ic.               !
    !----------------------------------------------------------!
    in = itop
    itop = -1

    do
      if(in .le. 0) then
        picknext = .true.
        exit
      end if
      i = in+ilo1
      ii = i+npts1
      if(icg(i) .eq. -2) then
        ifg(icg(ii)) = ifg(ii)
        icg(ifg(ii)) = icg(ii)
        jv = ifg(i)
        if (jv.gt.jcnbhi) jcnbhi = jv
        icg(ii) = icg(jv)
        ifg(ii) = jv
        icg(jv) = ii
        ifg(icg(ii)) = ii
      end if
      iip = in
      in = ir(in)
      ir(iip) = 0
    end do
    if(picknext) then
      cycle
    else
      exit
    end if
  end do

  end subroutine
