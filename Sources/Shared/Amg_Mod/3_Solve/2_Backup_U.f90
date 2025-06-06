!==============================================================================!
  subroutine backup_u(amg, level, icgr, u, m)
!------------------------------------------------------------------------------!
!   Makes a back-up of the current approx. on level "level" if icgr.ne.0.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(amg_type)  :: amg
  integer          :: level, icgr
  double precision :: u(:)
  integer          :: m
!-----------------------------------[locals]-----------------------------------!
  integer :: i, ishift
  integer :: ndu
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  ndu = size(u, 1)

  if(icgr .eq. 0) return

  ishift = amg % imax(m)+1-amg % imin(level)
  amg % mdu = max(amg % mdu, amg % imax(level)+ishift)

  if(amg % mdu .gt. ndu .or. amg % mdu .gt. ndu) then
    write (6, '(a,a,i9,a,i9)')                              &
      ' --- warng in usave: no cg because of storage ---',  &
      '     required: ndu =',  amg % mdu, ' ndu = ', amg % mdu
    amg % ierr = AMG_WARN_CG_STORAGE_U
    icgr = 0
    return
  end if

  call amg % timer_start()

  do i = amg % imin(level), amg % imax(level)
    u(i+ishift) = u(i)
  end do

  call amg % timer_stop(16)

  end subroutine
