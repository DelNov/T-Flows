!==============================================================================!
  subroutine Backup_U(Amg, level, icgr, u, m)
!------------------------------------------------------------------------------!
!   Makes a back-up of the current approx. on level "level" if icgr.ne.0.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type) :: Amg
  integer         :: level, icgr
  real            :: u(:)
  integer         :: m
!-----------------------------------[locals]-----------------------------------!
  integer :: i, ishift
  integer :: ndu
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  ndu = size(u, 1)

  if(icgr .eq. 0) return

  ishift = Amg % imax(m)+1-Amg % imin(level)
  Amg % mdu = max(Amg % mdu, Amg % imax(level)+ishift)

  if(Amg % mdu .gt. ndu .or. Amg % mdu .gt. ndu) then
    write (6, '(a,a,i9,a,i9)')                              &
      ' --- warng in usave: no cg because of storage ---',  &
      '     required: ndu =',  Amg % mdu, ' ndu = ', Amg % mdu
    Amg % ierr = AMG_WARN_CG_STORAGE_U
    icgr = 0
    return
  end if

  do i = Amg % imin(level), Amg % imax(level)
    u(i+ishift) = u(i)
  end do

  end subroutine
