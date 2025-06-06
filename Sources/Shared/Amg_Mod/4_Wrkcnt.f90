!==============================================================================!
  subroutine wrkcnt(amg, iout,  &
                    ia,         &
                    iw,         &
                    levels, ncyc0)
!-----------------------------------------------------------------------------!
!   Residuals / cp-times / complexity / dimensioning
!-----------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(amg_type)  :: amg
  integer          :: iout
  integer          :: ia(:)
  integer          :: iw(:)
  integer          :: levels
  integer          :: ncyc0
!-----------------------------------[locals]-----------------------------------!
  real             :: t(10), sum1, sum2
  integer          :: i, idima, level, mdta, mdtf, mdtia, mdtig, mdtja,  &
                      mdtu, nnu
  double precision :: acmplx, ocmplx, scmplx, tcmplx, cfac, cfpc
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  if(iout .lt. 1) return

  !-----------------------------!
  !   Residuals / convergence   !
  !-----------------------------!
  if(ncyc0 .gt. 0) then
    cfac = amg % res / (amg % res0+1.d-40)
    write(6, 9110) amg % res0, amg % res, cfac
    cfpc = cfac**(1.d0/dble(ncyc0))
    write(6, 9120) cfpc
  endif
  if(iout .le. 1) return
  write(6, 9000)
  nnu = amg % imax(1)

  !---------------------!
  !   Computing times   !
  !---------------------!
  sum1 = 0.0
  sum2 = 0.0
  do i = 1, 10
    t(i) = 0.0
    if(ncyc0 .gt. 0) t(i) = amg % time(i+10)/real(ncyc0)
    sum1 = sum1 + amg % time(i)
    sum2 = sum2+t(i)
  end do

  write(6, 9020)
  write(6, 9030)
  write(6, 9040) (amg % time(i), t(i), i = 1, 8)
  write(6, 9030)
  write(6, 9050) sum1,sum2
  write(6, 9030)

  !--------------------------------------------------!
  !   Space occupied by operators a(1) - a(levels)   !
  !--------------------------------------------------!
  idima = 0
  do level = 1, levels
    idima = idima+iw(amg % iminw(level))-ia(amg % imin(level))
  end do

  !--------------------------------------------!
  !   Theoretical minimal space requirements   !
  !--------------------------------------------!
  if(levels .lt. 2) return
  mdta  = iw(amg % iminw(levels))-1
  mdtja = mdta
  mdtia = amg % imax(levels)+1
  mdtu  = amg % imax(levels)
  mdtf  = mdtu
  mdtig = 2*mdtu+nnu
  write(6, 9100) amg % mda, mdta, amg % mda, mdtja,  &
                   amg % mdu, mdtia, amg % mdu, mdtu,  &
                   amg % mdu, mdtf,  amg % mdw, mdtig

  !------------------!
  !   Complexities   !
  !------------------!
  scmplx = dble(2*(  amg % mda + amg % mdu + amg % mdu)   &
                   + amg % mda + amg % mdu + amg % mdw)/  &
           dble(1+5*nnu+3*(iw(amg % iminw(1))-ia(amg % imin(1))))
  tcmplx = dble(2*(mdta+mdtu+mdtf)+mdtja+mdtia+mdtig)/  &
           dble(1+5*nnu+3*(iw(amg % iminw(1))-ia(amg % imin(1))))

  acmplx = dble(idima)/dble(iw(1)-1)
  ocmplx = dble(mdtu )/dble(nnu)
  write(6, 9080) acmplx,ocmplx,scmplx,tcmplx

9000  format (//' ************** work count ***************'/)
9020  format (  '   prep       sec       sol      sec/cycle')
9030  format (  ' -----------------------------------------')
9040  format (  ' 1 row_sort',f7.2,'   11 intadd   ',f7.2/  &
                ' 2 pre-col ',f7.2,'   12 rescal   ',f7.2/  &
                ' 3 chk-col ',f7.2,'   13 relax    ',f7.2/  &
                ' 4 interpol',f7.2,'   14 v-*      ',f7.2/  &
                ' 5 restrict',f7.2,'   15 others   ',f7.2/  &
                ' 6 opdfn   ',f7.2,'   16 conj-grad',f7.2/  &
                ' 7 trunc   ',f7.2,'   17 yale-smp ',f7.2/  &
                ' 8 others  ',f7.2,'   18 ------   ',f7.2)
9050  format (  '   sum     ',f7.2,'      sum      ',f7.2)
9080  format (/' ******************* complexities ********************'/  &
               ' space occupied by all operators / space of operator  '/  &
               ' on the finest grid   = ',f8.2,'   (a-complexity)     '/  &
               ' total number of grid points / number of points in    '/  &
               ' the  finest  grid    = ',f8.2,'   (o-complexity)     '/  &
               ' total space used by amg1r5 / space occupied by user- '/  &
               ' defined  problem     = ',f8.2,'   (s-complexity)     '/  &
               ' space used during solution phase / space occupied by '/  &
               ' user-defined problem = ',f8.2                         /  &
               ' *****************************************************')
9100  format (//' ********* space requirements *********'//  &
                ' vector      needed      theor. minimum'/   &
                ' --------------------------------------'/   &
                '   a ',        i14   ,2x,     i16       /   &
                '   ja',        i14   ,2x,     i16       /   &
                '   ia',        i14   ,2x,     i16       /   &
                '   u ',        i14   ,2x,     i16       /   &
                '   f ',        i14   ,2x,     i16       /   &
                '   ig',        i14   ,2x,     i16       /   &
                ' --------------------------------------'/)
9110  format (/' **************** convergence *****************'/  &
               ' l2-norm of residual before cycling =',1pe10.3/    &
               ' l2-norm of residual after  cycling =',1pe10.3/    &
               ' convergence factor                 =',1pe10.3)
9120  format ( ' convergence factor per cycle       =',1pe10.3)

  end subroutine
