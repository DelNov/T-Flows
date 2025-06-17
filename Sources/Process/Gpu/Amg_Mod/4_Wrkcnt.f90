!==============================================================================!
  subroutine Wrkcnt(Amg,      &
                    ia,       &
                    iw,       &
                    levels)
!-----------------------------------------------------------------------------!
!   Residuals / cp-times / complexity / dimensioning
!-----------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type) :: Amg
  integer         :: ia(:)
  integer         :: iw(:)
  integer         :: levels
!-----------------------------------[locals]-----------------------------------!
  real(SP) :: t(10), sum1, sum2
  integer  :: i, idima, level, mdta, mdtf, mdtia, mdtig, mdtja, mdtu, nnu
  real     :: acmplx, ocmplx, scmplx, tcmplx, cfac, cfpc
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  if(Amg % iout .lt. 1) return

  !-----------------------------!
  !   Residuals / convergence   !
  !-----------------------------!
  if(Amg % ncyc0 .gt. 0) then
    cfac = Amg % res / (Amg % res0+1.0e-40)
    write(6, 9110) Amg % res0, Amg % res, cfac
    cfpc = cfac**(1.0/real(Amg % ncyc0))
    write(6, 9120) cfpc
  end if

  if(Amg % iout .le. 1) return

  nnu = Amg % imax(1)

  !--------------------------------------------------!
  !   Space occupied by operators a(1) - a(levels)   !
  !--------------------------------------------------!
  idima = 0
  do level = 1, levels
    idima = idima+iw(Amg % iminw(level))-ia(Amg % imin(level))
  end do

  !--------------------------------------------!
  !   Theoretical minimal space requirements   !
  !--------------------------------------------!
  if(levels .lt. 2) return
  mdta  = iw(Amg % iminw(levels))-1
  mdtja = mdta
  mdtia = Amg % imax(levels)+1
  mdtu  = Amg % imax(levels)
  mdtf  = mdtu
  mdtig = 2*mdtu+nnu
  write(6, 9100) Amg % mda, mdta, Amg % mda, mdtja,  &
                 Amg % mdu, mdtia, Amg % mdu, mdtu,  &
                 Amg % mdu, mdtf,  Amg % mdw, mdtig

  !------------------!
  !   Complexities   !
  !------------------!
  scmplx = real(2*(  Amg % mda + Amg % mdu + Amg % mdu)   &
                   + Amg % mda + Amg % mdu + Amg % mdw)/  &
           real(1+5*nnu+3*(iw(Amg % iminw(1))-ia(Amg % imin(1))))
  tcmplx = real(2*(mdta+mdtu+mdtf)+mdtja+mdtia+mdtig)/  &
           real(1+5*nnu+3*(iw(Amg % iminw(1))-ia(Amg % imin(1))))

  acmplx = real(idima)/real(iw(1)-1)
  ocmplx = real(mdtu )/real(nnu)
  write(6, 9080) acmplx,ocmplx,scmplx,tcmplx

9080  format (/' ******************* complexities ********************'/  &
               ' space occupied by all operators / space of operator  '/  &
               ' on the finest grid   = ',f8.2,'   (a-complexity)     '/  &
               ' total number of grid points / number of points in    '/  &
               ' the  finest  grid    = ',f8.2,'   (o-complexity)     '/  &
               ' total space used by Amg1r5 / space occupied by user- '/  &
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
