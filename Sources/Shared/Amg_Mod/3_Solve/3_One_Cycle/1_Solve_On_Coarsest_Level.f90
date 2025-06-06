!==============================================================================!
  subroutine Solve_On_Coarsest_Level(Amg, m, ifac,     &
                                     a, u, f, ia, ja,  &
                                     iw, icg)
!------------------------------------------------------------------------------!
!   Solves on coarsest grid, either with Gauss-Seidel relaxation (nsc=1)
!   or with BiCG solver (nsc=3).  Yale8 (nsc=2) has been decomissioned.
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type)  :: Amg
  integer          :: nifac
  double precision :: a(:), u(:), f(:)
  integer          :: ia(:), ja(:)
  integer          :: iw(:), icg(:)
!-----------------------------------[locals]-----------------------------------!
  double precision :: fmax, resnew, resold
  integer          :: m,ifac,aaux,esp,flag,i,iaux,ihi,ii,ilo,is, &
                      iter,j,jhi,jj,jlo,jpos,js,np,npoint,nsp,path

  !---------------------------------------------------------------------!
  !   conv: if coarse grid solution is done with gs-relaxation and      !
  !   nrcx=0, as many gs-sweeps are performed as are necessary to re-   !
  !   duce the residual by the factor conv                              !
  !---------------------------------------------------------------------!
  double precision, parameter :: conv =1.d-2
  integer,          parameter :: gauss=1
  integer,          parameter :: yale =2
  integer,          parameter :: bicg =3
!#  logical                     :: yalefail
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

!#  !------------------------------------------!
!#  !                                          !
!#  !   Attempt a solution with Yale8 solver   !
!#  !                                          !
!#  !------------------------------------------!
!#  if(Amg % nsc .eq. yale) then
!#
!#    !------------------------------------!
!#    !   Attempt solution with yale-smp   !
!#    !------------------------------------!
!#    yalefail = .false.
!#    call Amg % timer_start()
!#    ilo = Amg % imin(m)
!#    jlo = ia(ilo)
!#
!#    !--------------------------------------------------!
!#    !   First call on grid m, first factorize matrix   !
!#    !--------------------------------------------------!
!#    if (ifac.eq.1) then
!#      ihi = Amg % imax(m)
!#      jhi = iw(iminw(m))-1
!#      np = ihi-ilo+1
!#      is = ilo-1
!#      js = jlo-1
!#
!#      !----------------------------------!
!#      !   Test of available work space   !
!#      !----------------------------------!
!#      if (jhi+3*np.gt.nda) then
!#        Amg % nsc = 1
!#        write(6, '(a,a)')                                        &
!#          ' --- warng in coarse: no yale-smp because nda too ',  &
!#          'small'
!#        Amg % ierr = AMG_WARN_YALE_STORAGE_JA
!#        yalefail = .true.
!#      end if
!#
!#      !------------------------------------------------!
!#      !   Initialisation of yale-smp pointer vectors   !
!#      !------------------------------------------------!
!#      if(.not. yalefail) then
!#        do i = 1, np
!#          ja(jhi+i) = i
!#          ja(jhi+np+i) = i
!#          ja(jhi+2*np+i)= i
!#        end do
!#        if (Amg % irow0.ne.1) then
!#
!#          !--------------------------------------------------!
!#          !   Coarse grid operator regular, shift contents   !
!#          !   of pointer vectors ia and ja                   !
!#          !--------------------------------------------------!
!#          do i = ilo, ihi
!#            ia(i) = ia(i)-js
!#          end do
!#          iaux = ia(ihi+1)
!#          ia(ihi+1) = iw(iminw(m))-js
!#          do i = jlo, jhi
!#            ja(i) = ja(i)-is
!#          end do
!#          npoint = np
!#        else
!#
!#          !-------------------------------------------------------------------!
!#          !   Coarse grid operator has rowsum zero, eliminate last solution   !
!#          !   component by setting u(ihi) to zero and cancelling the cor-     !
!#          !   responding entries in a and ja, respectively. The cancelled     !
!#          !   entries are stored before the last row in a and ja.             !
!#          !-------------------------------------------------------------------!
!#          iaux = ia(ihi)
!#
!#          !----------------------------------------------!
!#          !   jpos: pointer to position in a and ja to   !
!#          !         contain next eliminated entry        !
!#          !----------------------------------------------!
!#          jpos = iaux-1
!#          j = ia(ilo)
!#          do i = ilo, ihi-1
!#            ia(i) = j-js
!#            ja(j) = i-is
!#            j = j+1
!#            do
!#              if(j .eq. ia(i+1)) exit
!#              if(ja(j) .ne. ihi) then
!#                ja(j) = ja(j)-is
!#                j = j+1
!#              else
!#                aaux = a(j)
!#                do jj = j, jpos-1
!#                  a(jj) = a(jj+1)
!#                  ja(jj) = ja(jj+1)
!#                end do
!#                do ii = i+1, ihi
!#                  ia(ii) = ia(ii)-1
!#                end do
!#                a(jpos) = aaux
!#                ja(jpos) = i
!#                jpos = jpos-1
!#              end if
!#            end do
!#          end do
!#          ia(ihi) = ia(ihi)-js
!#
!#          !------------------------------------------!
!#          !   Decrease number of points by one and   !
!#          !   set last solution component to zero    !
!#          !------------------------------------------!
!#          npoint = np-1
!#          u(ihi) = 0.d0
!#        endif
!#        nsp = nda-jhi
!#        path = 1
!#
!#        !---------------------------------------------------!
!#        !  This line throws an error on modern compilers:   !
!#        !  Error: Type mismatch in argument ‘isp’ at (1);   !
!#        !  passed REAL(8) to INTEGER(4)                     !
!#        !---------------------------------------------------!
!#        call ndrv(npoint,ja(jhi+1),ja(jhi+np+1),ja(jhi+2*np+1),  &
!#                  ia(ilo),ja(jlo),a(jlo),f(ilo),u(ilo),nsp,      &
!#                  a(jhi+1),a(jhi+1),esp,path,flag)
!#
!#        !-----------------------------------------------------------------!
!#        !   Restore previous values for ia, ja and a, in particular put   !
!#        !   back eliminated matrix entries, if rowsum zero                !
!#        !-----------------------------------------------------------------!
!#        do i = ilo, ihi
!#          ia(i) = ia(i)+js
!#        end do
!#        if (Amg % irow0.ne.1) then
!#          ia(ihi+1) = iaux
!#          do i = jlo, jhi
!#            ja(i) = ja(i)+is
!#          end do
!#        else
!#          do i = jlo, jpos
!#            ja(i)=ja(i)+is
!#          end do
!#          do j = jpos+1, iaux-1
!#            aaux = a(j)
!#            i = ja(j)
!#            do ii = i+1, ihi
!#              ia(ii) = ia(ii)+1
!#            end do
!#            do jj = j, ia(i+1), -1
!#              a(jj) = a(jj-1)
!#              ja(jj) = ja(jj-1)
!#            end do
!#            a(ia(i+1)-1) = aaux
!#            ja(ia(i+1)-1) = ihi
!#          end do
!#        endif
!#
!#        !----------------------------------------------!
!#        !   If an error occured during execution of    !
!#        !   ndrv, solve with gauss-seidel relaxation   !
!#        !----------------------------------------------!
!#        if(flag.ne.0) then
!#          Amg % nsc = gauss
!#          if (esp.lt.0) then
!#            write(6, '(a,a)')                                        &
!#              ' --- warng in coarse: no yale-smp because nda too ',  &
!#              'small'
!#            Amg % ierr = AMG_WARN_YALE_STORAGE_A
!#          else
!#            write(6, '(a,a,i8)')                                     &
!#              ' --- warng in coarse: no yale-smp because of error',  &
!#              ' in factorization, code=',                            &
!#              flag
!#            Amg % ierr = AMG_ERR_YALE_FACTOR_FAIL
!#          endif
!#          yalefail = .true.
!#        else
!#
!#          !--------------------------------------------------------------!
!#          !   Factorization successfull, update limits of used storage   !
!#          !--------------------------------------------------------------!
!#          mda = max(mda,nda-esp)
!#          ifac = 0
!#        endif
!#      end if
!#
!#    !--------------------------------------------!
!#    !   Not the first call, factorization done   !
!#    !--------------------------------------------!
!#    else
!#      if(.not. yalefail) then
!#
!#        !---------------------------------!
!#        !   Factorization allready done   !
!#        !---------------------------------!
!#        if (Amg % irow0.ne.1) then
!#          npoint = np
!#        else
!#          npoint = np-1
!#          u(ihi) = 0.d0
!#        endif
!#        path = 3
!#
!#        !---------------------------------------------------!
!#        !  This line throws an error on modern compilers:   !
!#        !  Error: Type mismatch in argument ‘isp’ at (1);   !
!#        !  passed REAL(8) to INTEGER(4)                     !
!#        !---------------------------------------------------!
!#        call ndrv(npoint,ja(jhi+1),ja(jhi+np+1),ja(jhi+2*np+1),  &
!#                  ia(ilo),ja(jlo),a(jlo),f(ilo),u(ilo),nsp,      &
!#                  a(jhi+1),a(jhi+1),esp,path,flag)
!#
!#        !---------------------------!
!#        !   The end yale attempts   !
!#        !---------------------------!
!#      end if
!#    end if
!#
!#    !-------------------------!
!#    !   Update time counter   !
!#    !-------------------------!
!#    if(.not. yalefail) then
!#      call Amg % timer_stop(17)
!#      return
!#    end if
!#  end if

  !-------------------------------------------!
  !                                           !
  !   Solution with gauss-seidel relaxation   !
  !                                           !
  !-------------------------------------------!
  if(Amg % nsc .eq. gauss) then

    !---------------------------------------------!
    !   Just perform given number of iterations   !
    !---------------------------------------------!
    if(Amg % nrcx .ne. 0) then
      do iter = 1, Amg % nrcx
        call Amg % Gauss_Seidel_Sweep(m, 2,             &
                                      a, u, f, ia, ja,  &
                                      iw, icg)
      end do

    !-----------------------------------------------------!
    !   Reduce residual on coarsest grid by factor conv   !
    !   if not yet in the range of the truncation error   !
    !-----------------------------------------------------!
    else

      call Amg % timer_start()

     ! Calculate supremum norm of right hand side
      fmax = 0.d0
      do i = Amg % imin(m), Amg % imax(m)
        fmax = max(fmax,abs(f(i)))
      end do
      call Amg % Calculate_Residual(m, resold,  &
                                    a, u, f, ia, ja,  &
                                    iw)

      call Amg % timer_stop(15)
      resold = max(resold*conv,fmax*1.d-12)
      do i = 1, 10
        do j = 1, 10
          call Amg % Gauss_Seidel_Sweep(m, 2, a, u, f, ia, ja,  &
                                        iw, icg)
        end do
        call Amg % timer_start()
        call Amg % Calculate_Residual(m, resnew, a, u, f, ia, ja,  &
                                      iw)

        call Amg % timer_stop(15)
        if(resnew .le. resold) return
      end do
    end if

  !-------------------------------!
  !                               !
  !   Solution with BiCG solver   !
  !                               !
  !-------------------------------!
  else if(Amg % nsc .eq. bicg) then

    call Amg % Cg_On_Coarsest_Level(m, 2,             &
                                    a, u, f, ia, ja,  &
                                    iw, icg)
    call Amg % Bicg_On_Coarsest_Level(m, 2,             &
                                      a, u, f, ia, ja,  &
                                      iw,icg)

  end if

  end subroutine
