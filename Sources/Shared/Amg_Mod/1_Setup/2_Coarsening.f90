!==============================================================================!
  subroutine Coarsening(Amg, levelx,   &
                        a, u, ia, ja,  &  ! linear system
                        iw, icg, ifg,  &  ! work arrays
                        levels,        &
                        iwork, jtr)
!------------------------------------------------------------------------------!
!     Performs coarsening, defines interpolations and cgc-operators
!
!     ============== standard values of parameters =====================
!
!             ecg1=0.,    ecg2=0.25,    ewt2=0.35,   nwt=2
!
!     ================ description of parameters =======================
!
!     ecg1 --    defines criterion for determining diagonal
!                dominance in Row_Sort. i.e., if the absolute
!                value of the sum of the off-diagonals of row
!                i is smaller than ecg1 times the absolute
!                value of the diagonal entry, then point i is
!                immediately forced to be an f-point in the
!                pre-coloring algorithm Pre_Color.
!                in the second part (wint), no interpolation
!                is defined for point i, and no points use
!                i for interpolation. in addition, the weight
!                for point i is not distributed to other points
!                when defining the interpolation weights for
!                points which depend on point i. (this is
!                equivalent to completely ignoring such points
!                in determining the coarse grid and interpolation
!                weights.
!
!     ecg2 --    defines strong connections (alpha in the paper)
!
!     ewt2 --    defines strong dependence on a set (beta in the paper)
!
!     nwt  --    parameter controlling the definition of interpolation
!                formulas:
!                  =1 - checking of int-formula: off
!                  =2 - checking of int-formula: on
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[parameters]---------------------------------!
  class(Amg_Type)  :: Amg
  integer          :: levelx
  double precision :: a(:), u(:)
  integer          :: ia(:), ja(:)
  integer          :: iw(:), icg(:), ifg(:)
  integer          :: levels
  integer          :: iwork(:), jtr(:)
!-----------------------------------[locals]-----------------------------------!
  integer          :: i, iajas, icall, ichk, iias, iirs, irst, isaja, isia
  integer          :: jtrst, level, kerr, mdiw, mmax
  integer          :: first_level, ncolx
  integer          :: nda, ndu
  logical          :: exitout, recover
!------------------------------------[save]------------------------------------!
  save  ! this is included only as a precaution as Ruge-Stueben had it
!==============================================================================!

  nda = size(a, 1)
  ndu = size(u, 1)

  !-----------------------------------!
  !   Assign default values if zero   !
  !   (ecg1 and ntr are OK to be 0)   !
  !-----------------------------------!
  if(Amg % nwt  .eq. 0)      Amg % nwt  = 2
  if(Amg % ecg2 .eq. 0.0d0)  Amg % ecg2 = 0.25d0
  if(Amg % ewt2 .eq. 0.0d0)  Amg % ewt2 = 0.35d0

  !--------------------------!
  !   Decode parameter nwt   !
  !--------------------------!
  ichk = Amg % nwt

  !----------------!
  !   Coarsening   !
  !----------------!
  levels = min(levelx, AMG_MAX_LEVELS)
  mmax   = AMG_MAX_LEVELS

  !-----------------------------------------------------------------!
  !   Initialize parameters Amg % mda - mdiw, later set to actual   !
  !   storage requirements of corresponding vectors                 !
  !-----------------------------------------------------------------!
  Amg % mda = 0
  mdiw      = Amg % imax(1)

  !---------------------------------------------------------------------!
  !   first_level is the number of the first stored grid, iajas, iias   !
  !   and iirs are the shifts in vectors a, ja, ia and ir, respecti-    !
  !   vely, as compared to full storage of all grids                    !
  !---------------------------------------------------------------------!
  kerr        = 0
  first_level = 1
  iajas       = 0
  iias        = 0
  iirs        = 0

  do level = 2, levels

    !-------------------------------------------------------------------!
    !   jtrst: initial pointer for work space in vector a, to contain   !
    !          the strong transpose connections                         !
    !-------------------------------------------------------------------!
    do
      exitout = .false.
      do icall = 1, 4
        recover = .false.
        if(icall .eq. 1) then
          jtrst = ia(Amg % imax(level-1)+1)
          call Amg % Row_Sort(level-1,        &
                              a, ia, ja,      &
                              iw, ifg, jtr)
        end if
        if(icall .eq. 2) then

          !-------------------------------------------------------!
          !   irst: initial pointer for work space in vector u,   !
          !         to contain reset stack                        !
          !-------------------------------------------------------!
          irst = mdiw + 1 - iirs
          call Amg % Pre_Color(level-1,       &
                               ia, ja,        &
                               iw, icg, ifg,  &
                               iwork, jtr,    &
                               iias)
        end if
        if(icall .eq. 3) then
          call Amg % Interpolation_Weights(level-1, ichk, mmax,  &
                                           a, ia, ja,            &
                                           iw, ifg, icg, iwork,  &
                                           ncolx, iajas)
        end if
        if(icall .eq. 4) then
          call Amg % Define_Operators(level, mmax,   &
                                      a, ia, ja,     &
                                      iw, icg, ifg,  &
                                      iwork,         &
                                      ncolx, iajas)
        end if
        if(Amg % ierr .gt. 0) then
          ! OK, it will attempt a recovery
          if(Amg % ierr .ge. 1 .and. Amg % ierr .le. 6) then
            if(level .le. first_level + 1 .and. iirs .ne. 0) return
            kerr = Amg % ierr
            Amg % ierr = 0
            first_level = level-1
            iirs = mdiw
            isia = Amg % imin(level-1)-1
            iias = iias+isia
            isaja = ia(Amg % imin(level-1))-1
            iajas = iajas+isaja
            do i = ia(Amg % imin(level-1)), ia(Amg % imax(level-1)+1) - 1
              ja(i-isaja) = ja(i)-isia
              a(i-isaja)  = a(i)
            end do
            do i = Amg % imin(level-1), Amg % imax(level-1)+1
              ia(i-isia) = ia(i)-isaja
            end do
            Amg % imin(level-1) = 1
            Amg % imax(level-1) = Amg % imax(level-1)-isia
            iw(Amg % iminw(level-1)) = ia(Amg % imax(level-1)+1)
            recover = .true.
            exit
          end if
          return
        endif
      end do
      if(.not. recover) then
        if(level .gt. mmax) then
          levels = mmax
          exitout = .true.
          exit
        end if
        call Amg % Truncate_Operator(level, Amg % ntr, a, ia, ja)
        if(Amg % ierr .gt. 0) return

        ! In this AMG implementation, unknowns (u), right-hand side (f), and
        ! the row pointers (ia) are all stored in contiguous blocks, stacked
        ! level after level without gaps. This compact layout does not reserve
        ! space for the final CSR entry, i.e., ia(n+1), which in standard CSR
        ! format marks the end of the last row.
        !
        ! To work around this, the value of ia(imax+1) is temporarily stored
        ! during setup in the integer work array iw (here), at position
        ! iminw(level).  It is then copied back into ia(imax+1) during the
        ! solve phase to safely enable loops of the form:
        !   do j = ia(i), ia(i+1)-1.
        !
        ! While awkard, this workaround is well-contained and avoids
        ! restructuring the data layout.  Attempts to replace it with padding
        ! or offsetting the arrays lead to increased complexity elsewhere, so
        ! it is best left as-is.
        !
        ! For application of this trick (technique), have a look at:
        ! - Calculate_Residual.f90
        ! - Cg_On_Coarsest_Level.f90
        ! - Bicg_On_Coarsest_Level.f90
        ! - Gauss_Seidel_Sweep.f90
        ! - Scale_Solution.f90
        ! - Cg_Alpha.f90
        ! - Cg_Epsilon.f90
        iw(Amg % iminw(level)) = ia(Amg % imax(level)+1)
        if(level .ge. mmax) then
          levels = mmax
          exitout = .true.
          exit
        end if
        ! This will get you out of the inner loop
        exit
      end if
    end do
    if(exitout) exit
  end do

  Amg % mdu = Amg % imax(levels) + iias
  Amg % mdw = mdiw + 2
  if(kerr.ne.0 .or. Amg % mdu.gt.ndu .or. Amg % mdu.gt.ndu) then
    write(6, 1024) Amg % mda, Amg % mda, Amg % mdu,  &
                   Amg % mdu, Amg % mdu, Amg % mdw
    write(6, 2048) mdiw
    Amg % ierr = kerr
    return
  end if

  call Amg % Set_Inverse_Pointer(icg, ifg, levels)

1024  format (/' **************** space requirements ****************'//  &
               ' vector          needed                              '/   &
               ' ----------------------------------------------------'/   &
               '    a ',       i16,   '   adjust the dimensioning of '/   &
               '    ja',       i16,   '   vectors a - ig in the      '/   &
               '    ia',       i16,   '   calling program according  '/   &
               '    u ',       i16,   '   to the calculated space re-'/   &
               '    f ',       i16,   '   quirements and rerun the   '/   &
               '    ig',       i16,   '   program.                   '/   &
               ' ----------------------------------------------------')
2048  format (/' note: if you want to use cg-corrections in the solu-'/   &
               '       tion process (ncyc-subparameter icgr=1 or =2),'/   &
               '       provide for additional',i6,' storage locations'/   &
               '       in vectors u and f.                           '/   &
               '         similarly, usage of the yale-smp solver on  '/   &
               '       the coarsest grid (nsolco=2) will require     '/   &
               '       additional space in vector a during the solu- '/   &
               '       tion phase. in this case, however, its exact  '/   &
               '       amount isn''t predictable.                    '/)

  end subroutine
