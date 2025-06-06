!.......................................................................
!
!     inju                                                   subroutine
!
!.......................................................................
      subroutine inju(kc, u, imin, imax, ifg)
!
!     injects u-values from grid kc-1 to grid kc
!
      double precision u(*)
      integer imin(*),imax(*), ifg(*)
      save
      do ic = imin(kc), imax(kc)
        u(ic) = u(ifg(ic))
      end do

      end
