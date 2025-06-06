!.......................................................................
!
!     injf                                                   subroutine
!
!.......................................................................
      subroutine injf(kc,f,imin,imax,ifg)
!
!     injects f-values from grid kc-1 to grid kc
!
      double precision f(*)
      integer imin(*),imax(*),ifg(*)
      save
      do ic = imin(kc), imax(kc)
        f(ic) = f(ifg(ic))
      end do

      end
