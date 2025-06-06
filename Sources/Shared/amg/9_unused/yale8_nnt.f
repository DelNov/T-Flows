!       ----------------------------------------------------------------
!
!*** subroutine nnt
!*** numeric solution of the transpose of a sparse nonsymmetric system
!      of linear equations given ldu-factorization (uncompressed pointer
!      storage)
!
        subroutine  nnt8
     *     (n, r,c, il,jl,l, d, iu,ju,u, z, b, tmp)
!
!       input variables:   n, r,c, il,jl,l, d, iu,ju,u, b
!       output variables:  z
!
!       parameters used internally:
! fia     tmp   - holds new right-hand side b' for solution of the
!                   equation lx = b'.
!                   size = n.
!
        integer  r(*), c(*),  il(*), jl(*),  iu(*), ju(*)
!       real  l(*), d(*), u(*),  z(*), b(*),  tmp(*), tmpk
        double precision  l(*), d(*), u(*),  z(*), b(*),  tmp(*), tmpk
!
!  ******  solve  ut y = b  by forward substitution  *******************
        do k = 1, n
          tmp(k) = b(c(k))
        end do
        do k = 1, n
          tmpk = - tmp(k)
          jmin = iu(k)
          jmax = iu(k+1) - 1
          do j = jmin, jmax
            tmp(ju(j)) = tmp(ju(j)) + u(j) * tmpk
          end do
        end do
!
!  ******  solve  d lt x = y  by back substitution  ********************
        k = n
        do i = 1, n
          tmpk = - tmp(k) * d(k)
          jmin = il(k)
          jmax = il(k+1) - 1
          do j = jmin, jmax
            tmp(jl(j)) = tmp(jl(j)) + l(j) * tmpk
          end do
          z(r(k)) = -tmpk
          k = k-1
        end do

        end
