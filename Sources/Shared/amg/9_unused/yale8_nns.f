!       ----------------------------------------------------------------
!
!*** subroutine nns
!*** numeric solution of a sparse nonsymmetric system of linear
!      equations given ldu-factorization (uncompressed pointer storage)
!
        subroutine  nns8
     *     (n, r,c, il,jl,l, d, iu,ju,u, z, b, tmp)
!
!       input variables:   n, r,c, il,jl,l, d, iu,ju,u, b
!       output variables:  z
!
!       parameters used internally:
! fia     tmp   - holds new right-hand side b' for solution of the
!                   equation ux = b'.
!                   size = n.
!
        integer  r(*), c(*),  il(*), jl(*),  iu(*), ju(*)
!       real  l(*), d(*), u(*),  z(*), b(*),  tmp(*), sum
        double precision  l(*), d(*), u(*),  z(*), b(*),  tmp(*), sum
!
!  ******  solve ldy = b  by forward substitution  *********************
        do k = 1, n
          sum = b(r(k))
          jmin = il(k)
          jmax = il(k+1) - 1
          if(jmin .le. jmax) then
            do j = jmin, jmax
              sum = sum - l(j) * tmp(jl(j))
            end do
          endif
          tmp(k) = sum * d(k)
        end do
!
!  ******  solve  ux = y  by back substitution  ************************
        k = n
        do i=1,n
          sum = tmp(k)
          jmin = iu(k)
          jmax = iu(k+1) - 1
          if(jmin .le. jmax) then
            do j = jmin, jmax
              sum = sum - u(j) * tmp(ju(j))
            end do
          endif
          tmp(k) = sum
          z(c(k)) = sum
          k = k-1
        end do

        end
