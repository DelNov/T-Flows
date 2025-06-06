!       ----------------------------------------------------------------
!
!*** subroutine nnf
!*** numeric ldu-factorization of sparse nonsymmetric matrix and
!      solution of system of linear equations (uncompressed pointer
!      storage)
!
        subroutine  nnf8
     *     (n, r,c,ic, ia,ja,a, z, b, il,jl,l,lmax, d, iu,ju,u,umax,
     *      row, tmp, flag)
!
!       input variables:   n, r,c,ic, ia,ja,a, b, il,jl,lmax, iu,ju,umax
!       output variables:  z, l,d,u, flag
!
!       parameters used internally:
! fia     row   - holds intermediate values in calculation of l, d, u.
!                   size = n.
! fia     tmp   - holds new right-hand side b' for solution of the
!                   equation  ux = b'.
!                   size = n.
!
        integer  r(*), c(*), ic(*),  ia(*), ja(*),
     *     il(*), jl(*), lmax,  iu(*), ju(*), umax,  flag
!       real  a(*), z(*), b(*),  l(*), d(*), u(*),
!    *     row(*), tmp(*),  li, sum, dk
        double precision  a(*), z(*), b(*),  l(*), d(*), u(*),
     *     row(*), tmp(*),  li, sum, dk
!
!  ******  check storage  **********************************************
! ** error:  insufficient storage for l
        if(il(n+1)-1 .gt. lmax) then
          flag = 4*n + 1
          return
        end if
! ** error:  insufficient storage for u
        if(iu(n+1)-1 .gt. umax) then
          flag = 7*n + 1
          return
        end if
!
!  ******  for each row  ***********************************************
        do k = 1, n
!  ******  set the initial structure of row  ***************************
          jmin = il(k)
          jmax = il(k+1) - 1
!  ******  if l(k,m) .ne. 0, row(m)=0  *********************************
          do j = jmin, jmax
            row(jl(j)) = 0
          end do
          row(k) = 0
          jmin = iu(k)
          jmax = iu(k+1) - 1
!  ******  if u(k,m) .ne. 0, row(m)=0  *********************************
          do j = jmin, jmax
            row(ju(j)) = 0
          end do
          jmin = ia(r(k))
          jmax = ia(r(k)+1) - 1
!  ******  set row to kth row of reordered a  **************************
          do j = jmin, jmax
            row(ic(ja(j))) = a(j)
          end do
!  ******  initialize sum  *********************************************
          sum = b(r(k))
!
!  ******  assign the kth row of l and adjust row, sum  ****************
          imin = il(k)
          imax = il(k+1) - 1
          do i = imin, imax
            li = - row(jl(i))
!  ******  if l is not required, then comment out the following line  **
            l(i) = - li
            sum = sum + li * tmp(jl(i))
            jmin = iu(jl(i))
            jmax = iu(jl(i)+1) - 1
            do j = jmin, jmax
              row(ju(j)) = row(ju(j)) + li * u(j)
            end do
          end do
!
!  ******  assign diagonal d and kth row of u, set tmp(k)  *************
! ** error:  zero pivot
          if(row(k) .eq .0) then
            flag = 8*n + k
            return
          end if
          dk = 1 / row(k)
          d(k) = dk
          tmp(k) = sum * dk
          jmin = iu(k)
          jmax = iu(k+1) - 1
          do j = jmin, jmax
            u(j) = row(ju(j)) * dk
          end do
        end do
!
!  ******  solve  ux = tmp  by back substitution  **********************
        k = n
        do i = 1, n
          sum = tmp(k)
          jmin = iu(k)
          jmax = iu(k+1) - 1
          do j = jmin, jmax
            sum = sum - u(j) * tmp(ju(j))
          end do
          tmp(k) = sum
          z(c(k)) = sum
          k = k-1
        end do

        flag = 0

        end
