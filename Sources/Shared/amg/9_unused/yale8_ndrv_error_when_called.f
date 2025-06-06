!***********************************************************************
! $large
! $nofloatcalls
!                            appendix 1                          1/15/81
!
!        subroutines for solving sparse nonsymmetric systems
!        of linear equations  (uncompressed pointer storage)
!
!        real*8 version. note: the original subroutines
!
!            ndrv, nsf, nnf, nns and nnt
!
!        have been renamed to
!
!            yale8, nsf8, nnf8, nns8 and nnt8, respectively.
!
!*** subroutine yale8 (old name: ndrv)
!*** driver for subroutines for solving sparse nonsymmetric systems of
!       linear equations (uncompressed pointer storage)
!
!       subroutine  ndrv  (= old name)
!       subroutine  yale8 (= new name)
        subroutine  ndrv(n,    !  1
     *                   r,    !  2
     *                   c,    !  3
     *                   ic,   !  4
     *                   ia,   !  5
     *                   ja,   !  6
     *                   a,    !  7
     *                   b,    !  8
     *                   z,    !  9
     *                   nsp,  ! 10
! Error appears here: we send a real array, but it expects integer array
     *                   isp,  ! 11
     *                   rsp,
     *                   esp,
     *                   path,
     *                   flag)
!
!    parameters
!    class abbreviations are --
!       n - integer variable
!       f - real variable
!       v - supplies a value to the driver
!       r - returns a result from the driver
!       i - used internally by the driver
!       a - array
!
! class   parameter
! ------+----------
!
!         the nonzero entries of the coefficient matrix m are stored
!    row-by-row in the array a.  to identify the individual nonzero
!    entries in each row, we need to know in which column each entry
!    lies.  the column indices which correspond to the nonzero entries
!    of m are stored in the array ja;  i.e., if  a(k) = m(i,j),  then
!    ja(k) = j.  in addition, we need to know where each row starts and
!    how long it is.  the index positions in ja and a where the rows of
!    m begin are stored in the array ia;  i.e., if m(i,j) is the first
!    nonzero entry (stored) in the i-th row and a(k) = m(i,j),  then
!    ia(i) = k.  moreover, the index in ja and a of the first location
!    following the last element in the last row is stored in ia(n+1).
!    thus, the number of entries in the i-th row is given by
!    ia(i+1) - ia(i),  the nonzero entries of the i-th row are stored
!    consecutively in
!            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
!    and the corresponding column indices are stored consecutively in
!            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
!    for example, the 5 by 5 matrix
!                ( 1. 0. 2. 0. 0.)
!                ( 0. 3. 0. 0. 0.)
!            m = ( 0. 4. 5. 6. 0.)
!                ( 0. 0. 0. 7. 0.)
!                ( 0. 0. 0. 8. 9.)
!    would be stored as
!                 1  2  3  4  5  6  7  8  9
!            ---+--------------------------
!            ia   1  3  4  7  8 10
!            ja   1  3  2  2  3  4  4  4  5
!             a   1. 2. 3. 4. 5. 6. 7. 8. 9.         .
!
! nv      n     - number of variables/equations.
! fva     a     - nonzero entries of the coefficient matrix m, stored
!                   by rows.
!                   size = number of nonzero entries in m.
! nva     ia    - pointers to delimit the rows in a.
!                   size = n+1.
! nva     ja    - column numbers corresponding to the elements of a.
!                   size = size of a.
! fva     b     - right-hand side b;  b and z can the same array.
!                   size = n.
! fra     z     - solution x;  b and z can be the same array.
!                   size = n.
!
!         the rows and columns of the original matrix m can be
!    reordered (e.g., to reduce fillin or ensure numerical stability)
!    before calling the driver.  if no reordering is done, then set
!    r(i) = c(i) = ic(i) = i  for i=1,...,n.  the solution z is returned
!    in the original order.
!
! nva     r     - ordering of the rows of m.
!                   size = n.
! nva     c     - ordering of the columns of m.
!                   size = n.
! nva     ic    - inverse of the ordering of the columns of m;  i.e.,
!                   ic(c(i)) = i  for i=1,...,n.
!                   size = n.
!
!         the solution of the system of linear equations is divided into
!    three stages --
!      nsf -- the matrix m is processed symbolically to determine where
!              fillin will occur during the numeric factorization.
!      nnf -- the matrix m is factored numerically into the product ldu
!              of a unit lower triangular matrix l, a diagonal matrix d,
!              and a unit upper triangular matrix u, and the system
!              mx = b  is solved.
!      nns -- the linear system  mx = b  is solved using the ldu
!  or          factorization from nnf.
!      nnt -- the transposed linear system  mt x = b  is solved using
!              the ldu factorization from nnf.
!    for several systems whose coefficient matrices have the same
!    nonzero structure, nsf need be done only once (for the first
!    system);  then nnf is done once for each additional system.  for
!    several systems with the same coefficient matrix, nsf and nnf need
!    be done only once (for the first system);  then nns or nnt is done
!    once for each additional right-hand side.
!
! nv      path  - path specification;  values and their meanings are --
!                   1  perform nsf and nnf.
!                   2  perform nnf only  (nsf is assumed to have been
!                       done in a manner compatible with the storage
!                       allocation used in the driver).
!                   3  perform nns only  (nsf and nnf are assumed to
!                       have been done in a manner compatible with the
!                       storage allocation used in the driver).
!                   4  perform nnt only  (nsf and nnf are assumed to
!                       have been done in a manner compatible with the
!                       storage allocation used in the driver).
!                   5  perform nsf only.
!
!         various errors are detected by the driver and the individual
!    subroutines.
!
! nr      flag  - error flag;  values and their meanings are --
!                     0     no errors detected
!                     n+k   null row in a  --  row = k
!                    2n+k   duplicate entry in a  --  row = k
!                    3n+k   insufficient storage in nsf  --  row = k
!                    4n+1   insufficient storage in nnf
!                    5n+k   null pivot  --  row = k
!                    6n+k   insufficient storage in nsf  --  row = k
!                    7n+1   insufficient storage in nnf
!                    8n+k   zero pivot  --  row = k
!                   10n+1   insufficient storage in ndrv
!                   11n+1   illegal path specification
!
!         working storage is needed for the factored form of the matrix
!    m plus various temporary vectors.  the arrays isp and rsp should be
!    equivalenced;  integer storage is allocated from the beginning of
!    isp and real storage from the end of rsp.
!
! nv      nsp   - declared dimension of rsp;  nsp generally must
!                   be larger than  5n+3 + 2k  (where  k = (number of
!                   nonzero entries in m)).
! nvira   isp   - integer working storage divided up into various arrays
!                   needed by the subroutines;  isp and rsp should be
!                   equivalenced.
!                   size = lratio*nsp
! fvira   rsp   - real working storage divided up into various arrays
!                   needed by the subroutines;  isp and rsp should be
!                   equivalenced.
!                   size = nsp.
! nr      esp   - if sufficient storage was available to perform the
!                   symbolic factorization (nsf), then esp is set to the
!                   amount of excess storage provided (negative if
!                   insufficient storage was available to perform the
!                   numeric factorization (nnf)).
!
!
!  conversion to double precision
!
!    to convert these routines for double precision arrays, simply use
!    the double precision declarations in place of the real declarations
!    in each subprogram;  in addition, the data value of the integer
!    variable lratio must be set as indicated in subroutine ndrv
!
        integer  r(1), c(1), ic(1),  ia(1), ja(1),  isp(1), esp,
     *     path, flag,  q, im, d, u, row, tmp,  umax
!       real  a(1),  b(1),  z(1),  rsp(1)
        double precision  a(1),  b(1),  z(1),  rsp(1)
!
!  set lratio equal to the ratio between the length of floating point
!  and integer array data;  e. g., lratio = 1 for (real, integer),
!  lratio = 2 for (double precision, integer)
!
        data lratio/2/
!
        if(path.lt.1 .or. 5.lt.path) then
! ** error:  illegal path specification
          flag = 11*n + 1
          return
        end if
!  ******  initialize and divide up temporary storage  *****************
        il = 1
        iu = il + n+1
        jl = iu + n+1
!
!  ******  call nsf if flag is set  ************************************
        if((path-1) * (path-5) .eq. 0) then
          max = (lratio*nsp + 1 - jl) - (n+1) - n
          jlmax = max/2
          q     = jl  + jlmax
          im    = q   + (n+1)
          jutmp = im  +   n
          jumax = lratio*nsp + 1 - jutmp
          esp = max/lratio
! ** error:  insufficient storage
          if(jlmax.le.0 .or. jumax.le.0) then
            flag = 10*n + 1
            return
          end if
          call  nsf8
     *       (n,  r, ic,  ia, ja,
     *        isp(il), isp(jl), jlmax,  isp(iu), isp(jutmp), jumax,
     *        isp(q),  isp(im),  flag)
! ** error: error detected in nsf, nnf, nns, or nnt
          if(flag .ne. 0) return
!  ******  move ju next to jl  *****************************************
          jlmax = isp(il+n)-1
          ju    = jl + jlmax
          jumax = isp(iu+n)-1
          if(jumax .gt. 0) then
            do j = 1, jumax
              isp(ju+j-1) = isp(jutmp+j-1)
            end do
          end if
        end if
!
!  ******  call remaining subroutines  *********************************
        jlmax = isp(il+n)-1
        ju    = jl  + jlmax
        jumax = isp(iu+n)-1
        l     = (ju + jumax - 2 + lratio)  /  lratio    +    1
        lmax  = jlmax
        d     = l   + lmax
        u     = d   + n
        row   = nsp + 1 - n
        tmp   = row - n
        umax  = tmp - u
        esp = umax - jumax
!
        if((path-1) * (path-2) .eq. 0) then
! ** error:  insufficient storage
          if(umax .le. 0) then
            flag = 10*n + 1
            return
          end if
          call  nnf8
     *       (n,  r, c, ic,  ia, ja, a,  z,  b,
     *        isp(il), isp(jl), rsp(l), lmax,   rsp(d),
     *           isp(iu), isp(ju), rsp(u), umax,
     *        rsp(row),  rsp(tmp),  flag)
! ** error: error detected in nsf, nnf, nns, or nnt
          if(flag .ne. 0) return
          return
!
        end if
        if((path-3) .eq. 0) then
          call  nns8
     *       (n,  r, c,
     *        isp(il), isp(jl), rsp(l),  rsp(d),
     *          isp(iu), isp(ju), rsp(u),
     *        z,  b,  rsp(tmp))
!
        end if
        if((path-4) .eq. 0) then
          call  nnt8
     *       (n,  r, c,
     *        isp(il), isp(jl), rsp(l),  rsp(d),
     *          isp(iu), isp(ju), rsp(u),
     *        z,  b,  rsp(tmp))
        end if

        end
