!       ----------------------------------------------------------------
!
!               yale sparse matrix package - nonsymmetric codes
!                    solving the system of equations mx = b
!                        (uncompressed pointer storage)
!
!    i.   calling sequences
!         the coefficient matrix can be processed by an ordering routine
!    (e.g., to reduce fillin or ensure numerical stability) before using
!    the remaining subroutines.  if no reordering is done, then set
!    r(i) = c(i) = ic(i) = i  for i=1,...,n.  the calling sequence is --
!        (      (matrix ordering))
!         nsf   (symbolic factorization to determine where fillin will
!                 occur during numeric factorization)
!         nnf   (numeric factorization into product ldu of unit lower
!                 triangular matrix l, diagonal matrix d, and unit upper
!                 triangular matrix u, and solution of linear system)
!         nns   (solution of linear system for additional right-hand
!     or          side using ldu factorization from nnf)
!         nnt   (solution of transposed linear system for additional
!                 right-hand side using ldu factorization from nnf)
!
!    ii.  storage of sparse matrices
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
!         the strict triangular portions of the matrices l and u are
!    stored in the same fashion using the arrays  il, jl, l  and
!    iu, ju, u  respectively.  the diagonal entries of l and u are
!    assumed to be equal to one and are not stored.  the array d
!    contains the reciprocals of the diagonal entries of the matrix d.
!
!    iii. additional storage savings
!         in nsf, r and ic can be the same array in the calling
!    sequence if no reordering of the coefficient matrix has been done.
!         in nnf, r, c and ic can all be the same array if no reordering
!    has been done.  if only the rows have been reordered, then c and ic
!    can be the same array.  if the row and column orderings are the
!    same, then r and c can be the same array.  z and row can be the
!    same array.
!         in nns or nnt, r and c can be the same array if no reordering
!    has been done or if the row and column orderings are the same.  z
!    and b can be the same array;  however, then b will be destroyed.
!
!    iv.  parameters
!         following is a list of parameters to the programs.  names are
!    uniform among the various subroutines.  class abbreviations are --
!       n - integer variable
!       f - real variable
!       v - supplies a value to a subroutine
!       r - returns a result from a subroutine
!       i - used internally by a subroutine
!       a - array
!
! class   parameter
! ------+----------
! fva     a     - nonzero entries of the coefficient matrix m, stored
!                   by rows.
!                   size = number of nonzero entries in m.
! fva     b     - right-hand side b.
!                   size = n.
! nva     c     - ordering of the columns of m.
!                   size = n.
! fvra    d     - reciprocals of the diagonal entries of the matrix d.
!                   size = n.
! nr      flag  - error flag;  values and their meanings are --
!                    0     no errors detected
!                    n+k   null row in a  --  row = k
!                   2n+k   duplicate entry in a  --  row = k
!                   3n+k   insufficient storage for jl  --  row = k
!                   4n+1   insufficient storage for l
!                   5n+k   null pivot  --  row = k
!                   6n+k   insufficient storage for ju  --  row = k
!                   7n+1   insufficient storage for u
!                   8n+k   zero pivot  --  row = k
! nva     ia    - pointers to delimit the rows in a.
!                   size = n+1.
! nva     ic    - inverse of the ordering of the columns of m;  i.e.,
!                   ic(c(i) = i  for i=1,...n.
!                   size = n.
! nvra    il    - pointers to delimit the rows in l.
!                   size = n+1.
! nvra    iu    - pointers to delimit the rows in u.
!                   size = n+1.
! nva     ja    - column numbers corresponding to the elements of a.
!                   size = size of a.
! nvra    jl    - column numbers corresponding to the elements of l.
!                   size = jlmax.
! nv      jlmax - declared dimension of jl;  jlmax must be larger than
!                   the number of nonzero entries in the strict lower
!                   triangle of m plus fillin  (il(n+1)-1 after nsf).
! nvra    ju    - column numbers corresponding to the elements of u.
!                   size = jumax.
! nv      jumax - declared dimension of ju;  jumax must be larger than
!                   the number of nonzero entries in the strict upper
!                   triangle of m plus fillin  (iu(n+1)-1 after nsf).
! fvra    l     - nonzero entries in the strict lower triangular portion
!                   of the matrix l, stored by rows.
!                   size = lmax
! nv      lmax  - declared dimension of l;  lmax must be larger than
!                   the number of nonzero entries in the strict lower
!                   triangle of m plus fillin  (il(n+1)-1 after nsf).
! nv      n     - number of variables/equations.
! nva     r     - ordering of the rows of m.
!                   size = n.
! fvra    u     - nonzero entries in the strict upper triangular portion
!                   of the matrix u, stored by rows.
!                   size = umax.
! nv      umax  - declared dimension of u;  umax must be larger than
!                   the number of nonzero entries in the strict upper
!                   triangle of m plus fillin  (iu(n+1)-1 after nsf).
! fra     z     - solution x.
!                   size = n.
!
!
!       ----------------------------------------------------------------
!
!*** subroutine nsf
!*** symbolic ldu-factorization of a nonsymmetric sparse matrix
!      (uncompressed pointer storage)
!
        subroutine  nsf8
     *     (n, r,ic, ia,ja, il,jl,jlmax, iu,ju,jumax, q, im, flag)
!
!       input variables:   n, r,ic, ia,ja, jlmax, jumax.
!       output variables:  il,jl, iu,ju, flag.
!
!       parameters used internally:
! nia     q     - suppose m' is the result of reordering m;  if
!                   processing of the kth row of m' (hence the kth rows
!                   of l and u) is being done, then q(j) is initially
!                   nonzero if m'(k,j) is nonzero;  since values need
!                   not be stored, each entry points to the next
!                   nonzero;  for example, if  n=9  and the 5th row of
!                   m' is
!                           0 x x 0 x 0 0 x 0,
!                   then q will initially be
!                           a 3 5 a 8 a a 10 a 2        (a - arbitrary);
!                   q(n+1) points to the first nonzero in the row and
!                   the last nonzero points to  n+1;  as the algorithm
!                   proceeds, other elements of q are inserted in the
!                   list because of fillin.
!                   size = n+1.
! nia     im    - at each step in the factorization, im(i) is the last
!                   element in the ith row of u which needs to be
!                   considered in computing fillin.
!                   size = n.
!
!  internal variables--
!    jlptr - points to the last position used in  jl.
!    juptr - points to the last position used in  ju.
!
        integer  r(1), ic(1),  ia(1), ja(1),  il(1), jl(1),
     *     iu(1), ju(1),  q(1),  im(1),  flag,  qm, vj
!
!  ******  initialize pointers  ****************************************
        jlptr = 0
        il(1) = 1
        juptr = 0
        iu(1) = 1
!
!  ******  for each row of l and u  ************************************
        do k = 1, n
!  ******  set q to the reordered row of a  ****************************
          q(n+1) = n+1
          jmin = ia(r(k))
          jmax = ia(r(k)+1) - 1
! ** error:  null row in a
          if(jmin .gt. jmax) then
            flag = n + r(k)
            return
          end if
          do j = jmin, jmax
            vj = ic(ja(j))
            qm = n+1
            do
              m = qm
              qm = q(m)
              if(qm .ge. vj) exit
            end do
! ** error:  duplicate entry in a
            if(qm .eq. vj) then
              flag = 2*n + r(k)
              return
            end if
            q(m) = vj
            q(vj) = qm
          end do
!
!  ******  for each entry in the lower triangle  ***********************
          i = n+1
          do
            i = q(i)
            if(i .ge. k) exit
!  ******  l(k,i) will be nonzero, so add it to jl  ********************
              jlptr = jlptr+1
! ** error:  insufficient storage for jl
              if(jlptr .gt. jlmax) then
                flag = 3*n + k
                return
              end if
              jl(jlptr) = i
              qm = i
!  ******  inspect ith row for fillin, adjust im if possible  **********
              jmin = iu(i)
              jmax = im(i)
              do j = jmin, jmax
                vj = ju(j)
                if (vj.eq.k)  im(i) = j
                do
                  m = qm
                  qm = q(m)
                  if(qm .ge. vj) exit
                end do
                if(qm .ne. vj) then
                  q(m) = vj
                  q(vj) = qm
                  qm = vj
                end if
              end do
            end do
!
! ******  check for null pivot  ***************************************
! ** error:  null pivot
          if(i .ne. k) then
            flag = 5*n + k
            return
          end if
!  ******  remaining elements of q define structure of u(k, )  *********
          do
            i = q(i)
            if(i .gt. n) exit
            juptr = juptr+1
! ** error:  insufficient storage for ju
            if(juptr .gt. jumax) then
              flag = 6*n + k
              return
            end if
            ju(juptr) = i
          end do

!  ******  get ready for next row  *************************************
          im(k) = juptr
          il(k+1) = jlptr+1
          iu(k+1) = juptr+1
        end do

        flag = 0

        end
