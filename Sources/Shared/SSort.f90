      subroutine SSORT (X, IY, N, KFLAG)
      implicit none
!
!    Example of a Selection Sort   Using a Fortran 90 Intrinsic Function
!
!***BEGIN PROLOGUE  SSORT
!***PURPOSE  Sort an array and make the same interchanges in
!            an auxiliary array.  The array is sorted in
!            decreasing order.
!***TYPE      SINGLE PRECISION
!***KEYWORDS  SORT, SORTING
!
!   Description of Parameters
!      X - array of values to be sorted   (usually abscissas)
!      IY - array to be carried with X (all swaps of X elements are
!          matched in IY .  After the sort IY(J) contains the original
!          postition of the value X(J) in the unsorted X array.
!      N - number of values in array X to be sorted
!      KFLAG - Not used in this implementation
!
!***REVISION HISTORY  (YYMMDD!
!   950310  DATE WRITTEN
!   John Mahaffy
!***end PROLOGUE  SSORT
!     .. Scalar Arguments ..
      integer KFLAG, N
!     .. Array Arguments ..  -----NOTE the 2 new ways of declaring array size
      real X(1:N)
      integer IY(N)
!     .. Local Scalars ..
      real TEMP
      integer I, ISWAP(1), ITEMP, ISWAP1
!     .. External Subroutines ..
!     None
!     .. Intrinsic Functions ..
      INTRINSIC MAXLOC
!
!
!    MAXLOC is a FORTRAN 90 function that returns the index value for the
!    maximum element in the array
!***FIRST EXECUTABLE STATEMENT  SSORT
!
      DO 200 I=1,N-1
!<a name=1><font color=FF0000>
         ISWAP=MAXLOC(X(I:N))

!</font>
         ISWAP1=ISWAP(1)+I-1
         IF(ISWAP1.NE.I) THEN
           TEMP=X(I)
            X(I)=X(ISWAP1)
            X(ISWAP1)=TEMP
            ITEMP=IY(I)
            IY(I)=IY(ISWAP1)
            IY(ISWAP1)=ITEMP
         ENDIF
  200 CONTINUE
      RETURN
      end
