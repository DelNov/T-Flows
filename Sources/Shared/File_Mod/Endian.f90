!   https://atmos.washington.edu/~salathe/osx_unix/endian.html   !

!==============================================================================!
  logical function Is_Machine_Little_Endian()
!------------------------------------------------------------------------------!
!>  Checks if the machine on which it is running uses little-endian byte order.
!>  If it is little-endian, it returns true, otherwise it returns false.
!------------------------------------------------------------------------------!
!   Functionality:                                                             !
!                                                                              !
!   * The function sets i to 1. In a little-endian system, the least           !
!     significant byte (which is 1) will be stored in the first byte of i.     !
!     The equivalence statement allows j(1) to access this first byte.         !
!   * If j(1) equals 1, it concludes that the machine is little-endian;        !
!     otherwise, it's not.                                                     !
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer*1   :: j(2)
  integer*2   :: i
  equivalence :: (i,j)
!==============================================================================!

  i = 1
  if(j(1).eq.1) then
    Is_Machine_Little_Endian = .true.
  else
    Is_Machine_Little_Endian = .false.
  end if

  end function

!==============================================================================!
  subroutine Byte_Swap_Int_2(k)
!------------------------------------------------------------------------------!
!>  The subroutine swaps the bytes of a 2-byte integer.
!>  (I don't think this function is ever used which indicates that we are
!>  working in environments where endianness consistency isn't a concern.)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  integer(2) :: k
!-----------------------------------[Locals]-----------------------------------!
  integer*1   :: ii(2), jj(2)
  integer*2   :: i, j
  equivalence :: (i,ii)
  equivalence :: (j,jj)
!==============================================================================!

  i = k

  jj(1) = ii(2)
  jj(2) = ii(1)

  k = j

  end subroutine


!==============================================================================!
  subroutine Byte_Swap_Real_4(k)
!------------------------------------------------------------------------------!
!>  The subroutine swaps the bytes of a 4-byte real number.
!>  (I don't think this function is ever used which indicates that we are
!>  working in environments where endianness consistency isn't a concern.)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  real(4) :: r
!-----------------------------------[Locals]-----------------------------------!
  integer(1)  :: ii(4), jj(4)
  real(4)     :: s, t
  equivalence :: (s,ii)
  equivalence :: (t,jj)
!==============================================================================!

  s = r

  jj(1) = ii(4)
  jj(2) = ii(3)
  jj(3) = ii(2)
  jj(4) = ii(1)

  r = t

  end subroutine
