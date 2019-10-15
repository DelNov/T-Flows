!==============================================================================!
  logical function Math_Mod_Approx_String(s, t, tol)
!------------------------------------------------------------------------------!
!   For all i and j, d[i,j] will hold the Levenshtein distance between the     !
!   first "i" characters of "s" and the first "j" characters of "t".           !
!                                                                              !
!   Check for details here: https://tinyurl.com/y857p2hf                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(*)      :: s, t  ! source and target strings
  integer, optional :: tol   ! number of tolerated misstyped characters
!-----------------------------------[Locals]-----------------------------------!
  integer              :: i, j, m, n, cost, tolerance
  integer, allocatable :: d(:,:)
!==============================================================================!

  if( .not. present(tol) ) then
    tolerance = 3
  else
    tolerance = tol
  end if

  m = len_trim(s)
  n = len_trim(t)

  ! Declare int d[0..m, 0..n], and set each element in d to zero
  allocate(d(0:m, 0:n));  d(:,:) = 0

  ! Source prefixes can be transformed into ...
  ! ... empty string by dropping all characters
  do i = 1, m
    d(i, 0) = i
  end do

  ! Target prefixes can be reached from empty ...
  ! ... source prefix by inserting every character
  do j = 1, n
    d(0, j) = j
  end do

  do j = 1, n
    do i = 1, m
      if(s(i:i) .eq. t(j:j)) then
        cost = 0
      else
        cost = 1
      end if
      d(i, j) = min(d(i-1, j  ) + 1,           &   ! deletion
                    d(i,   j-1) + 1,           &   ! insertion
                    d(i-1, j-1) + cost)            ! substitution
    end do
  end do

  if( d(m,n) > 0 .and. d(m,n) <= tolerance ) then
    Math_Mod_Approx_String = .true.
  else
    Math_Mod_Approx_String = .false.
  end if

  end function
