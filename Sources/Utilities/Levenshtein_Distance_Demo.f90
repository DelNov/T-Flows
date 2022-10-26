!==============================================================================!
  subroutine Levenshtein_Distance(s, t)
!------------------------------------------------------------------------------!
!   For all i and j, d[i,j] will hold the Levenshtein distance between the     !
!   first "i" characters of "s" and the first "j" characters of "t".           !
!                                                                              !
!   Check for details here: https://tinyurl.com/y857p2hf                       !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  character(*) :: s, t  ! source and target strings
!-----------------------------------[Locals]-----------------------------------!
  integer              :: i, j, m, n, cost
  integer, allocatable :: d(:,:)
!==============================================================================!

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

  ! do j = 0, n
  !   write(*,'(999i4)'), d(0:m, j)
  ! end do

  print *, '# Distance between ', trim(s) , ' and ', trim(t), ' is ', d(m,n)

  end subroutine

!==============================================================================!
  program Levenshtein_Distance_Demo
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  character(len=80) :: from_file, from_database
!==============================================================================!

  from_file     = "KITTEN"
  from_database = "SITTING"
  call Levenshtein_Distance(from_file, from_database)

  from_file     = "JUMBO"
  from_database = "DUMBO"
  call Levenshtein_Distance(from_file, from_database)

  from_file     = "NUMBER_OF_TIME_STEPS"
  from_database = "NUBER_OF_TIME_STEP"
  call Levenshtein_Distance(from_file, from_database)

  from_file     = "PROBLEM_NAME"
  from_database = "PRUBLEM_NAME"
  call Levenshtein_Distance(from_file, from_database)

  end program
