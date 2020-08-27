subroutine Vof_Init_Random_Seed(problem_name)

  integer :: i, n, clock, int_p_name
  integer, dimension(:), allocatable :: seed
  character(len=20)  :: problem_name

  call random_seed(size = n)
  allocate(seed(n))
  call Convert_Problem_Name_To_Integer(problem_name,int_p_name)
  !call system_clock(count=clock)

  seed = int_p_name * (/ (i - 1, i = 1, n) /)
  call random_seed(put = seed)

  deallocate(seed)

end subroutine
