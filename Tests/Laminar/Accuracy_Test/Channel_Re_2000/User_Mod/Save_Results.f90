!==============================================================================!
  subroutine User_Mod_Save_Results(flow, turb, mult, swarm, ts)
!------------------------------------------------------------------------------!
!   This subroutine reads name.1d file created by Convert or Generator and     !
!   averages the results in homogeneous directions.                            !
!                                                                              !
!   The results are then writen in files name_res.dat and name_res_plus.dat    !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),      target :: flow
  type(Turb_Type),       target :: turb
  type(Multiphase_Type), target :: mult
  type(Swarm_Type),      target :: swarm
  integer                       :: ts
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer      :: grid
  type(Bulk_Type), pointer      :: bulk
  type(Var_Type),  pointer      :: u
  real, pointer                 :: ub
  integer                       :: n, i, c, fu
  character(len=80)             :: coord_name, res_name
  real, allocatable             :: y_p(:), u_p(:), y_f(:)
  integer, allocatable          :: n_count(:), ind(:)
  logical                       :: exist
!==============================================================================!

  ! Take aliases
  grid   => flow % pnt_grid
  u      => flow % u
  ub     => flow % bulk % u

  ! Set the name for coordinate file
  call File_Mod_Set_Name(coord_name, extension='.1d')

  ! Set file name for results
  call File_Mod_Set_Name(res_name,          &
                         time_step=ts,      &
                         appendix='-res',   &
                         extension='.dat')

  !------------------!
  !   Read 1d file   !
  !------------------!
  inquire(file=coord_name, exist=exist)
  if(.not. exist) then
    if(this_proc < 2) then
      print *, '# Critial error: chan.1d file does not exist'
    end if
    stop
  end if

  call File_Mod_Open_File_For_Writing(coord_name, fu)

  ! Write the number of searching intervals
  call File_Mod_Read_Line(fu)
  read(line % tokens(1),*) n
  allocate(ind(n))
  allocate(y_f(n))

  ! Read the intervals positions
  do i = 1, n
    call File_Mod_Read_Line(fu)
    read(line % tokens(1),*) ind(i)
    read(line % tokens(2),*) y_f(i)
  end do
  close(fu)

  allocate(n_count(n)); n_count = 0
  allocate(y_p(n));     y_p     = 0.
  allocate(u_p(n));     u_p     = 0.

  !---------------------------!
  !   Summarize the results   !
  !---------------------------!
  do i = 1, n - 1
    do c = 1, grid % n_cells - grid % comm % n_buff_cells
      if(grid % yc(c) > y_f(i) .and. grid % yc(c) < y_f(i+1)) then

        n_count(i) = n_count(i) + 1

        y_p(i) = y_p(i) + grid % yc(c)
        u_p(i) = u_p(i) + u % n(c)

      end if
    end do
  end do

  ! Average over all processors
  do i = 1, n-1
    call Comm_Mod_Global_Sum_Int (n_count(i))
    call Comm_Mod_Global_Sum_Real(y_p(i))
    call Comm_Mod_Global_Sum_Real(u_p(i))
  end do

  call Comm_Mod_Wait

  do i = 1, n-1
    if(n_count(i) .ne. 0) then
      y_p(i) = y_p(i) / n_count(i)
      u_p(i) = u_p(i) / n_count(i)
    end if
  end do

  call File_Mod_Open_File_For_Writing(res_name, fu)

  write(fu,'(a,1es17.7e3)') '# U_bulk:', ub

  do i = 1, n
    if(n_count(i) .ne. 0) then
      write(fu,'(2es17.7e3)')  y_p(i),   &  !  1
                               u_p(i)/ub    !  2
    end if
  end do

  close(fu)

  deallocate(n_count)
  deallocate(y_p)
  deallocate(u_p)
  deallocate(ind)
  deallocate(y_f)

  if(this_proc < 2)  write(*,*) '# Finished with User_Mod_Save_Results.f90.'

  end subroutine
