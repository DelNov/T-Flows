!==============================================================================!
  subroutine Comm_Mod_Create_Buffers(grid)
!------------------------------------------------------------------------------!
!   Reads: name.buf                                                            !
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Grid_Mod
  use Tokenizer_Mod
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: sub, subo
  character(len=80) :: name_in
  integer           :: c1, c2, s, buf_cnt
!==============================================================================!
!   Each subdomain needs two buffers: a send buffer and a receive buffer.      !
!   A receive buffer will be stored as aditional cells for each subdomain.     !
!   So each subdomain will have grid % n_cells cells, which entials physical   !
!   and buffer cells. It is handy to do it that way, because most of the       !
!   algorythms can remain the same as in sequential run. On the other hand,    !
!   a sending buffer has to be allocated in a new separate array called        !
!   simply buffer(). An additional array is needed to keep track of all the    !
!   indexes. That one is called buffind().                                     !
!------------------------------------------------------------------------------!

  if(n_proc .eq. 0) return

  allocate (grid % comm % nbb_s(0:n_proc))  ! why from zero?
  allocate (grid % comm % nbb_e(0:n_proc))

  ! Initialize
  do sub = 0, n_proc
    grid % comm % nbb_s(sub) = grid % n_cells - grid % comm % n_buff_cells
    grid % comm % nbb_e(sub) = grid % n_cells - grid % comm % n_buff_cells
  end do

  ! Browse through other sub-domains
  do subo = 1, n_proc

    ! Connection with new domain will start ...
    ! ... where the domain with previous ended
    grid % comm % nbb_s(subo) = grid % comm % nbb_e(subo-1) + 1
    grid % comm % nbb_e(subo) = grid % comm % nbb_e(subo-1)

    ! Initialize buffer connections with this subdomain
    buf_cnt = 0

    if(subo .ne. this_proc) then

      do s = 1, grid % n_faces
        c1 = grid % faces_c(1,s)
        c2 = grid % faces_c(2,s)
        if(c2 > 0) then

          ! If c2 is a buffer face
          if( (grid % comm % proces(c1) .eq. this_proc) .and.  &
              (grid % comm % proces(c2) .eq. subo) ) then
            buf_cnt = buf_cnt + 1                ! increase buffer cell count
            grid % comm % buffer_index(c2) = c1  ! buffer send index
          end if
          if( (grid % comm % proces(c2) .eq. this_proc) .and.  &
              (grid % comm % proces(c1) .eq. subo) ) then
            print *, '# FATAL ERROR! Shouldn''t be here at all!'
            print *, '# Exiting now!'
            stop
          end if
        end if  ! c2 > 0
      end do    ! through faces

    end if

    ! Set the end of the current buffer
    grid % comm % nbb_e(subo) = grid % comm % nbb_e(subo) + buf_cnt

  end do

  end subroutine
