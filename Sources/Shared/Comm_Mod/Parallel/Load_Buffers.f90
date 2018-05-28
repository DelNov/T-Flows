!==============================================================================!
  subroutine Comm_Mod_Load_Buffers
!------------------------------------------------------------------------------!
!   Reads: name.buf                                                            !
!----------------------------------[Modules]-----------------------------------!
  use Tokenizer_Mod
!------------------------------------------------------------------------------!
  implicit none
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, dummy 
  integer           :: sub, subo, n_bnd_cells_sub
  character(len=80) :: name_in
!==============================================================================!
!   Each subdomain needs two buffers: a send buffer and a receive buffer.
!   A receive buffer will be stored as aditional boundary cells for each
!   subdomain. So each subdomain will have NBC physical boundary faces
!   and nbb_C-NBC buffer bounndary cells. It is handy to do it that way,
!   because most of the algorythms can remain the same as they are now.
!   They won't even "know" that they use values from other processors.
!   On the other hand, a sending buffer has to be allocated in a new 
!   separate array called simply buffer(). An additional array is needed 
!   to keep track of all the indexes. That one is called buffer_index().
!   buffer_index() has stored cell numbers from it's own subdomain so that
!   later they can be copied with (well, something like that):
!   do i=1,BUFFSIZ
!     buffer(i) = U(buffer_index(i))
!   end do
!------------------------------------------------------------------------------!

  if(n_proc .eq. 0) return

  call Name_File(this_proc, name_in, '.buf')
  open(9, file=name_in)
  if(this_proc < 2) print *, '# Now reading the file:', name_in

  allocate (nbb_s(0:n_proc))
  allocate (nbb_e(0:n_proc))

  ! Number of physical boundary cells
  call Tokenizer_Mod_Read_Line(9)
  read(line % whole,*) n_bnd_cells_sub

  ! Initialize 
  do sub=0,n_proc
    nbb_s(sub) = -(n_bnd_cells_sub) 
    nbb_e(sub) = -(n_bnd_cells_sub)
  end do

  ! Fill the indexes and the buffers
  do sub=1,n_proc
    if(sub  /=  this_proc) then

      ! Connections with subdomain          
      call Tokenizer_Mod_Read_Line(9)
      read(line % whole,*) subo 

      ! Number of local connections with subdomain sub 
      call Tokenizer_Mod_Read_Line(9)
      read(line % whole,*) nbb_e(sub)

      nbb_s(sub) = nbb_e(sub-1) - 1  
      nbb_e(sub) = nbb_s(sub) - nbb_e(sub) + 1

      do c=nbb_s(sub),nbb_e(sub),-1
        call Tokenizer_Mod_Read_Line(9)
        read(line % whole,*) dummy, buffer_index(c) 
      end do 
    else
      nbb_s(sub) = nbb_e(sub-1)-1  ! just to become "sloppy" 
      nbb_e(sub) = nbb_e(sub-1)    ! this_proc will be needed for next 
    end if
  end do   ! through subdomains

  close(9)

  ! Correct the "sloppy" indexes
  do sub=1,n_proc
    if(nbb_e(sub)  > nbb_s(sub)) then  
      nbb_s(sub) = -1 
      nbb_e(sub) = 0 
    end if
  end do 

  call Comm_Mod_Wait

  end subroutine
