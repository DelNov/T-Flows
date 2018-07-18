!==============================================================================!
  subroutine Save_Ini(grid)
!------------------------------------------------------------------------------!
!   
!------------------------------------------------------------------------------!
!----------------------------------[Modules]-----------------------------------!
  use Name_Mod, only: problem_name
  use Const_Mod
  use Flow_Mod
  use Les_Mod
  use Comm_Mod, only: this_proc
  use Rans_Mod
  use Tokenizer_Mod
  use Grid_Mod
  use Control_Mod
!------------------------------------------------------------------------------!
  implicit none
  type(Grid_Type) :: grid
!-----------------------------------[Locals]-----------------------------------!
  integer           :: c, nn
  character(len=80) :: name_out, answer
  character(len=5)  :: ext
!==============================================================================!

  call Control_Mod_Save_Initial_Solution_Name(name_out)

  answer=name_out
  call To_Upper_Case(answer)
  if(answer .eq. 'SKIP') return

  ! save the name
  answer = problem_name
  problem_name = name_out
  nn = 0
  ext = '.xyz'
  call Name_File(this_proc, name_out, ext)
  open(9,file=name_out)
  do c = 1, grid % n_cells
    nn = nn + 1
  end do    ! through centers 
  write(9,'(I10)') nn
  do c = 1, grid % n_cells
    write(9,'(3E25.8)') grid % xc(c),grid % yc(c),grid % zc(c)
  end do    ! through centers 
  close(9)

  ext(1:4) = '.  '
  ext(2:4) = u % name
  call Name_File(this_proc, name_out, ext)
  open(9,file=name_out)
  do c = 1, grid % n_cells
    write(9,'(7E18.8)') U % n(c), U % o(c), U % a(c), U % a_o(c),  &
                        U % d_o(c), U % c(c), U % c_o(c)
  end do    ! through centers 
  close(9)

  ext(2:4) = v % name
  call Name_File(this_proc, name_out, ext)
  open(9,file=name_out)
  do c = 1, grid % n_cells
    write(9,'(7E18.8)') V % n(c), V % o(c), V % a(c), V % a_o(c),  &
                        V % d_o(c), V % c(c), V % c_o(c)
  end do    ! through centers 
  close(9)

  ext(2:4) = w % name
  call Name_File(this_proc, name_out, ext)
  open(9,file=name_out)
  do c = 1, grid % n_cells
    write(9,'(7E18.8)') W % n(c), W % o(c), W % a(c), W % a_o(c),  &
                        W % d_o(c), W % c(c), W % c_o(c)
  end do    ! through centers 
  close(9)

  ext(2:4) = p % name
  call Name_File(this_proc, name_out, ext)
  open(9,file=name_out)
  do c = 1, grid % n_cells
    write(9,'(5E18.8)') P % n(c), PP % n(c), p % x(c), p % y(c), p % z(c)
  end do    ! through centers 
  close(9)
 
  if(heat_transfer .eq. YES) then 
    ext(2:4) = t % name
    call Name_File(this_proc, name_out, ext)
    open(9,file=name_out)
    do c = 1, grid % n_cells
      write(9,'(7E18.8)') T % n(c), T % o(c), T % a(c), T % a_o(c),  &
                          T % d_o(c), T % c(c), T % c_o(c)
    end do    ! through centers 
    close(9)
  end if 
 
  if(turbulence_model .eq. K_EPS_ZETA_F) then
    ext(2:4) = kin % name
    call Name_File(this_proc, name_out, ext)
    open(9,file=name_out)
    do c = 1, grid % n_cells
      write(9,'(7E18.8)') kin % n(c), kin % o(c), kin % a(c), kin % a_o(c),  &
                          kin % d_o(c), kin % c(c), kin % c_o(c)
    end do    ! through centers 
    close(9)

    ext(2:4) = eps % name
    call Name_File(this_proc, name_out, ext)
    open(9,file=name_out)
    do c = 1, grid % n_cells
      write(9,'(7E18.8)') eps % n(c), eps % o(c), eps % a(c), eps % a_o(c),  &
                          eps % d_o(c), eps % c(c), eps % c_o(c)
    end do    ! through centers 
    close(9)

    ext(2:5) = zeta % name
    call Name_File(this_proc, name_out, ext)
    open(9,file=name_out)
    do c = 1, grid % n_cells
      write(9,'(7E18.8)') zeta % n(c), zeta % o(c), zeta % a(c), zeta % a_o(c),  &
                          zeta % d_o(c), zeta % c(c), zeta % c_o(c)
    end do    ! through centers 
    close(9)

    ext(2:4) = f22 % name
    call Name_File(this_proc, name_out, ext)
    open(9,file=name_out)
    do c = 1, grid % n_cells
      write(9,'(7E18.8)') f22 % n(c), f22 % o(c),  &
                          f22 % d_o(c), f22 % c(c), f22 % c_o(c)
    end do    ! through centers 
    close(9)
  end if

  if(turbulence_model .eq. K_EPS) then
    ext(2:4) = kin % name
    call Name_File(this_proc, name_out, ext)
    open(9,file=name_out)
    do c = 1, grid % n_cells
      write(9,'(7E18.8)') kin % n(c), kin % o(c), kin % a(c), kin % a_o(c),  &
                          kin % d_o(c), kin % c(c), kin % c_o(c)
    end do    ! through centers 
    close(9)

    ext(2:4) = eps % name
    call Name_File(this_proc, name_out, ext)
    open(9,file=name_out)
    do c = 1, grid % n_cells
      write(9,'(7E18.8)') eps % n(c), eps % o(c), eps % a(c), eps % a_o(c),  &
                          eps % d_o(c), eps % c(c), eps % c_o(c)
    end do    ! through centers 
    close(9)
  end if

  problem_name = answer 

  end subroutine
