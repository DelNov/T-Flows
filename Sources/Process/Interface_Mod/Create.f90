!==============================================================================!
  subroutine Interface_Mod_Create(inter, grid, n_dom)
!------------------------------------------------------------------------------!
!   Create interface between two grids.                                        !
!                                                                              !
!   Each processors saves the whole faces_1 and faces_2                        !
!   Cells surrounding these interfaces are stored by their global values (at   !
!   first) and interfaces are sorted according to them.  This sorting in done  !
!   in such a way to preserve mappings to local cells                          !
!   The values are later filled from each processor to the global face-based   !
!   array                                                                      !
!   Then they are summed up over all processors                                !
!   Then distributed back to local cells                                       !
!                                                                              !
!                                                                              !
!                                                                              !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Interface_Type)    :: inter(MD, MD)
  type(Grid_Type), target :: grid(MD)
  integer                 :: n_dom

  integer           :: c1, c2, d1, d2, k, nks, s, n1, n2, n1_tot, n2_tot, n_tot
  integer           :: off_1, off_2
  logical           :: found
  character(len=80) :: keys(128)
!==============================================================================!

  !--------------------------------------------!
  !                                            !
  !   Browse through all domain combinations   !
  !                                            !
  !--------------------------------------------!
  do d1 = 1, n_dom
    do d2 = 1, n_dom
      inter(d1, d2) % n_faces = 0
      if(d1 .ne. d2) then

        !-------------------------------------------------------------------!
        !   Try to find specified interface condition between d1 and d2   !
        !-------------------------------------------------------------------!
        call Control_Mod_Position_At_Three_Keys('INTERFACE_CONDITION',    &
                                                trim(problem_name(d1)),  &
                                                trim(problem_name(d2)),  &
                                                found,                    &
                                                verbose=.false.)

        !-----------------------------------------------------!
        !   Found specification between domains d1 and d2   !
        !-----------------------------------------------------!
        if(found) then
          call Control_Mod_Read_Strings_On('BOUNDARY_CONDITIONS',  &
                                           keys, nks, .false.)
          do k = 1, nks
            call To_Upper_Case(keys(k))
          end do

          !---------------------------------------!
          !   Count boundary cells on each side   !
          !---------------------------------------!
          n1 = 0
          do s = 1, grid(d1) % n_faces
            c2 = grid(d1) % faces_c(2,s)
            if(c2 < 0) then
              if(Grid_Mod_Bnd_Cond_Name(grid(d1), c2) .eq. keys(1)) then
                n1 = n1 + 1
              end if
            end if
          end do

          n2 = 0
          do s = 1, grid(d2) % n_faces
            c2 = grid(d2) % faces_c(2,s)
            if(c2 < 0) then
              if(Grid_Mod_Bnd_Cond_Name(grid(d2), c2) .eq. keys(2)) then
                n2 = n2 + 1
              end if
            end if
          end do

          n1_tot = n1
          n2_tot = n2
          call Comm_Mod_Global_Sum_Int(n1_tot)
          call Comm_Mod_Global_Sum_Int(n2_tot)

          if(n1_tot .ne. n2_tot) then
            print *, '# Number of cells at the interface between ',  &
                      trim(problem_name(d1)), ' and ',               &
                      trim(problem_name(d2)), ' is not the same!'
            print *, '# Only conformal mappings are supported.  Exiting!'
            stop
          else
            n_tot = n1_tot
            inter(d1, d2) % n_faces = n1_tot
            print '(a,a,a,a,a,i6,a)', ' # Domains ', trim(problem_name(d1)),  &
                                            ' and ', trim(problem_name(d2)),  &
                      ' are connected with ', n_tot, ' interface cells!'
          end if

          !---------------------------------------!
          !   Store boundary cells on each side   !
          !---------------------------------------!
          inter(d1, d2) % n_faces = n1  ! the same as n2
          allocate(inter(d1, d2) % glo_1(n_tot)); inter(d1, d2) % glo_1 = 0
          allocate(inter(d1, d2) % glo_2(n_tot)); inter(d1, d2) % glo_2 = 0
          allocate(inter(d1, d2) % bnd_1(n_tot)); inter(d1, d2) % bnd_1 = 0
          allocate(inter(d1, d2) % bnd_2(n_tot)); inter(d1, d2) % bnd_2 = 0
          allocate(inter(d1, d2) % phi_1(n_tot)); inter(d1, d2) % phi_1 = 0.0
          allocate(inter(d1, d2) % phi_2(n_tot)); inter(d1, d2) % phi_2 = 0.0

          off_1 = 0
          off_2 = 0

          ! On the side of d1
          n1 = 0
          do s = 1, grid(d1) % n_faces
            c1 = grid(d1) % faces_c(1,s)
            c2 = grid(d1) % faces_c(2,s)
            if(c2 < 0) then
              if(Grid_Mod_Bnd_Cond_Name(grid(d1), c2) .eq. keys(1)) then
                n1 = n1 + 1
                inter(d1, d2) % glo_1(n1 + off_1) = grid(d1) % comm % cell_glo(c1)
                inter(d1, d2) % bnd_1(n1 + off_1) = grid(d1) % comm % cell_glo(c2)
                xf_1 (n1 + off_1) = grid(d1) % xf(s)
                yf_1 (n1 + off_1) = grid(d1) % yf(s)
                zf_1 (n1 + off_1) = grid(d1) % zf(s)
              end if
            end if
          end do

          ! On the side of d2
          n2 = 0
          do s = 1, grid(d2) % n_faces
            c1 = grid(d2) % faces_c(1,s)
            c2 = grid(d2) % faces_c(2,s)
            if(c2 < 0) then
              if(Grid_Mod_Bnd_Cond_Name(grid(d2), c2) .eq. keys(2)) then
                n2 = n2 + 1
                inter(d1, d2) % glo_2(n2 + off_2) = grid(d2) % comm % cell_glo(c1)
                inter(d1, d2) % bnd_2(n2 + off_2) = grid(d2) % comm % cell_glo(c2)
                xf_2 (n2 + off_2) = grid(d2) % xf(s)
                yf_2 (n2 + off_2) = grid(d2) % yf(s)
                zf_2 (n2 + off_2) = grid(d2) % zf(s)
              end if
            end if
          end do

          ! Sort interfaces from domain 1
          call Sort_Mod_3_Real_Carry_2_Int(xf_1(1:n_tot),  &
                                           yf_1(1:n_tot),  &
                                           zf_1(1:n_tot),  &
                                           inter(d1, d2) % glo_1(1:n_tot), &
                                           inter(d1, d2) % bnd_1(1:n_tot))
          ! Sort interfaces from domain 2
          call Sort_Mod_3_Real_Carry_2_Int(xf_2(1:n_tot),  &
                                           yf_2(1:n_tot),  &
                                           zf_2(1:n_tot),  &
                                           inter(d1, d2) % glo_2(1:n_tot), &
                                           inter(d1, d2) % bnd_2(1:n_tot))
          ! Write some debugging information
          do s = 1, inter(d1, d2) % n_faces
            write(100+this_proc, '(i4, 2i8, 3f10.5, 2i8, 3f10.5)')   &
                s,                         &
                inter(d1, d2) % glo_1(s),  &
                inter(d1, d2) % bnd_1(s),  &
                xf_1(s),                   &
                yf_1(s),                   &
                zf_1(s),                   &
                inter(d1, d2) % glo_2(s),  &
                inter(d1, d2) % bnd_2(s),  &
                xf_2(s),                   &
                yf_2(s),                   &
                zf_2(s)
          end do
        end if  ! if this interface is found in control mod
      end if

    end do
  end do

  end subroutine
