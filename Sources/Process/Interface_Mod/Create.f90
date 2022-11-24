!==============================================================================!
  subroutine Interface_Mod_Create(inter, Grid, n_dom)
!------------------------------------------------------------------------------!
!   Create interface between two grids.                                        !
!                                                                              !
!   The subroutine works in the following way:                                 !
!   1. It allocates the entire interface on all processors.                    !
!   2. For both sides of the interface (1 and 2) it stores face coordinates    !
!      (xf_1 to zf_2), cell inside (ic_1, ic_2), bundary cells (ib_1, ib_2).   !
!      and processor ids (ip_1, ip_2).                                         !
!   3. While doing the step above, it takes care of the offsets of each        !
!      processor (off_1 and off_2).                                            !
!   4. Then takes global summs of all the arrays mentioned in step 2.          !
!   5. Sorts the above arrays by their x, y and z coordinates, carrying        !
!      information on cells inside (ic_1, ic), boundary cells (ib_1, ib_2)     !
!      and processor ids (ip_1, ip_2)                                          !
!   6. The sorting from step 5 gives mapping of local cells inside, and local  !
!      boundary cells to the global interface.                                 !
!   7. The mapping information is finally stored in arrays cell_1, cell_2      !
!      (for inside cells), bcel_1, bcel_2 for boundary cells and face_1        !
!      and face_2 for global faces at the interface.                           !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Interface_Type)    :: inter(MD, MD)
  type(Grid_Type), target :: Grid(MD)
  integer                 :: n_dom
!-----------------------------------[Locals]-----------------------------------!
  integer,         allocatable :: off_1(:)
  integer,         allocatable :: off_2(:)
  integer                      :: n, n1, n2, n1_tot, n2_tot, n_tot
  integer                      :: pos, c, c1, c2, d1, d2, k, nks, s, p
  logical                      :: found
  character(SL)                :: keys(128)
  real,    contiguous, pointer :: xf_1(:), yf_1(:), zf_1(:)
  real,    contiguous, pointer :: xf_2(:), yf_2(:), zf_2(:)
  integer, contiguous, pointer :: ic_1(:), ib_1(:), ip_1(:)
  integer, contiguous, pointer :: ic_2(:), ib_2(:), ip_2(:)
!==============================================================================!

  call Work % Connect_Real_Face(xf_1, yf_1, zf_1, xf_2, yf_2, zf_2)
  call Work % Connect_Int_Face(ic_1, ib_1, ip_1, ic_2, ib_2, ip_2)

  ! Allocate memory for offsets
  if(n_proc > 1) then
    allocate(off_1(n_proc)); off_1 = 0
    allocate(off_2(n_proc)); off_2 = 0
  end if

  !--------------------------------------------!
  !                                            !
  !   Browse through all domain combinations   !
  !                                            !
  !--------------------------------------------!
  do d1 = 1, n_dom
    do d2 = 1, n_dom

      ! Innitialize number of faces at the interface
      inter(d1, d2) % n_tot = 0

      ! Store pointers to the grids surrounding the interface
      inter(d1, d2) % pnt_grid1 => Grid(d1)
      inter(d1, d2) % pnt_grid2 => Grid(d2)

      if(d1 .ne. d2) then

        xf_1(:) = 0;  yf_1(:) = 0; zf_1(:) = 0
        xf_2(:) = 0;  yf_2(:) = 0; zf_2(:) = 0
        ic_1(:) = 0;  ib_1(:) = 0; ip_1(:) = 0
        ic_2(:) = 0;  ib_2(:) = 0; ip_2(:) = 0

        !-----------------------------------------------------------------!
        !   Try to find specified interface condition between d1 and d2   !
        !-----------------------------------------------------------------!
        call Control_Mod_Position_At_Three_Keys('INTERFACE_CONDITION',    &
                                                trim(problem_name(d1)),   &
                                                trim(problem_name(d2)),   &
                                                found,                    &
                                                verbose=.false.)

        !---------------------------------------------------!
        !   Found specification between domains d1 and d2   !
        !---------------------------------------------------!
        if(found) then
          call Control_Mod_Read_Strings_On('BOUNDARY_CONDITIONS',  &
                                           keys, nks, .false.)
          do k = 1, nks
            call String % To_Upper_Case(keys(k))
          end do

          !---------------------------------------!
          !   Count boundary cells on each side   !
          !---------------------------------------!

          ! Innitialize counter in the first domain
          n1 = 0

          ! On physical boundary cells
          do s = 1, Grid(d1) % n_faces
            c2 = Grid(d1) % faces_c(2,s)
            if(c2 < 0) then
              if(Grid(d1) % Bnd_Cond_Name(c2) .eq. keys(1)) then
                n1 = n1 + 1
              end if
            end if
          end do

          ! On periodic faces of domain 1
          do s = 1, Grid(d1) % n_faces
            if(Grid(d1) % faces_s(s) > 0) then  ! only if it has a shadow
              c1 = Grid(d1) % faces_c(1,s)
              if(Grid(d1) % Bnd_Cond_Name(s) .eq. keys(1) .and.  &
                 Grid(d1) % comm % cell_proc(c1) .eq. this_proc) then
                n1 = n1 + 1
              end if
            end if
          end do
          inter(d1, d2) % n1_sub = n1

          ! Innitialize counter in the second domain
          n2 = 0

          ! On physical boundary cells
          do s = 1, Grid(d2) % n_faces
            c2 = Grid(d2) % faces_c(2,s)
            if(c2 < 0) then
              if(Grid(d2) % Bnd_Cond_Name(c2) .eq. keys(2)) then
                n2 = n2 + 1
              end if
            end if
          end do

          ! On periodic faces
          do s = 1, Grid(d2) % n_faces
            if(Grid(d2) % faces_s(s) > 0) then  ! only if it has a shadow
              c1 = Grid(d2) % faces_c(1,s)
              if(Grid(d2) % Bnd_Cond_Name(s) .eq. keys(2) .and.  &
                 Grid(d2) % comm % cell_proc(c1) .eq. this_proc) then
                n2 = n2 + 1
              end if
            end if
          end do
          inter(d1, d2) % n2_sub = n2

          !------------------------------------------!
          !     Check if the mapping is conformal    !
          !   (Simply by comparing number of faces   !
          !      on each side of the interface)      !
          !------------------------------------------!
          n1_tot = n1
          n2_tot = n2
          call Comm_Mod_Global_Sum_Int(n1_tot)
          call Comm_Mod_Global_Sum_Int(n2_tot)

          if(n1_tot .ne. n2_tot) then
            if(this_proc < 2) then
              print *, '# Number of cells at the interface between ',  &
                        trim(problem_name(d1)), ' and ',               &
                        trim(problem_name(d2)), ' is not the same!'
              print *, '# Only conformal mappings are supported.  Exiting!'
            end if
            call Comm_Mod_End
            stop
          else
            n_tot = n1_tot
            inter(d1, d2) % n_tot = n_tot
            if(this_proc < 2) then
              print '(5a,i6,a)', ' # Domains ', trim(problem_name(d1)),  &
                                       ' and ', trim(problem_name(d2)),  &
                        ' are connected with ', n_tot, ' interface cells!'
            end if
          end if

          !-----------------------------------------------------!
          !          Store boundary cells on each side          !
          !   (Inside cells in case periodicity is specified)   !
          !-----------------------------------------------------!
          ic_1(1:n_tot) = 0; ic_2(1:n_tot) = 0
          ib_1(1:n_tot) = 0; ib_2(1:n_tot) = 0
          ip_1(1:n_tot) = 0; ip_2(1:n_tot) = 0
          allocate(inter(d1, d2) % phi_1(n_tot, MAX_VARS_INTERFACE))
          allocate(inter(d1, d2) % phi_2(n_tot, MAX_VARS_INTERFACE))
          inter(d1, d2) % phi_1(:,:) = 0.0
          inter(d1, d2) % phi_2(:,:) = 0.0

          ! Work out offsets for each interface for each processor
          if(n_proc > 1) then
            off_1(1:n_proc)  = 0
            off_2(1:n_proc)  = 0
            off_1(this_proc) = n1
            off_2(this_proc) = n2
            call Comm_Mod_Global_Sum_Int_Array(n_proc, off_1)
            call Comm_Mod_Global_Sum_Int_Array(n_proc, off_2)

            do p = n_proc, 2, -1
              off_1(p) = sum(off_1(1:p-1))
              off_2(p) = sum(off_2(1:p-1))
            end do
            off_1(1) = 0
            off_2(1) = 0
          end if

          ! On the side of d1
          n1 = 0
          do s = 1, Grid(d1) % n_faces
            c1 = Grid(d1) % faces_c(1,s)
            c2 = Grid(d1) % faces_c(2,s)
            if(c2 < 0) then
              if(Grid(d1) % Bnd_Cond_Name(c2) .eq. keys(1)) then
                n1 = n1 + 1
                pos = n1
                if(n_proc > 1) pos = pos + off_1(this_proc)
                ic_1(pos) = c1
                ib_1(pos) = c2
                ip_1(pos) = this_proc
                xf_1(pos) = Grid(d1) % xf(s)
                yf_1(pos) = Grid(d1) % yf(s)
                zf_1(pos) = Grid(d1) % zf(s)
              end if
            end if
          end do
          do s = 1, Grid(d1) % n_faces
            c1 = Grid(d1) % faces_c(1,s)
            c2 = Grid(d1) % faces_c(2,s)
            if(Grid(d1) % faces_s(s) > 0) then  ! only if it has a shadow
              if(Grid(d1) % Bnd_Cond_Name(s) .eq. keys(1) .and.  &
                 Grid(d1) % comm % cell_proc(c1) .eq. this_proc) then
                n1 = n1 + 1
                pos = n1
                if(n_proc > 1) pos = pos + off_1(this_proc)
                ic_1(pos) = c1
                ib_1(pos) = c2
                ip_1(pos) = this_proc
                xf_1(pos) = Grid(d1) % xf(s)
                yf_1(pos) = Grid(d1) % yf(s)
                zf_1(pos) = Grid(d1) % zf(s)
              end if
            end if
          end do

          ! On the side of d2
          n2 = 0
          do s = 1, Grid(d2) % n_faces
            c1 = Grid(d2) % faces_c(1,s)
            c2 = Grid(d2) % faces_c(2,s)
            if(c2 < 0) then
              if(Grid(d2) % Bnd_Cond_Name(c2) .eq. keys(2)) then
                n2 = n2 + 1
                pos = n2
                if(n_proc > 1) pos = pos + off_2(this_proc)
                ic_2(pos) = c1
                ib_2(pos) = c2
                ip_2(pos) = this_proc
                xf_2 (pos) = Grid(d2) % xf(s)
                yf_2 (pos) = Grid(d2) % yf(s)
                zf_2 (pos) = Grid(d2) % zf(s)
              end if
            end if
          end do
          do s = 1, Grid(d2) % n_faces
            c1 = Grid(d2) % faces_c(1,s)
            c2 = Grid(d2) % faces_c(2,s)
            if(Grid(d2) % faces_s(s) > 0) then  ! only if it has a shadow
              if(Grid(d2) % Bnd_Cond_Name(s) .eq. keys(2) .and.  &
                 Grid(d2) % comm % cell_proc(c1) .eq. this_proc) then
                n2 = n2 + 1
                pos = n2
                if(n_proc > 1) pos = pos + off_2(this_proc)
                ic_2(pos) = c1
                ib_2(pos) = c2
                ip_2(pos) = this_proc
                xf_2 (pos) = Grid(d2) % xf(s)
                yf_2 (pos) = Grid(d2) % yf(s)
                zf_2 (pos) = Grid(d2) % zf(s)
              end if
            end if
          end do

          !----------------------------------------------------------!
          !   Distribute interface coordinates over all processors   !
          !----------------------------------------------------------!
          call Comm_Mod_Global_Sum_Real_Array(n_tot, xf_1(1:n_tot))
          call Comm_Mod_Global_Sum_Real_Array(n_tot, yf_1(1:n_tot))
          call Comm_Mod_Global_Sum_Real_Array(n_tot, zf_1(1:n_tot))
          call Comm_Mod_Global_Sum_Real_Array(n_tot, xf_2(1:n_tot))
          call Comm_Mod_Global_Sum_Real_Array(n_tot, yf_2(1:n_tot))
          call Comm_Mod_Global_Sum_Real_Array(n_tot, zf_2(1:n_tot))
          call Comm_Mod_Global_Sum_Int_Array (n_tot, ic_1(1:n_tot))
          call Comm_Mod_Global_Sum_Int_Array (n_tot, ib_1(1:n_tot))
          call Comm_Mod_Global_Sum_Int_Array (n_tot, ip_1(1:n_tot))
          call Comm_Mod_Global_Sum_Int_Array (n_tot, ic_2(1:n_tot))
          call Comm_Mod_Global_Sum_Int_Array (n_tot, ib_2(1:n_tot))
          call Comm_Mod_Global_Sum_Int_Array (n_tot, ip_2(1:n_tot))

          ! Sort interfaces from domain 1 carrying
          ! information of cells surrounding it along
          call Sort % Three_Real_Carry_Three_Int(xf_1(1:n_tot),  &
                                                 yf_1(1:n_tot),  &
                                                 zf_1(1:n_tot),  &
                                                 ic_1(1:n_tot),  &
                                                 ib_1(1:n_tot),  &
                                                 ip_1(1:n_tot))

          ! Sort interfaces from domain 2 carrying
          ! information of cells surrounding it along
          call Sort % Three_Real_Carry_Three_Int(xf_2(1:n_tot),  &
                                                 yf_2(1:n_tot),  &
                                                 zf_2(1:n_tot),  &
                                                 ic_2(1:n_tot),  &
                                                 ib_2(1:n_tot),  &
                                                 ip_2(1:n_tot))
!         ! Write some debugging information
!         do n = 1, inter(d1, d2) % n_tot
!           write(100*this_proc+d1*10+d2, '(i4, 6f10.5)')   &
!                     n,                                    &
!                     xf_1(n), yf_1(n), zf_1(n),            &
!                     xf_2(n), yf_2(n), zf_2(n)
!         end do

          !------------------------------!
          !   Store buffer information   !
          !------------------------------!
          allocate(inter(d1, d2) % cell_1(inter(d1, d2) % n1_sub))
          allocate(inter(d1, d2) % face_1(inter(d1, d2) % n1_sub))
          allocate(inter(d1, d2) % bcel_1(inter(d1, d2) % n1_sub))
          allocate(inter(d1, d2) % cell_2(inter(d1, d2) % n2_sub))
          allocate(inter(d1, d2) % face_2(inter(d1, d2) % n2_sub))
          allocate(inter(d1, d2) % bcel_2(inter(d1, d2) % n2_sub))
          inter(d1, d2) % cell_1 = 0
          inter(d1, d2) % face_1 = 0
          inter(d1, d2) % bcel_1 = 0
          inter(d1, d2) % cell_2 = 0
          inter(d1, d2) % face_2 = 0
          inter(d1, d2) % bcel_2 = 0

          n1 = 0
          n2 = 0
          do n = 1, inter(d1, d2) % n_tot  ! browse through interface now

            ! Handle domain 1
            if(ip_1(n) .eq. this_proc) then
              n1 = n1 + 1
              do c = 1, Grid(d1) % n_cells
                if(c .eq. ic_1(n)) then
                  inter(d1, d2) % cell_1(n1) = c
                  inter(d1, d2) % face_1(n1) = n
                  inter(d1, d2) % bcel_1(n1) = ib_1(n)
                  goto 1
                end if
              end do
1             continue
            end if

            ! Handle domain 2
            if(ip_2(n) .eq. this_proc) then
              n2 = n2 + 1
              do c = 1, Grid(d2) % n_cells
                if(c .eq. ic_2(n)) then
                  inter(d1, d2) % cell_2(n2) = c
                  inter(d1, d2) % face_2(n2) = n
                  inter(d1, d2) % bcel_2(n2) = ib_2(n)
                  goto 2
                end if
              end do
2             continue
            end if

          end do  ! next face in the interface

        end if  ! if this interface is found in control mod
      end if

    end do
  end do

  if(n_proc > 1) then
    deallocate(off_1)
    deallocate(off_2)
  end if

  call Work % Disconnect_Real_Face(xf_1, yf_1, zf_1, xf_2, yf_2, zf_2)
  call Work % Disconnect_Int_Face(ic_1, ib_1, ip_1, ic_2, ib_2, ip_2)

  end subroutine
