!==============================================================================!
  subroutine User_Mod_Backstep_Profiles(grid, save_name) 
!------------------------------------------------------------------------------!
!   Description
!------------------------------------------------------------------------------!
  use Grid_Mod
  use Flow_Mod
  use Rans_Mod
  use Comm_Mod                       ! parallel stuff
  use Name_Mod,  only: problem_name
  use Const_Mod                      ! constants
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Grid_Type)  :: grid
  character(len=*) :: save_name
!----------------------------------[Calling]-----------------------------------!
  integer              :: n_prob, pl, c, idumm, i, count, k, c1, c2, s, n_hor, l
  character(len=80)    :: coord_name, result_name
  real, allocatable    :: x1_p(:), x2_p(:), lnum(:), z_p(:), &
                          um_p(:), vm_p(:), wm_p(:), & 
                          uu_p(:), vv_p(:), ww_p(:), &
                          uv_p(:), uw_p(:), vw_p(:), &
                          tm_p(:), tt_p(:),          &
                          ut_p(:), vt_p(:), wt_p(:), &
                          v1_p(:), v2_p(:), v3_p(:), &
                          v4_p(:), v5_p(:)
  integer, allocatable :: n_p(:), n_count(:)
  real                 :: z_coor
  logical              :: there
!==============================================================================!

  ! Set the name for coordinate file
  call Name_File(0, coord_name, ".1d")

  inquire( file='horizontal_positions.dat', exist=there ) 

  if(.not.there) then
    if(this_proc < 2) then
      print *, "#============================================================"
      print *, "# In order to extract profiles and write them in ascii files "
      print *, "# the code has to read ascii file in which positions of the  "
      print *, "# profiles extraction are defined. Name of file is           "
      print *, "# 'horizontal_positions.dat.'                                "
      print *, "# The file format should be as follows:                      "
      print *, "# 1  ! number of positions                                   "
      print *, "# 0.03 0.02 0.35 first two numbers defined the range in the  "
      print *, "#                streamwise direction in which the profile   "
      print *, "#                will be extracted and the last number       "
      print *, "#                represents x/h, where h = 0.038, which will "
      print *, "#                appears in the file name.                   "
      print *, "#------------------------------------------------------------"
    end if
    return
  end if

  open(9, file='horizontal_positions.dat ')
  read(9,*) n_hor
  allocate(x1_p(n_hor))
  allocate(x2_p(n_hor))
  allocate(lnum(n_hor))
  do pl=1, n_hor
    read(9,*) x1_p(pl), x2_p(pl), lnum(pl)
  end do
  close(9)

  !------------------!
  !   Read 1d file   !
  !------------------!
  if(this_proc < 2) print *, '# Now reading the file:', coord_name

  inquire( file=coord_name, exist=there ) 
  if(.not.there) then
    if(this_proc < 2) then
      print *, "==========================================================="
      print *, "In order to extract profiles and write them in ascii files "
      print *, "the code has to read cell-faces coordinates                "
      print *, "in wall-normal direction in the ascii file 'case_name'.1D. "
      print *, "The file format should be as follows:                      "
      print *, "N   ! N is number of cells + 1                             "
      print *, "1 0.0                                                      "
      print *, "2 0.1                                                      "
      print *, "3 0.2                                                      "
      print *, "...                                                        "
      print *, "==========================================================="
    end if
    return
  end if

  open(9, file=coord_name)

  ! Read the number of probes 
  read(9,*) n_prob
  allocate(z_p(n_prob))

  ! Read the probe coordinates out
  do pl=1,n_prob
    read(9,*) idumm, z_p(pl)
  end do
  close(9)
  
  allocate(n_p(n_prob));     n_p     = 0 
  allocate(um_p(n_prob));    um_p    = 0.0
  allocate(vm_p(n_prob));    vm_p    = 0.0
  allocate(wm_p(n_prob));    wm_p    = 0.0
  allocate(uu_p(n_prob));    uu_p    = 0.0
  allocate(vv_p(n_prob));    vv_p    = 0.0
  allocate(ww_p(n_prob));    ww_p    = 0.0
  allocate(uv_p(n_prob));    uv_p    = 0.0
  allocate(uw_p(n_prob));    uw_p    = 0.0
  allocate(vw_p(n_prob));    vw_p    = 0.0
  allocate(v1_p(n_prob));    v1_p    = 0.0
  allocate(v2_p(n_prob));    v2_p    = 0.0
  allocate(v3_p(n_prob));    v3_p    = 0.0
  allocate(v4_p(n_prob));    v4_p    = 0.0
  allocate(v5_p(n_prob));    v5_p    = 0.0
  allocate(n_count(n_prob)); n_count = 0
  count = 0

  if(heat_transfer) then
    allocate(tm_p(n_prob));   tm_p = 0.0
    allocate(tt_p(n_prob));   tt_p = 0.0
    allocate(ut_p(n_prob));   ut_p = 0.0
    allocate(vt_p(n_prob));   vt_p = 0.0
    allocate(wt_p(n_prob));   wt_p = 0.0
  end if  

  !-------------------------!
  !   Average the results   !
  !-------------------------!
  do k = 1, n_hor
    do i = 1, n_prob-1
      do c = 1, grid % n_cells
        z_coor = grid % zc(c)
        if(grid % xc(c) < x1_p(k) .and. grid % xc(c) > x2_p(k)) then
          if(z_coor > z_p(i) .and. z_coor < z_p(i+1)) then
            um_p(i) = um_p(i) + u % n(c)
            vm_p(i) = vm_p(i) + v % n(c)
            wm_p(i) = wm_p(i) + w % n(c)
            uu_p(i) = uu_p(i) + kin % n(c)
            vv_p(i) = vv_p(i) + eps % n(c)
            v1_p(i) = v1_p(i) + vis_t(c)*(U % y(c) + V % x(c))/11.3**2 
            v2_p(i) = v2_p(i) + t % n(c) - 20.0
            n_count(i) = n_count(i) + 1
            if(turbulence_model == K_EPS_ZETA_F) then      
              ww_p(i) = ww_p(i) + zeta % n(c)
              uv_p(i) = uv_p(i) + f22 % n(c)
            end if
          end if
        end if
      end do 
    end do 

    ! Average over all processors
    do pl=1, n_prob
      call Comm_Mod_Global_Sum_Int(n_count(pl))
      call Comm_Mod_Global_Sum_Real(um_p(pl))
      call Comm_Mod_Global_Sum_Real(vm_p(pl))
      call Comm_Mod_Global_Sum_Real(wm_p(pl))
      call Comm_Mod_Global_Sum_Real(uu_p(pl))
      call Comm_Mod_Global_Sum_Real(vv_p(pl))
      call Comm_Mod_Global_Sum_Real(ww_p(pl))
      call Comm_Mod_Global_Sum_Real(uv_p(pl))
      call Comm_Mod_Global_Sum_Real(uw_p(pl))
      call Comm_Mod_Global_Sum_Real(vw_p(pl))
      call Comm_Mod_Global_Sum_Real(v1_p(pl))
      call Comm_Mod_Global_Sum_Real(v2_p(pl))
      call Comm_Mod_Global_Sum_Real(v3_p(pl))
      call Comm_Mod_Global_Sum_Real(v4_p(pl))
      call Comm_Mod_Global_Sum_Real(v5_p(pl))

      count =  count + n_count(pl) 

      if(heat_transfer) then
        call Comm_Mod_Global_Sum_Real(tm_p(pl))
        call Comm_Mod_Global_Sum_Real(tt_p(pl))
        call Comm_Mod_Global_Sum_Real(ut_p(pl))
        call Comm_Mod_Global_Sum_Real(vt_p(pl))
        call Comm_Mod_Global_Sum_Real(wt_p(pl))
      end if
    end do

    result_name = save_name ! problem_name
    l = len_trim(result_name)
    write(result_name(l+1:l+14),'(a5,f4.2,a5)') '-prof', lnum(k), 'h.dat'

    open(3,file=result_name)
    write(3,*) '# z, U, TKE, EPS, zeta, f, uv, T'
    write(3,*) '# all data are normalized by Ub = 11.3 and step height h = 0.038'
    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        wm_p(i) = wm_p(i) / n_count(i)
        um_p(i) = um_p(i) / n_count(i)
        vm_p(i) = vm_p(i) / n_count(i)
        uu_p(i) = uu_p(i) / n_count(i)
        vv_p(i) = vv_p(i) / n_count(i)
        ww_p(i) = ww_p(i) / n_count(i)
        uv_p(i) = uv_p(i) / n_count(i)
        uw_p(i) = uw_p(i) / n_count(i)
        vw_p(i) = vw_p(i) / n_count(i)
        v1_p(i) = v1_p(i) / n_count(i)
        v2_p(i) = v2_p(i) / n_count(i)
        v3_p(i) = v3_p(i) / n_count(i)
        v4_p(i) = v4_p(i) / n_count(i)
        v5_p(i) = v5_p(i) / n_count(i)

        write(3,'(8E15.7)') (z_p(i)+z_p(i+1))/(2.*0.038),  &
                           um_p(i)/11.3,                   &
                           uu_p(i)/11.3**2,                &
                           vv_p(i)*0.038/11.3**3,          &
                           ww_p(i)/11.3**2,                &
                           uv_p(i)*0.038/11.3**2,          &
                           v1_p(i),                        &
                           v2_p(i) 

        wm_p(i)    = 0.0   
        um_p(i)    = 0.0 
        vm_p(i)    = 0.0 
        uu_p(i)    = 0.0 
        vv_p(i)    = 0.0 
        ww_p(i)    = 0.0 
        uv_p(i)    = 0.0 
        uw_p(i)    = 0.0 
        vw_p(i)    = 0.0 
        v1_p(i)    = 0.0 
        v2_p(i)    = 0.0 
        v3_p(i)    = 0.0
        v4_p(i)    = 0.0
        v5_p(i)    = 0.0
        n_count(i) = 0
      end if
    end do 
    close(3)
  end do
 
  deallocate(n_p)
  deallocate(z_p)
  deallocate(um_p)
  deallocate(vm_p)
  deallocate(wm_p)
  deallocate(uu_p)
  deallocate(vv_p)
  deallocate(ww_p)
  deallocate(uv_p)
  deallocate(uw_p)
  deallocate(vw_p)
  deallocate(v1_p)
  deallocate(v2_p)
  deallocate(v3_p)
  deallocate(v4_p)
  deallocate(v5_p)
  deallocate(n_count)
  if(heat_transfer) then
    deallocate(tm_p)
    deallocate(tt_p)
    deallocate(ut_p)
    deallocate(vt_p)
    deallocate(wt_p)
  end if

  if(this_proc < 2) print *, '# Finished with User_Backstep_Profiles'

  end subroutine
