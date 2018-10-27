!==============================================================================!
  subroutine User_Mod_Impinging_Jet_Nu(grid, save_name)     
!------------------------------------------------------------------------------!
!   The subroutine creates ASCII file with Nusselt number averaged             !
!   in azimuthal direction.                                                    !    
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
!-----------------------------------[Locals]-----------------------------------!
  integer             :: n_prob, pl, i, count, s, c1, c2
  character(len=80)   :: res_name, store_name
  real, allocatable   :: z_p(:),                              &
                         um_p(:), vm_p(:), wm_p(:), tm_p(:),  &
                         uu_p(:), vv_p(:), ww_p(:),           &
                         uv_p(:), uw_p(:), vw_p(:),           &
                         v1_p(:), v2_p(:), v3_p(:),           &
                         v4_p(:), v5_p(:), v6_p(:),           &
                         rm_p(:), rad_1(:),                   &
                         ind(:)
  integer,allocatable :: n_p(:), n_count(:)
  real                :: r
  logical             :: there
!==============================================================================!

  inquire( file='rad_coordinate.dat', exist=there ) 
  if(.not.there) then
    if(this_proc < 2) then
      print *, "#=========================================================="
      print *, "# In order to extract Nusselt number profile               "
      print *, "# an ascii file with cell-faces coordinates has to be read."
      print *, "# The name of the file is rad_coordinate.dat.              "
      print *, "# The file format should be as follows:                    "
      print *, "# 10  ! number of cells + 1                                "
      print *, "# 0.0                                                      "
      print *, "# 0.1                                                      "
      print *, "# 0.2                                                      "
      print *, "# ...                                                      "
      print *, "#----------------------------------------------------------"
    end if
    return
  end if

  open(9, file='rad_coordinate.dat')

  ! Read the number of searching intervals 
  read(9,*) n_prob
  allocate(z_p(n_prob*2))
  allocate(ind(n_prob*2))

  ! Read the intervals positions
  do pl=1,n_prob
    read(9,*) ind(pl), z_p(pl) 
  end do
  close(9)

  allocate(n_p(n_prob));   n_p  = 0
  allocate(um_p(n_prob));  um_p = 0.0
  allocate(vm_p(n_prob));  vm_p = 0.0
  allocate(wm_p(n_prob));  wm_p = 0.0
  allocate(rm_p(n_prob));  rm_p = 0.0
  allocate(v1_p(n_prob));  v1_p = 0.0
  allocate(v2_p(n_prob));  v2_p = 0.0
  allocate(v3_p(n_prob));  v3_p = 0.0
  allocate(v4_p(n_prob));  v4_p = 0.0
  allocate(v5_p(n_prob));  v5_p = 0.0
  allocate(v6_p(n_prob));  v6_p = 0.0

  allocate(rad_1(n_prob));   rad_1   = 0.0
  allocate(n_count(n_prob)); n_count = 0

  count = 0

  if(heat_transfer) then
    allocate(tm_p(n_prob));   tm_p=0.0
  end if  

  !-------------------------!
  !   Average the results   !
  !-------------------------!
  do i = 1, n_prob
    rad_1(i) = abs(z_p(i))
  end do

  do i = 1, n_prob-1
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 < 0) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2) .eq. WALLFL) then
          r = sqrt(grid % xc(c1)*grid % xc(c1)  + &
                   grid % yc(c1)*grid % yc(c1)) + TINY
          if(r < rad_1(i+1) .and.  &
             r > rad_1(i)   .and.  &
             grid % zc(c1) < 0.5) then
            rm_p(i) = rm_p(i) + sqrt(grid % xc(c1)*grid % xc(c1)  + &
                                     grid % yc(c1)*grid % yc(c1))
            um_p(i) = um_p(i) +   u % n(c1) * grid % xc(c1) / r  + &
                                  v % n(c1) * grid % yc(c1) / r
            vm_p(i) = vm_p(i) + (-u % n(c1) * grid % yc(c1) / r  + &
                                  v % n(c1) * grid % xc(c1) / r)
            wm_p(i) = wm_p(i) +   w % n(c1)
            tm_p(i) = tm_p(i) + t % n(c2) 
            v1_p(i) = v1_p(i) + grid % zc(c1)
            v2_p(i) = v2_p(i) + sqrt(tau_wall(c1))
            v3_p(i) = v3_p(i) + (c_mu**0.25 * kin % n(c1)**0.5) 
            v4_p(i) = v4_p(i) + kin % n(c1)
            v5_p(i) = v5_p(i) + eps % n(c1)
            v6_p(i) = v6_p(i) + t % q(c2)
            n_count(i)= n_count(i) + 1
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

    call Comm_Mod_Global_Sum_Real(v1_p(pl))
    call Comm_Mod_Global_Sum_Real(v2_p(pl))
    call Comm_Mod_Global_Sum_Real(v3_p(pl))
    call Comm_Mod_Global_Sum_Real(v4_p(pl))
    call Comm_Mod_Global_Sum_Real(v5_p(pl))
    call Comm_Mod_Global_Sum_Real(v6_p(pl))

    call Comm_Mod_Global_Sum_Real(rm_p(pl))
    call Comm_Mod_Global_Sum_Real(tm_p(pl))

    count = count + n_count(pl)
  end do

  do i = 1, n_prob
    if(n_count(i) .ne. 0) then
      wm_p(i) = wm_p(i) / n_count(i)
      vm_p(i) = vm_p(i) / n_count(i)
      um_p(i) = um_p(i) / n_count(i)
      tm_p(i) = tm_p(i) / n_count(i)
      v1_p(i) = v1_p(i) / n_count(i)
      v2_p(i) = v2_p(i) / n_count(i)
      v3_p(i) = v3_p(i) / n_count(i)
      v4_p(i) = v4_p(i) / n_count(i)
      v5_p(i) = v5_p(i) / n_count(i)
      v6_p(i) = v6_p(i) / n_count(i)
      rm_p(i) = rm_p(i) / n_count(i)
    end if
  end do
  call Comm_Mod_Wait

    ! Set the file name
    store_name = problem_name
    problem_name = save_name
    call Name_File(0, res_name, '-nu-utau.dat')
    open(3, file = res_name)
    problem_name = store_name

  write(3,*) '# Xrad, Nu, Utau, Yplus, Temp, Numb of points '
  do i = 1, n_prob
    if(n_count(i) .ne. 0) then
      write(3,'(5e11.3,i6)') rm_p(i)/2.0,                                &
                             2.0*v6_p(i)/(conductivity*(tm_p(i)-20.0)),  & 
                             v2_p(i),                                    &
                             v2_p(i) * v1_p(i)/viscosity,                &
                             tm_p(i),                                    &
                             n_count(i)
    end if
  end do
  close(3)

  deallocate(n_p)
  deallocate(um_p)
  deallocate(vm_p)
  deallocate(wm_p)

  deallocate(v1_p)
  deallocate(v2_p)
  deallocate(v3_p)
  deallocate(v4_p)
  deallocate(v5_p)
  deallocate(v6_p)
  deallocate(rm_p)
  deallocate(rad_1)
  deallocate(n_count)

  if(heat_transfer) then
    deallocate(tm_p)
  end if

  if(this_proc < 2) print *, '# Finished with Impinging_Jet_Nu'

  end subroutine
