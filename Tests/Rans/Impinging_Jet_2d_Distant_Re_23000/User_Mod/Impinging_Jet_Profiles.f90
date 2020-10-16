!==============================================================================!
  subroutine User_Mod_Impinging_Jet_Profiles(turb)
!------------------------------------------------------------------------------!
!   Subroutine reads ".1D" file created by the "Generator" or "Convert"        !
!   and extracts profiles on several locations that corresponds with the       !
!   experimental measurements.                                                 !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Turb_Type), target :: turb
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),  pointer :: grid
  type(Field_Type), pointer :: flow
  type(Var_Type),   pointer :: u, v, w, t
  type(Var_Type),   pointer :: kin, eps, zeta, f22
  integer                   :: n_prob, pl, i, count, s, c, idumm, k, fu
  character(len=80)         :: coord_name, res_name, store_name, ext
  real, allocatable         :: z_p(:),                              &
                               um_p(:), vm_p(:), wm_p(:), tm_p(:),  &
                               v1_p(:), v2_p(:), v3_p(:),           &
                               v4_p(:), v5_p(:), v6_p(:),           &
                               zm_p(:), rad_1(:),                   &
                               ind(:)
  integer, allocatable      :: n_p(:), n_count(:)
  real                      :: r, r1, r2, u_aver, u_rad, u_tan, lnum
  logical                   :: there
!==============================================================================!

  ! Take aliases
  flow => turb % pnt_flow
  grid => flow % pnt_grid
  u    => flow % u
  v    => flow % v
  w    => flow % w
  t    => flow % t
  kin  => turb % kin
  eps  => turb % eps
  zeta => turb % zeta
  f22  => turb % f22

  u_aver = 1.14

  inquire(file='rad_coordinate.dat', exist=there ) 
  if(.not.there) then
    if(this_proc < 2) then
      print *, "#==========================================================="
      print *, "# In order to extract profiles and write them in ascii files"
      print *, "# the code has to read cell-faces coordinates "
      print *, "# in wall-normal direction in the ascii file 'case_name'.1D."
      print *, "# The file format should be as follows:"
      print *, "# 10  ! number of cells + 1"
      print *, "# 1 0.0"
      print *, "# 2 0.1"
      print *, "# 3 0.2"
      print *, "# ... "
      print *, "#-----------------------------------------------------------"
    end if
    return
  end if

  ! Set the name for coordinate file
  call File_Mod_Set_Name(coord_name, extension='.1d')
  call File_Mod_Open_File_For_Writing(coord_name, fu)

  ! Read the number of searching intervals
  read(fu,*) n_prob
  allocate(z_p(n_prob*2))
  allocate(ind(n_prob*2))

  ! Read the intervals positions
  do pl=1,n_prob
    read(fu,*) ind(pl), z_p(pl)
  end do
  close(fu)
!  call Sort_Mod_Real_Carry_Int(z_p(1:n_prob), ind(1:n_prob))

  allocate(n_p(n_prob));   n_p  = 0 
  allocate(um_p(n_prob));  um_p = 0.0
  allocate(vm_p(n_prob));  vm_p = 0.0
  allocate(wm_p(n_prob));  wm_p = 0.0
  allocate(v1_p(n_prob));  v1_p = 0.0
  allocate(v2_p(n_prob));  v2_p = 0.0
  allocate(v3_p(n_prob));  v3_p = 0.0
  allocate(v4_p(n_prob));  v4_p = 0.0
  allocate(v5_p(n_prob));  v5_p = 0.0
  allocate(zm_p(n_prob));  zm_p = 0.0

  allocate(n_count(n_prob)); n_count=0
  count = 0

  if(heat_transfer) then
    allocate(tm_p(n_prob));   tm_p = 0.0
  end if

  !-------------------------!
  !   Average the results   !
  !-------------------------!
  do k = 0, 6
    if(k .eq. 0) then
      r1 = 0.0
      r2 = 0.04
      lnum = 0.0
    else if(k .eq. 1) then
      r1 = 0.992
      r2 = 1.0
      lnum = 0.5
    else if(k .eq. 2) then
      r1 = 2.0
      r2 = 2.1500
      lnum = 1.0
    else if(k .eq. 3) then
      r1 = 2.9744
      r2 = 3.0684
      lnum = 1.5
    else if(k .eq. 4) then
      r1 = 3.9098
      r2 = 4.1433
      lnum = 2.0
    else if(k .eq. 5) then
      r1 = 0.4803200E+01
      r2 = 0.5347000E+01
      lnum = 2.5
    else if(k .eq. 6) then
      r1 = 0.5876600E+01
      r2 = 0.6000000E+01
      lnum = 3.0
    end if

    do i = 1, n_prob-1
      do c = 1, grid % n_cells
        r = sqrt(grid % xc(c)**2 + grid % yc(c)**2) + TINY
        if(r > r1 .and. r < r2) then
          if(grid % zc(c) > z_p(i) .and.  &
             grid % zc(c) < z_p(i+1)) then
            u_rad   = (  u % n(c)*grid % xc(c)/r   &
                       + v % n(c)*grid % yc(c)/r)
            u_tan   = (- u % n(c)*grid % yc(c)/r   &
                       + v % n(c)*grid % xc(c)/r)
            um_p(i)   = um_p(i) + sqrt(  u % n(c)**2   &
                                       + v % n(c)**2   &
                                       + w % n(c)**2)
            vm_p(i)   = vm_p(i) + u_rad
            wm_p(i)   = wm_p(i) + w % n(c)

            if(turb % model .eq. K_EPS) then
              v1_p(i) = v1_p(i) + kin % n(c)
              v2_p(i) = v2_p(i) + eps % n(c)
              v3_p(i) = v3_p(i) + turb % vis_t(c) / flow % viscosity(c)
            end if

            if(turb % model .eq. K_EPS_ZETA_F) then
              v1_p(i)   = v1_p(i) + kin % n(c)
              v2_p(i)   = v2_p(i) + eps % n(c)
              v3_p(i)   = v3_p(i) + turb % vis_t(c) / flow % viscosity(c)
              v4_p(i)   = v4_p(i) + zeta % n(c)
              v5_p(i)   = v5_p(i) + f22 % n(c)
            end if

            if(heat_transfer) then
              tm_p(i)   = tm_p(i) + T % n(c)
            end if

            zm_p(i) = zm_p(i) + grid % zc(c)
            n_count(i) = n_count(i) + 1
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
      call Comm_Mod_Global_Sum_Real(zm_p(pl))

      call Comm_Mod_Global_Sum_Real(v1_p(pl))
      call Comm_Mod_Global_Sum_Real(v2_p(pl))
      call Comm_Mod_Global_Sum_Real(v3_p(pl))
      call Comm_Mod_Global_Sum_Real(v4_p(pl))
      call Comm_Mod_Global_Sum_Real(v5_p(pl))

      count = count + n_count(pl) 

      if(heat_transfer) then
        call Comm_Mod_Global_Sum_Real(tm_p(pl))
      end if
    end do

    do i = 1, n_prob
      if(n_count(i) .ne. 0) then
        wm_p(i) = wm_p(i)/n_count(i)
        um_p(i) = um_p(i)/n_count(i)
        vm_p(i) = vm_p(i)/n_count(i)
        v1_p(i) = v1_p(i)/n_count(i)
        v2_p(i) = v2_p(i)/n_count(i)
        v3_p(i) = v3_p(i)/n_count(i)
        v4_p(i) = v4_p(i)/n_count(i)
        v5_p(i) = v5_p(i)/n_count(i)
        tm_p(i) = tm_p(i)/n_count(i)
        zm_p(i) = zm_p(i)/n_count(i)
      end if
    end do

    if(this_proc < 2) then

      ! Set the file name
      call File_Mod_Set_Name(res_name, appendix='-prof', extension='.dat')
      call File_Mod_Open_File_For_Writing(res_name, fu)

      write(fu,'(a1,2x,a101)') '#', ' 1:Xrad, ' // &
                                    ' 2:Umag, ' // &
                                    ' 3:Urad, ' // &
                                    ' 4:Uaxi, ' // &
                                    ' 5:Kin,  ' // &
                                    ' 6:Eps,  ' // &
                                    ' 7:Temp, ' // &
                                    ' 8:vis_t/viscosity, '//  &
                                    ' 9:zeta, ' // &
                                    '10:f22   '

      do i = 1, n_prob
        if(n_count(i) .ne. 0) then
          write(fu,'(10e11.3)') zm_p(i)/2.0,           &
                                um_p(i)/u_aver,        &
                                vm_p(i)/u_aver,        &
                                wm_p(i)/u_aver,        &
                                v1_p(i)/u_aver**2,     &
                                v2_p(i),               &
                                tm_p(i),               &
                                v3_p(i),               &
                                v4_p(i),               &
                                v5_p(i)
        end if
      end do
      close(fu)
    end if

    do i = 1, n_prob
      n_count(i) = 0
      wm_p(i)    = 0.0
      um_p(i)    = 0.0
      vm_p(i)    = 0.0
      v1_p(i)    = 0.0
      v2_p(i)    = 0.0
      v3_p(i)    = 0.0
      v4_p(i)    = 0.0
      v5_p(i)    = 0.0
      tm_p(i)    = 0.0
      zm_p(i)    = 0.0
    end do
    if(this_proc < 2) print *, 'Finished with profile r/D =  ', lnum
  end do   ! end number of radius

  deallocate(n_p)
  deallocate(z_p)
  deallocate(um_p)
  deallocate(vm_p)
  deallocate(wm_p)
  deallocate(v1_p)
  deallocate(v2_p)
  deallocate(v3_p)
  deallocate(v4_p)
  deallocate(v5_p)
  deallocate(zm_p)
  deallocate(n_count)
  if(heat_transfer) then
    deallocate(tm_p)
  end if

  if(this_proc < 2) write(*,*) 'Finished with User_Impinging_Jet_Profiles'

  end subroutine 
