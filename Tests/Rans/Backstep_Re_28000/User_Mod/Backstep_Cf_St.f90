!==============================================================================!
  subroutine User_Mod_Backstep_Cf_St(flow, turb, save_name)
!------------------------------------------------------------------------------!
!   Subroutine extracts skin friction coefficient and Stanton number for       !
!   backstep case.                                                             !
!------------------------------------------------------------------------------!
  use Grid_Mod,  only: Grid_Type
  use Field_Mod, only: Field_Type,  &
                       viscosity, capacity, density, conductivity, heat_transfer
  use Var_Mod,   only: Var_Type
  use Turb_Mod
  use Comm_Mod                       ! parallel stuff
  use Name_Mod,  only: problem_name
  use Const_Mod                      ! constants
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type), target :: flow
  type(Turb_Type),  target :: turb
  character(len=*)         :: save_name
!----------------------------------[Calling]-----------------------------------!
  real :: Y_Plus_Low_Re
!-----------------------------------[Locals]-----------------------------------!
  type(Var_Type),  pointer :: u, v, w, t
  type(Grid_Type), pointer :: grid
  integer                  :: n_prob, pl, c, dummy, i, count, k, c1, c2, s, l
  character(len=80)        :: result_name
  real, allocatable        :: r1_p(:), r2_p(:), z_p(:),  &
                              um_p(:), vm_p(:), wm_p(:), & 
                              uu_p(:), vv_p(:), ww_p(:), &
                              uv_p(:), uw_p(:), vw_p(:), &
                              tm_p(:), tt_p(:),          &
                              ut_p(:), vt_p(:), wt_p(:), &
                              v1_p(:), v2_p(:), v3_p(:), &
                              v4_p(:), v5_p(:)  
  integer, allocatable     :: n_p(:), n_count(:)
  real                     :: kin_vis, u_tan, u_tau, tau_wall
!==============================================================================!

  ! Take aliases
  grid => flow % pnt_grid
  u    => flow % u
  v    => flow % v
  w    => flow % w
  t    => flow % t

  !----------------------------------!
  !   Read "x_coordinate.dat" file   !
  !----------------------------------!
  if(this_proc < 2) print *, '# Now reading the file: x_coordinate.dat ' 
  open(9, file='x_coordinate.dat')

  ! Read the number of probes 
  read(9,*) n_prob
  allocate(z_p(n_prob))

  ! Read the probe coordinates 
  do pl=1,n_prob
    read(9,*) z_p(pl)
  end do
  close(9)

  allocate(n_p(n_prob));   n_p  = 0 
  allocate(um_p(n_prob));  um_p = 0.0
  allocate(vm_p(n_prob));  vm_p = 0.0
  allocate(wm_p(n_prob));  wm_p = 0.0
  allocate(uu_p(n_prob));  uu_p = 0.0
  allocate(vv_p(n_prob));  vv_p = 0.0
  allocate(ww_p(n_prob));  ww_p = 0.0
  allocate(uv_p(n_prob));  uv_p = 0.0
  allocate(uw_p(n_prob));  uw_p = 0.0
  allocate(vw_p(n_prob));  vw_p = 0.0
  allocate(v1_p(n_prob));  v1_p = 0.0
  allocate(v2_p(n_prob));  v2_p = 0.0
  allocate(v3_p(n_prob));  v3_p = 0.0
  allocate(v4_p(n_prob));  v4_p = 0.0
  allocate(v5_p(n_prob));  v5_p = 0.0

  allocate(n_count(n_prob)); n_count = 0
  count = 0

  if(heat_transfer) then
    allocate(tm_p(n_prob));  tm_p = 0.0
    allocate(tt_p(n_prob));  tt_p = 0.0
    allocate(ut_p(n_prob));  ut_p = 0.0
    allocate(vt_p(n_prob));  vt_p = 0.0
    allocate(wt_p(n_prob));  wt_p = 0.0
  end if  

  !-------------------------!
  !   Average the results   !
  !-------------------------!
  do i = 1, n_prob-1
    do s = 1, grid % n_faces
      c1 = grid % faces_c(1,s)
      c2 = grid % faces_c(2,s)
      if(c2 < 0) then
        if(Grid_Mod_Bnd_Cond_Type(grid,c2).eq.WALLFL.and.t % q(c2) > 1.e-8) then
          if(grid % xc(c1) > z_p(i) .and. grid % xc(c1) < z_p(i+1)) then
            um_p(i)   = um_p(i) + u % n(c1)
            vm_p(i)   = vm_p(i) + v % n(c1)
            wm_p(i)   = wm_p(i) + w % n(c1)
            if(turb % y_plus(c1) < 4.0) then
              v1_p(i) = v1_p(i)  &
                      + (2.0 * viscosity * u % n(c1)   &
                             / grid % wall_dist(c1))   &
                      / (density * 11.3**2)
            else
              kin_vis = viscosity / density
              u_tan = Field_Mod_U_Tan(flow, s)
              u_tau = c_mu25 * sqrt(turb % kin % n(c1))
              turb % y_plus(c1) = Y_Plus_Low_Re(u_tau,                 &
                                                grid % wall_dist(c1),  &
                                                kin_vis)
              tau_wall = density * kappa * u_tau * u_tan    &
                       / log(e_log*max(turb % y_plus(c1), 1.05))

              v1_p(i) = v1_p(i)  &
                      + 0.015663 * tau_wall * u % n(c1) / abs(u % n(c1))
            end if
            v2_p(i) = v2_p(i) + turb % y_plus(c1)
            v3_p(i) = v3_p(i) + t % q(c2)  &
                    / (density * capacity * (t % n(c2) - 20) * 11.3)
            v5_p(i) = v5_p(i) + t % n(c2) 
            n_count(i) = n_count(i) + 1
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
  write(result_name(l+1:l+10),'(a10)') '-cf-st.dat'

  open(3,file=result_name)

  write(3,*) '# x, Cf, St, U, T, yPlus'
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

      write(3,'(6es15.5e3)') (z_p(i)+z_p(i+1))/(2.*0.038),  &
                             v1_p(i),                       &
                             v3_p(i),                       &
                             um_p(i),                       &
                             v5_p(i),                       &
                             v2_p(i) 
    end if
  end do 
  close(3)

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

  if(this_proc < 2) write(*,*) '# Finished with User_Backstep_Cf_St'

  end subroutine
