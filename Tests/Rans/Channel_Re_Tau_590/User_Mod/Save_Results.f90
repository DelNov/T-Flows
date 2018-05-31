!==============================================================================!
  subroutine User_Mod_Save_Results(grid, save_name) 
!------------------------------------------------------------------------------!
!   This subroutine reads name.1d file created by Convert or Generator and     !
!   averages the results in homogeneous directions.                            !
!                                                                              !
!   The results are then writen in files name_res.dat and name_res_plus.dat    !
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
  integer             :: n_prob, pl, c, i, count, s, c1, c2
  character(len=80)   :: coord_name, res_name, res_name_plus
  character(len=80)   :: store_name
  real,allocatable    :: z_p(:),                    &  ! probe coordinates
                         um_p(:), vm_p(:), wm_p(:), &  ! mean velocties 
                         uu_p(:), vv_p(:), ww_p(:), &  ! stresses
                         uv_p(:), uw_p(:), vw_p(:), &  ! stresses
                         tm_p(:), tt_p(:),          &  ! mean temp and fluct.
                         ut_p(:), vt_p(:), wt_p(:), &  ! heat fluxes
                         v1_p(:), v2_p(:), v3_p(:), &  ! helping variables
                         ind(:),                    & 
                         wall_p(:), u_tau_p(:)
  integer,allocatable :: n_p(:), n_count(:)
  real                :: length_scale, t_wall, t_tau, d_wall, nu_max, u_tau_max
  logical             :: there
!==============================================================================!

  ! Set the name for coordinate file
  call Name_File(0, coord_name, ".1d")

  ! Store the name
  store_name = problem_name
  problem_name = save_name

  call Name_File(0, res_name,      "-res.dat")
  call Name_File(0, res_name_plus, "-res-plus.dat")

  !------------------!
  !   Read 1d file   !
  !------------------!
  inquire(file=coord_name, exist=there) 
  if(.not. there) then
    if(this_proc < 2) then
      print *, '==============================================================='
      print *, 'In order to extract profiles and write them in ascii files'
      print *, 'the code has to read cell-faces coordinates '
      print *, 'in wall-normal direction in the ascii file ''case_name.1d.'''
      print *, 'The file format should be as follows:'
      print *, '10  ! number of cells + 1'
      print *, '1 0.0'
      print *, '2 0.1'
      print *, '3 0.2'
      print *, '... '
      print *, '==============================================================='
    end if
    return
  end if

  open(9, file=coord_name)

  ! Write the number of searching intervals 
  read(9,*) n_prob
  allocate(z_p(n_prob*2))
  allocate(ind(n_prob*2))

  ! Read the intervals positions
  do pl=1,n_prob
    read(9,*) ind(pl), z_p(pl) 
  end do
  close(9)

  ! Sort z_p and carry ind along
  call Sort_Real_Carry_Int(z_p, ind, n_prob, 2)
   
  allocate(n_p   (n_prob));  n_p     = 0 
  allocate(wall_p(n_prob));  wall_p  = 0.0
  allocate(um_p  (n_prob));  um_p    = 0.0
  allocate(vm_p  (n_prob));  vm_p    = 0.0
  allocate(wm_p  (n_prob));  wm_p    = 0.0
  allocate(uu_p  (n_prob));  uu_p    = 0.0
  allocate(vv_p  (n_prob));  vv_p    = 0.0
  allocate(ww_p  (n_prob));  ww_p    = 0.0
  allocate(uv_p  (n_prob));  uv_p    = 0.0
  allocate(uw_p  (n_prob));  uw_p    = 0.0
  allocate(vw_p  (n_prob));  vw_p    = 0.0
  allocate(v1_p  (n_prob));  v1_p    = 0.0
  allocate(v2_p  (n_prob));  v2_p    = 0.0
  allocate(v3_p  (n_prob));  v3_p    = 0.0
  allocate(u_tau_p(n_prob)); u_tau_p = 0.0

  allocate(n_count(n_prob)); n_count=0
  count = 0
  if(heat_transfer == YES) then
    allocate(tm_p(n_prob));  tm_p = 0.0
    allocate(tt_p(n_prob));  tt_p = 0.0
    allocate(ut_p(n_prob));  ut_p = 0.0
    allocate(vt_p(n_prob));  vt_p = 0.0
    allocate(wt_p(n_prob));  wt_p = 0.0
  end if  

  length_scale = 0.0

  !-------------------------!
  !   Average the results   !
  !-------------------------!
  do i = 1, n_prob-1
    do c=1, grid % n_cells
      length_scale = max(length_scale, grid % wall_dist(c))
      if(grid % zc(c) > (z_p(i)) .and.  &
         grid % zc(c) < (z_p(i+1))) then
!        if(ROT == YES) then
!          wall_p(i) = wall_p(i) + grid % yc(c) 
!        else
          wall_p(i) = wall_p(i) + grid % wall_dist(c) 
!        end if 
        if(turbulence_model == LES      .or.  &
           turbulence_model == DNS      .or.  &
           turbulence_model == HYBRID_K_EPS_ZETA_F .or.  &
           turbulence_model == DES_SPALART) then       
          um_p(i)   = um_p(i) + u % mean(c)
          vm_p(i)   = vm_p(i) + v % mean(c)
          wm_p(i)   = wm_p(i) + w % mean(c)
        else
          um_p(i)   = um_p(i) + u % n(c)
          vm_p(i)   = vm_p(i) + v % n(c)
          wm_p(i)   = wm_p(i) + w % n(c)
        end if

        if(turbulence_model == K_EPS_ZETA_F .or.  &
           turbulence_model == K_EPS) then
          if(grid % cell_near_wall(c)) then
            if(y_plus(c) > 5.0) then        
              u_tau_p(i) = u_tau_p(i) +                       &
                           sqrt(max(abs(bulk(1) % p_drop_x),  &
                                    abs(bulk(1) % p_drop_y),  &
                                    abs(bulk(1) % p_drop_z))/density) 
            else  
              u_tau_p(i) = u_tau_p(i) +                           &
                           sqrt( (viscosity*(u % n(c)**2 +        &
                                             v % n(c)**2 +        &
                                             w % n(c)**2) ** 0.5  &
                                  / grid % wall_dist(c))          & 
                                  / density)
            end if 
          end if  
          uu_p(i) = uu_p(i) + kin % n(c)
          vv_p(i) = vv_p(i) + eps % n(c)
          ww_p(i) = ww_p(i) + vis_t(c) *(u % z(c) + w % x(c)) 
          uv_p(i) = uv_p(i) + vis_t(c) / viscosity
        end if

        if(turbulence_model == K_EPS_ZETA_F) then
          uw_p(i) = uw_p(i) + f22  % n(c)
          vw_p(i) = vw_p(i) + zeta % n(c) 
        end if

        if(turbulence_model == REYNOLDS_STRESS .or.  &
           turbulence_model == HANJALIC_JAKIRLIC) then
          uu_p(i) = uu_p(i) + uu % n(c) 
          vv_p(i) = vv_p(i) + vv % n(c) 
          ww_p(i) = ww_p(i) + ww % n(c)
          uv_p(i) = uv_p(i) + uv % n(c)
          uw_p(i) = uw_p(i) + uw % n(c)
          vw_p(i) = vw_p(i) + vw % n(c)
          v1_p(i) = v1_p(i) + kin % n(c)
          v2_p(i) = v2_p(i) + eps % n(c)
          if(turbulence_model == REYNOLDS_STRESS)  &
            v3_p(i) = v3_p(i) + f22 % n(c) 

          if(grid % cell_near_wall(c)) then
            u_tau_p(i) = u_tau_p(i)                         &
                       + sqrt((viscosity*sqrt(u % n(c)**2 + &
                                              v % n(c)**2 + &
                                              w % n(c)**2)  &
                                /grid % wall_dist(c))       &
                                /density)
          end if
        end if 

        if(turbulence_model == SPALART_ALLMARAS) then
          uu_p(i) = uu_p(i) + 2.0*vis_t(c)*u % x(c) 
          vv_p(i) = vv_p(i) + 2.0*vis_t(c)*v % y(c)
          ww_p(i) = ww_p(i) + 2.0*vis_t(c)*w % z(c)
          uv_p(i) = uv_p(i) + vis_t(c)*(u % y(c) + v % x(c))
          uw_p(i) = uw_p(i) + vis_t(c)*(u % z(c) + w % x(c))
          vw_p(i) = vw_p(i) + vis_t(c)*(v % z(c) + w % y(c))
          v1_p(i) = v1_p(i) + vis_t(c)/viscosity
          if(grid % cell_near_wall(c)) then
            u_tau_p(i) = u_tau_p(i)                         &
                       + sqrt((viscosity*sqrt(u % n(c)**2 + &
                                              v % n(c)**2 + &
                                              w % n(c)**2)  &
                                /grid % wall_dist(c))       &
                                /density)
          end if 
        end if

        if(turbulence_model == LES         .or.  &
           turbulence_model == DES_SPALART .or.  &
           turbulence_model == DNS) then
          uu_p(i)   = uu_p(i) + (uu % mean(c)- u % mean(c) * u % mean(c))
          vv_p(i)   = vv_p(i) + (vv % mean(c)- v % mean(c) * v % mean(c))
          ww_p(i)   = ww_p(i) + (ww % mean(c)- w % mean(c) * w % mean(c))
          uv_p(i)   = uv_p(i) + (uv % mean(c)- u % mean(c) * v % mean(c))
          uw_p(i)   = uw_p(i) + (uw % mean(c)- u % mean(c) * w % mean(c))
          vw_p(i)   = vw_p(i) + (vw % mean(c)- v % mean(c) * w % mean(c))
          if(grid % cell_near_wall(c)) then
            u_tau_p(i) = u_tau_p(i)                         &
                       + sqrt((viscosity*sqrt(u % n(c)**2 + &
                                              v % n(c)**2 + &
                                              w % n(c)**2)  &
                                /grid % wall_dist(c))       &
                                /density)
          end if 
        end if

        if(turbulence_model == HYBRID_K_EPS_ZETA_F) then
          v1_p(i) = v1_p(i) + kin % n(c) 
          v2_p(i) = v2_p(i) + eps % n(c) !vis_t(c)/viscosity 
          v3_p(i) = v3_p(i) + vis_t_eff(c)/viscosity 
          if(grid % cell_near_wall(c)) then
            u_tau_p(i) = u_tau_p(i)                         &
                       + sqrt((viscosity*sqrt(u % n(c)**2 + &
                                              v % n(c)**2 + &
                                              w % n(c)**2)  &
                                /grid % wall_dist(c))       &
                                /density)
          end if 
        end if 
 
        if(heat_transfer == YES) then
          if(turbulence_model == LES         .or.  &
             turbulence_model == DES_SPALART .or.  &
             turbulence_model == DNS         .or.  &
             turbulence_model == HYBRID_K_EPS_ZETA_F) then
            tm_p(i)   = tm_p(i) + t % mean(c)
            tt_p(i)   = tt_p(i) + (tt % mean(c) - t % mean(c) * t % mean(c))
            ut_p(i)   = ut_p(i) + (ut % mean(c) - u % mean(c) * t % mean(c))
            vt_p(i)   = vt_p(i) + (vt % mean(c) - v % mean(c) * t % mean(c))
            wt_p(i)   = wt_p(i) + (wt % mean(c) - w % mean(c) * t % mean(c))
          else
            tm_p(i)   = tm_p(i) + t % n(c)
            ut_p(i)   = ut_p(i) + ut % n(c) 
            vt_p(i)   = vt_p(i) + vt % n(c)
            wt_p(i)   = wt_p(i) + wt % n(c)
          end if
        end if
        n_count(i) = n_count(i) + 1
      end if
    end do 
  end do 

  if(heat_transfer == YES) then
    d_wall = 0.0
    do c = 1, grid % n_cells
      if(grid % wall_dist(c) > d_wall) then
        d_wall = grid % wall_dist(c)
        t_inf  = t % n(c)
      end if
    end do

    call Comm_Mod_Wait

    if(heat_flux> 0.0) then
      call Comm_Mod_Global_Min_Real(t_inf)
    else
      call Comm_Mod_Global_Max_Real(t_inf)
    end if

    do s = 1, grid % n_faces
      c1=grid % faces_c(1,s)
      c2=grid % faces_c(2,s)
      if(c2  < 0) then
        if( Grid_Mod_Bnd_Cond_Type(grid, c2) == WALL .or.  &
            Grid_Mod_Bnd_Cond_Type(grid, c2) == WALLFL) then
          t_wall = t % n(c2) 
          nu_max = t % q(c2)/(conductivity*(t_wall-t_inf))
        end if
      end if
    end do
  end if

  ! Average over all processors
  do pl=1, n_prob-1
    call Comm_Mod_Global_Sum_Int(n_count(pl))

    call Comm_Mod_Global_Sum_Real(wall_p(pl))

    call Comm_Mod_Global_Sum_Real(um_p(pl))
    call Comm_Mod_Global_Sum_Real(vm_p(pl))
    call Comm_Mod_Global_Sum_Real(wm_p(pl))

    call Comm_Mod_Global_Sum_Real(v1_p(pl))
    call Comm_Mod_Global_Sum_Real(v2_p(pl))
    call Comm_Mod_Global_Sum_Real(v3_p(pl))

    call Comm_Mod_Global_Sum_Real(uu_p(pl))
    call Comm_Mod_Global_Sum_Real(vv_p(pl))
    call Comm_Mod_Global_Sum_Real(ww_p(pl))

    call Comm_Mod_Global_Sum_Real(uv_p(pl))
    call Comm_Mod_Global_Sum_Real(uw_p(pl))
    call Comm_Mod_Global_Sum_Real(vw_p(pl))
    call Comm_Mod_Global_Sum_Real(u_tau_p(pl))

    count =  count + n_count(pl) 

    if(heat_transfer == YES) then
      call Comm_Mod_Global_Sum_Real(tm_p(pl))
      call Comm_Mod_Global_Sum_Real(tt_p(pl))
      call Comm_Mod_Global_Sum_Real(ut_p(pl))
      call Comm_Mod_Global_Sum_Real(vt_p(pl))
      call Comm_Mod_Global_Sum_Real(wt_p(pl))
    end if
  end do

  call Comm_Mod_Wait

  do i = 1, n_prob-1
    if(n_count(i) /= 0) then
      wall_p(i)  = wall_p(i)/n_count(i)
      um_p(i)    = um_p(i)/n_count(i)
      vm_p(i)    = vm_p(i)/n_count(i)
      wm_p(i)    = wm_p(i)/n_count(i)
      uu_p(i)    = uu_p(i)/n_count(i)
      vv_p(i)    = vv_p(i)/n_count(i)
      ww_p(i)    = ww_p(i)/n_count(i)
      uv_p(i)    = uv_p(i)/n_count(i)
      uw_p(i)    = uw_p(i)/n_count(i)
      vw_p(i)    = vw_p(i)/n_count(i)
      v1_p(i)    = v1_p(i)/n_count(i)
      v2_p(i)    = v2_p(i)/n_count(i)
      v3_p(i)    = v3_p(i)/n_count(i)
      u_tau_p(i) = u_tau_p(i)/n_count(i)
      if(heat_transfer == YES) then
        tm_p(i) = tm_p(i)/n_count(i)
        tt_p(i) = tt_p(i)/n_count(i)
        ut_p(i) = ut_p(i)/n_count(i)
        vt_p(i) = vt_p(i)/n_count(i)
        wt_p(i) = wt_p(i)/n_count(i)
      end if
    end if
  end do 

  if(turbulence_model == LES         .or.  &
     turbulence_model == DES_SPALART .or.  &
     turbulence_model == DNS) then
    do i = 1, (n_prob-1)/2
      um_p(i)    = (um_p(i)+um_p(n_prob-i))/2.0
      vm_p(i)    = (vm_p(i)+vm_p(n_prob-i))/2.0 
      wm_p(i)    = (wm_p(i)+wm_p(n_prob-i))/2.0
      uu_p(i)    = (uu_p(i)+uu_p(n_prob-i))/2.0
      vv_p(i)    = (vv_p(i)+vv_p(n_prob-i))/2.0
      ww_p(i)    = (ww_p(i)+ww_p(n_prob-i))/2.0
      uv_p(i)    = (uv_p(i)+abs(uv_p(n_prob-i)))/2.0
      uw_p(i)    = (uw_p(i)+abs(uw_p(n_prob-i)))/2.0
      vw_p(i)    = (vw_p(i)+vw_p(n_prob-i))/2.0
      u_tau_p(i) = (u_tau_p(i) + u_tau_p(n_prob-i))/2
    end do
  end if

  u_tau_max = 0.0
  do i = 1, n_prob
    u_tau_max = max(u_tau_p(i), u_tau_max)
  end do
  t_tau = heat_flux / (density * capacity * u_tau_max)

  if(u_tau_max == 0.0) then
    if(this_proc < 2) then
      write(*,*) 'Friction velocity is zero in User_Channel_profiles.f90 !'
    end if
    return
  end if

  open(3, file = res_name)
  open(4, file = res_name_plus)

  do i = 3, 4
  write(i,'(A1,2(A8,F12.5, 3X))') '#', 'Utau = ', &
  u_tau_max, 'Re_tau = ', u_tau_max/viscosity
  if(heat_transfer == YES) write(i,'(A1,3(A10, F10.5, 2X))') '#', 'heat_flux = ', abs(heat_flux), &
  'K+ = q/Uf =', abs(heat_flux)/u_tau_max*(viscosity/conductivity)**ONE_THIRD,'Nu max = ', &
  nu_max
  if(turbulence_model == DNS                 .or.  &
     turbulence_model == LES                 .or.  &
     turbulence_model == HYBRID_K_EPS_ZETA_F .or.  &
     turbulence_model == DES_SPALART) then
    if(heat_transfer == YES) then 
      write(i,'(A1,2X,A100)') '#', ' 1:Xrad, '                    //  &
                                   ' 2:u,    3:v,   4:w,  5:t, '  //  &    
                                   ' 6:uu,   7:vv,  8:ww, '       //  &
                                   ' 9:uv,  10:uw,  11:vw, '      //  &
                                   '12:tt,  '                     //  &
                                   '13:ut,  14:vt, 15:wt, '       //  &
                                   '16:kin, 17:vis_t/viscosity' 
    else
      write(i,'(A1,2X,A64)')  '#', ' 1:Xrad, '                    //  &
                                   ' 2:u,   3:v,    4:w, '        //  &
                                   ' 5:uu,  6:vv,   7:ww, '       //  &
                                   ' 8:uv,  9:uw,  10:vw, '       //  &
                                   '11:kin 12:vis_t/viscosity' 
    end if
  else if(turbulence_model == K_EPS) then
    if(heat_transfer == YES) then 
      write(i,'(A1,2X,A100)') '#', ' 1:Xrad, '                    //  &
                                   ' 2:u,   3:v,    4:w,   5:t, ' //  &
                                   ' 6:kin  7:eps,  8:uw, '       //  &
                                   ' 9:vis_t/viscosity, '         //  &
                                   '10:ut, 11:vt,  12:wt, '       //  &
                                   '13:vis_t/viscosity'
    else
      write(i,'(A1,2X,A70)')  '#', ' 1:Xrad, '                    //  &
                                   ' 2:u,   3:v,   4:w, '         //  &
                                   ' 5:kin, 6:eps, 7:uw, 8:vis_t/viscosity' 
    end if
  else if(turbulence_model == K_EPS_ZETA_F) then
    if(heat_transfer == YES) then 
      write(i,'(A1,2X,A100)') '#', ' 1:Xrad, '                    //  &
                                   ' 2:u,   3:v,    4:w,   5:t, ' //  &
                                   ' 6:kin, 7:eps,  8:uw, '       //  &
                                   ' 9:vis_t/viscosity, '         //  &
                                   '10:f22, 11:v2, '              //  &
                                   '12:ut,  13:vt, 14:wt ' 
    else
      write(i,'(A1,2X,A64)') '#', ' 1:Xrad, '                     //  & 
                                  ' 2:u,    3:v,    4:w, '        //  &
                                  ' 5:kin,  6:eps,  7:uw, '       //  &
                                  ' 8:vis_t/viscosity, '          //  &
                                  ' 9:f22, 10:v2' 
    end if
  else if(turbulence_model == HANJALIC_JAKIRLIC .or.  &
          turbulence_model == REYNOLDS_STRESS) then
    if(heat_transfer == YES) then 
      write(i,'(A1,2X,A100)') '#', ' 1:Xrad, '                    //  &
                                   ' 2:u,    3:v,   4:w,  5:t, '  //  &
                                   ' 6:uu,   7:vv,  8:ww, '       //  &
                                   ' 9:uv,  10:uw, 11:vw, '       //  &
                                   '12:ut,  13:vt, 14:wt, '       //  &
                                   '15:kin, 16:vis_t/viscosity' 
    else
      write(i,'(A1,2X,A64)') '#', ' 1:Xrad, '                     //  &
                                  ' 2:u,    3:v,    4:w, '        //  &
                                  ' 5:uu,   6:vv,   7:ww, '       //  &
                                  ' 8:uv,   9:uw,  10:vw, '       //  &
                                  '11:kin' 
    end if
  end if
  end do ! end i

  if(heat_transfer == YES) then
    if(turbulence_model == LES         .or.  &
       turbulence_model == DES_SPALART .or.  &
       turbulence_model == DNS         .or.  &
       turbulence_model == HYBRID_K_EPS_ZETA_F) then
      do i = 1, n_prob
        if(n_count(i) /= 0) then
          write(3,'(18e15.7)') wall_p(i),                           &
                               um_p(i), vm_p(i), wm_p(i), tm_p(i),  &
                               uu_p(i), vv_p(i), ww_p(i),           & 
                               uv_p(i), uw_p(i), vw_p(i),           &
                               tt_p(i),                             &
                               ut_p(i), vt_p(i), wt_p(i),           &
                               v1_p(i), v2_p(i), v3_p(i)
        end if
      end do 
    else
      do i = 1, n_prob
        if(n_count(i) /= 0) then
          write(3,'(17e15.7)') wall_p(i),                           &
                               um_p(i), vm_p(i), wm_p(i), tm_p(i),  &
                               uu_p(i), vv_p(i), ww_p(i),           &
                               uv_p(i), uw_p(i), vw_p(i),           &
                               ut_p(i), vt_p(i), wt_p(i),           &
                               v1_p(i), v2_p(i), v3_p(i)
        end if
      end do 
    end if
  else 
    do i = 1, n_prob
      if(n_count(i) /= 0) then
        write(3,'(13e15.7)') wall_p(i),                             &
                             um_p(i), vm_p(i), wm_p(i),             &
                             uu_p(i), vv_p(i), ww_p(i),             &
                             uv_p(i), uw_p(i), vw_p(i),             &
                             v1_p(i), v2_p(i), v3_p(i)
      end if
    end do 
  end if
  close(3)

  do i = 1, n_prob-1
    wall_p(i)= wall_p(i)*u_tau_max/viscosity 
    um_p(i) = um_p(i)/u_tau_max 
    vm_p(i) = vm_p(i)/u_tau_max 
    wm_p(i) = wm_p(i)/u_tau_max 

    if(turbulence_model == K_EPS_ZETA_F .or.  &
       turbulence_model == K_EPS        .or.  &
       turbulence_model == HYBRID_K_EPS_ZETA_F) then
      uu_p(i) = uu_p(i)/u_tau_max**2              ! kin%n(c)
      vv_p(i) = vv_p(i)*viscosity/u_tau_max**4.0  ! eps%n(c)
      ww_p(i) = ww_p(i)/u_tau_max**2              ! vis_t(c)*(u%z(c)+w%x(c)) 
      uv_p(i) = uv_p(i)                           ! vis_t(c)/viscosity

    else if(turbulence_model == K_EPS_ZETA_F) then
      uw_p(i) = uw_p(i)*viscosity/u_tau_max**2.0   ! f22%n(c)
      vw_p(i) = vw_p(i)                   ! v_2%n(c) 

    else if(turbulence_model == REYNOLDS_STRESS .or.  &
            turbulence_model == HANJALIC_JAKIRLIC) then
      uu_p(i) = uu_p(i)/u_tau_max**2             ! uu%n(c) 
      vv_p(i) = vv_p(i)/u_tau_max**2             ! vv%n(c) 
      ww_p(i) = ww_p(i)/u_tau_max**2             ! ww%n(c)
      uv_p(i) = uv_p(i)/u_tau_max**2             ! uv%n(c)
      uw_p(i) = uw_p(i)/u_tau_max**2             ! uw%n(c)
      vw_p(i) = vw_p(i)/u_tau_max**2             ! vw%n(c)
      v1_p(i) = v1_p(i)/u_tau_max**2             ! kin%n(c)
      v2_p(i) = v2_p(i)*viscosity/u_tau_max**4.0 ! eps%n(c)

    else if(turbulence_model==SPALART_ALLMARAS) then
      uu_p(i) = uu_p(i)/u_tau_max**2  ! 2.0*vis_t(c)*u % x(c) 
      vv_p(i) = vv_p(i)/u_tau_max**2  ! 2.0*vis_t(c)*v % y(c)
      ww_p(i) = ww_p(i)/u_tau_max**2  ! 2.0*vis_t(c)*w % z(c)
      uv_p(i) = uv_p(i)/u_tau_max**2  ! vis_t(c)*(u % y(c) + v % x(c))
      uw_p(i) = uw_p(i)/u_tau_max**2  ! vis_t(c)*(u % z(c) + w % x(c))
      vw_p(i) = vw_p(i)/u_tau_max**2  ! vis_t(c)*(v % z(c) + w % y(c))
      v1_p(i) = v1_p(i)               ! vis_t(c)/viscosity

    else if(turbulence_model == LES         .or. &
            turbulence_model == DES_SPALART .or. &
            turbulence_model == DNS         .or. &
            turbulence_model == HYBRID_K_EPS_ZETA_F) then
      uu_p(i) = uu_p(i)/u_tau_max**2  !(uu%mean(c)-u%mean(c)*u%mean(c))
      vv_p(i) = vv_p(i)/u_tau_max**2  !(vv%mean(c)-v%mean(c)*v%mean(c))
      ww_p(i) = ww_p(i)/u_tau_max**2  !(ww%mean(c)-w%mean(c)*w%mean(c))
      uv_p(i) = uv_p(i)/u_tau_max**2  !(uv%mean(c)-u%mean(c)*v%mean(c))
      uw_p(i) = uw_p(i)/u_tau_max**2  !(uw%mean(c)-u%mean(c)*w%mean(c))
      vw_p(i) = vw_p(i)/u_tau_max**2  !(vw%mean(c)-v%mean(c)*w%mean(c))

    else if(turbulence_model == HYBRID_K_EPS_ZETA_F) then
      v1_p(i) = v1_p(i)/u_tau_max**2    ! kin % n(c) 
      v2_p(i) = v2_p(i)                 ! vis_t(c)/viscosity 
      v3_p(i) = v3_p(i)                 ! vis_t_sgs(c)/viscosity 
    end if 
 
    if(heat_transfer == YES) then
      if(turbulence_model == LES         .or.  &
         turbulence_model == DES_SPALART .or.  &
         turbulence_model == DNS         .or.  &
         turbulence_model == HYBRID_K_EPS_ZETA_F) then
        tm_p(i) = (t_wall - tm_p(i))/t_tau   ! t % mean(c)
        tt_p(i) = tt_p(i)/t_tau**2           !(tt%mean(c)-t%mean(c)*t%mean(c))
        ut_p(i) = ut_p(i)/(u_tau_max*t_tau)  !(ut%mean(c)-u%mean(c)*t%mean(c))
        vt_p(i) = vt_p(i)/(u_tau_max*t_tau)  !(vt%mean(c)-v%mean(c)*t%mean(c))
        wt_p(i) = wt_p(i)/(u_tau_max*t_tau)  !(wt%mean(c)-w%mean(c)*t%mean(c))
      else
        tm_p(i) = (t_wall - tm_p(i))/t_tau   ! t % n(c)
        ut_p(i) = ut_p(i)/(u_tau_max*t_tau)  ! ut % n(c) 
        vt_p(i) = vt_p(i)/(u_tau_max*t_tau)  ! vt % n(c)
        wt_p(i) = wt_p(i)/(u_tau_max*t_tau)  ! wt % n(c)
      end if
    end if
  end do 

  if(heat_transfer == YES) then
    if(turbulence_model == LES         .or.  &
       turbulence_model == DES_SPALART .or.  &
       turbulence_model == DNS         .or.  &
       turbulence_model == HYBRID_K_EPS_ZETA_F) then
      do i = 1, n_prob
        if(n_count(i) /= 0) then
          write(4,'(18e15.7)') wall_p(i),                           &
                               um_p(i), vm_p(i), wm_p(i), tm_p(i),  &
                               uu_p(i), vv_p(i), ww_p(i),           &
                               uv_p(i), uw_p(i), vw_p(i),           &
                               tt_p(i), ut_p(i), vt_p(i), wt_p(i),  &
                               v1_p(i), v2_p(i), v3_p(i)
        end if
      end do 
    else
      do i = 1, n_prob
        if(n_count(i) /= 0) then
          write(4,'(17e15.7)') wall_p(i),                           &
                               um_p(i), vm_p(i), wm_p(i), tm_p(i),  &
                               uu_p(i), vv_p(i), ww_p(i),           &
                               uv_p(i), uw_p(i), vw_p(i),           &
                               ut_p(i), vt_p(i), wt_p(i),           &
                               v1_p(i), v2_p(i), v3_p(i)
        end if
      end do 
    end if
  else 
    do i = 1, n_prob
      if(n_count(i) /= 0) then
        write(4,'(13e15.7)') wall_p(i),                             &
                             um_p(i), vm_p(i), wm_p(i),             &
                             uu_p(i), vv_p(i), ww_p(i),             &
                             uv_p(i), uw_p(i), vw_p(i),             &
                             v1_p(i), v2_p(i), v3_p(i)
      end if
    end do 
  end if

  close(4)

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
  if(heat_transfer == YES) then
    deallocate(tm_p)
    deallocate(tt_p)
    deallocate(ut_p)
    deallocate(vt_p)
    deallocate(wt_p)
  end if

  if(this_proc < 2)  write(6, *) '# Finished with User_Mod_Save_Results.f90.'

  ! Restore the name
  problem_name = store_name

  end subroutine
