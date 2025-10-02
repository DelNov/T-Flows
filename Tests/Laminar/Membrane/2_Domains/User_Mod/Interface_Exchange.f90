# ifdef __INTEL_COMPILER
#   include "User_Mod/Pv_Sat_Salt.f90"
#   include "User_Mod/Jump_Cond.f90"
#   include "User_Mod/Brent_For_Jump_Cond.f90"
# else
#   include "Pv_Sat_Salt.f90"
#   include "Jump_Cond.f90"
#   include "Brent_For_Jump_Cond.f90"
# endif

!==============================================================================!
  subroutine User_Mod_Interface_Exchange(inter, Flow, Turb, Vof, Swarm, n_dom)
!------------------------------------------------------------------------------!
!   Create interface between two grids.                                        !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Interface_Type)        :: inter(MD, MD)
  type(Field_Type),    target :: Flow(MD)
  type(Turb_Type),     target :: Turb(MD)
  type(Vof_Type),      target :: Vof(MD)
  type(Swarm_Type),    target :: Swarm(MD)
  integer, intent(in)         :: n_dom
!------------------------------[Local parameters]------------------------------!
  integer, parameter :: T  = 1,  &  ! store temperature as the first ...
                        K  = 2,  &  ! ... conductivity as the second ...
                        WD = 3,  &  ! ... and wall distance as the third ...
                        C  = 4,  &  ! ... and scalar as the fourth variable
                        P  = 5,  &  ! ... and pressure as the fifth
                        CB = 6      ! ... boundary value of scalar
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid1, Grid2
  integer :: d1, d2     ! counters for domains 1 and 2
  integer :: n1, n2, n  ! local counters and number at interfaces 1 and 2
  integer :: ic1, bc1   ! internal (i) and boundary (b) cells in domain 1
  integer :: ic2, bc2   ! internal (i) and boundary (b) cells in domain 2
  integer :: s1, s2     ! surface index in domains 1 and 2
  real    :: t1, t2     ! temperatures in domains 1 and 2
  real    :: k1, k2     ! conductivities in domains 1 and 2
  real    :: wd1, wd2   ! wall distances in domains 1 and 2
  real    :: sc1, sc2   ! scalars in domains 1 and 2
  real    :: sc2b       ! scalar boundary value dom 2
  real    :: p1, p2     ! pressure in domains 1 and 2
  real    :: M, m_h2o, m_air, m_salt
  real    :: p_v_air, p_v_h2o    ! partial vapor pressure in air domain
  real    :: t_tmp, t_int, t_int_acc, t_int_avg
  real    :: mem_eps, mem_r, mem_tau, mem_d, mem_t, mem_matkap, mem_p, mem_pv
  real    :: mem_kap, k_res, area_acc
  real    :: diff, h_d  ! diffusivity of vapor in air and latent heat
  real    :: c_k, c_m, c_t ! DGM coefficients
  real    :: lhs_lin, lhs_fun, rhs ! variables for jump cond at membrane
  real    :: res        ! residual of jump condition
  real    :: mem_j_diff, mem_j_heat
  real    :: mem_j_diff_acc, mem_j_diff_avg, mem_j_heat_acc, mem_j_heat_avg
  real    :: sc1_new    ! accumulation of salt after evaporation helping variabl
  real    :: p_test
!------------------------------[Local parameters]------------------------------!
  real, parameter :: R = 8.31446
!==============================================================================!

  !-------------------------------------------!
  !                                           !
  !   Store the desired values to interface   !
  !                                           !
  !-------------------------------------------!
  do d1 = 1, n_dom
    do d2 = 1, n_dom

      ! Send temperature to interface
      call Interface_Mod_Exchange(inter(d1, d2),     &
                                  Flow(d1) % t % n,  &
                                  Flow(d2) % t % n,  &
                                  T)

      ! Send conductivities as well
      call Interface_Mod_Exchange(inter(d1, d2),            &
                                  Flow(d1) % conductivity,  &
                                  Flow(d2) % conductivity,  &
                                  K)

      ! Send wall distance to the interface
      call Interface_Mod_Exchange(inter(d1, d2),                    &
                                  Flow(d1) % pnt_grid % wall_dist,  &
                                  Flow(d2) % pnt_grid % wall_dist,  &
                                  WD)

      ! Send inside scalar 1 as well
      call Interface_Mod_Exchange(inter(d1, d2),             &
                                  Flow(d1) % scalar(1) % n,  &
                                  Flow(d2) % scalar(1) % n,  &
                                  C)

      ! Send pressure as well
      call Interface_Mod_Exchange(inter(d1, d2),     &
                                  Flow(d1) % p % n,  &
                                  Flow(d2) % p % n,  &
                                  P)

      ! Send boundary scalar 1 as well
      call Interface_Mod_Exchange(inter(d1, d2),             &
                                  Flow(d1) % scalar(1) % n,  &
                                  Flow(d2) % scalar(1) % n,  &
                                  CB,                        &
                                  boundary = .true.)

    end do
  end do

  !------------------------------------------------------!
  !                                                      !
  !   Use the values you sent to the interface buffers   !
  !    to impose boundary conditions for each domain.    !
  !                                                      !
  !------------------------------------------------------!
  mem_j_diff_acc = 0.0
  mem_j_heat_acc = 0.0
  area_acc       = 0.0
  t_int_acc      = 0.0
  do d1 = 1, n_dom

    ! Take the pointer to Grid1
    Grid1 => Flow(d1) % pnt_grid

    do d2 = 1, n_dom

      ! Take the pointer to Grid1
      Grid2 => Flow(d2) % pnt_grid

      !-----------------------------!
      !   On the side of domain 1   !
      !-----------------------------!

      ! (ic1 and bc1 here mean cell inside and on the boundary of domain 2)
      do n1 = 1, inter(d1, d2) % n1_sub

        ! Fetch indexes
        n   = inter(d1, d2) % face_1(n1)   ! interface index
        ic1 = inter(d1, d2) % cell_1(n1)   ! domain 1, cell inside the domain
        bc1 = inter(d1, d2) % bcel_1(n1)   ! domain 1, cell on the boundary

        ! Fetch dependent variables from domain 1 (this domain)
        t1  = Flow(d1) % t % n(ic1)                 ! temperature in dom 1
        k1  = Flow(d1) % conductivity(ic1)          ! conductivity in dom 1
        wd1 = Flow(d1) % pnt_grid % wall_dist(ic1)  ! wall distance in dom 1
        sc1 = Flow(d1) % scalar(1) % n(ic1)         ! scalar in dom 1
        p1  = Flow(d1) % p % n(ic1)                 ! pressure in dom 1

        ! Fetch values from buffers (other domain)
        t2  = inter(d1, d2) % phi_2(n, T)           ! temperature in dom 2
        k2  = inter(d1, d2) % phi_2(n, K)           ! conductivity in dom 2
        wd2 = inter(d1, d2) % phi_2(n, WD)          ! wall distance in dom 2
        sc2 = inter(d1, d2) % phi_2(n, C)           ! scalar in dom 2
        p2  = inter(d1, d2) % phi_2(n, P)           ! pressure in dom 2
        sc2b= inter(d1, d2) % phi_2(n, CB)          ! scalar boundary in dom 2

        !---------------------------------------------!
        !   Implementation of your model comes here   !
        !---------------------------------------------!

        ! Molar mass and partial vapor pressure
        m_h2o  = 18E-3 ! kg/mol
        m_air  = 28E-3 ! kg/mol
        m_salt = 58.4428e-3 ! kg/mol

        ! partial vapor pressure on air side of membrane
        M = 1.0/((1-sc2b)/m_air + sc2b/m_h2o)
        p_v_air = sc2b * M/m_h2o *1e5

        call Pv_Sat_Salt(t_int_mem_prev(n), m_h2o, m_salt, sc1, p_v_h2o)

        ! Calc_membrane
        h_d  = 2383e3  ! J/kg K latent heat
        diff = 3.5e-5  ! diffusivity of vapor in air

        ! Membrane values, hopefully one day read from
        ! membrane file; until then hard-coded
        mem_eps = 0.85
        mem_r =  0.1e-6 ! m
        mem_tau = 1.5  
        mem_d = 65e-6 ! m
        mem_matkap = 0.25 ! W/mK 

        mem_p = (p1 + p2) * 0.5 + 1e5
        mem_pv = (p_v_air + p_v_h2o) * 0.5
        mem_kap = mem_matkap * (1-mem_eps) + k2 * mem_eps

        k_res = k2*mem_d/mem_kap/wd2 ! in pyns const_mem_2
        mem_t = ((t_int_mem_prev(n)  &
              + (t_int_mem_prev(n) + k_res * t2)/(1+k_res)) * 0.5) + 273.15 ! K

        ! DGM coefficients
        c_k = 2.0 * mem_eps * mem_r  &
            / (3.0*mem_tau*mem_d) * (8.0*m_h2o/(mem_t*R*PI))**0.5
        c_m = mem_eps*mem_p*diff / (mem_d*mem_tau*R*mem_t*(mem_p-mem_pv))
        c_t = 1.0/(1.0/c_k + 1.0/c_m)

        ! lhs_lin * t_int + lhs_fun * pv_sat(t_int) = rhs --> solve for t_int
        lhs_lin = (k1/wd1 + 1.0/(wd2/k2 + mem_d/mem_kap)) * mem_eps / h_d
        lhs_fun = c_t * (1-m_h2o/m_salt * sc1)
        rhs     = c_t * p_v_air + mem_eps/h_d  &
                      * (k1 / wd1 * t1 + 1.0 / (wd2/k2+mem_d/mem_kap)*t2)

        ! Arguments: function to solve:
        !    initial guess, tolerance, resulting temperature, lower boundary,
        !    upper boundary, verbose mode (0 or 1), args for function to solve
        call Brent(jump_cond, t_int_mem_prev(n), 0.0000001,  &
                   t_int, 0.0, 90.0, 0, lhs_lin, lhs_fun, rhs)

        ! Test for function fzero
        ! call Jump_Cond(t_int, res, lhs_lin, lhs_fun, rhs)
        ! print * , res ! should give values close to zero

        ! Update partial vapor pressure on membrane top
        call Pv_Sat_Salt(t_int, m_h2o, m_salt, sc1, p_v_h2o)

        ! Calculation of membrane flux
        mem_j_diff = c_t * (p_v_h2o - p_v_air) ! kg/(s m2)
        mem_j_heat = (  mem_eps*k1/wd1*(t1-t_int)  &
                      - mem_eps/(mem_d/mem_kap + wd2/k2) * (t_int-t2)) / h_d
        ! Test correct implementation of jump condition:
        ! print *, mem_j_diff-mem_j_heat  ! should give values close to zero

        ! Salt water: boundary velocity in z direction due to membrane flux
        ! in upper loop

        ! Salt water: source term in 1st domain cell for scalar
        sc1_new = sc1 * Flow(d1) % density(ic1)  &
                / (Flow(d1) % density(ic1) - mem_j_heat / (2*wd1)* Flow(d1) % dt)
        ! Flow(d1) % scalar(1) % q(bc1) = Flow(d1) % density(ic1)  &
        !                               / Flow(d1) % dt * (sc1_new - sc1)

        ! Interface temperature
        ! using the interface index n is very ugly but the easiest
        ! way to have t_int for dom 2 in the loop below!!
        t_int_mem_prev(n) = t_int

        ! Set temperature at the boundary of domain 1
        Flow(d1) % t % n(bc1) = t_int

        ! If not in a buffer, update accumulated variables
        if(Grid1 % Comm % cell_proc(ic1) .eq. This_Proc()) then
          mem_j_heat_acc = mem_j_heat_acc  + mem_j_heat * Grid1 % s(n)
          mem_j_diff_acc = mem_j_diff_acc  + mem_j_diff * Grid1 % s(n)
          t_int_acc      = t_int_acc       + t_int      * Grid1 % s(n)
          area_acc       = area_acc        + Grid1 % s(n)
        end if

      end do

      !-----------------------------!
      !   On the side of domain 2   !
      !-----------------------------!

      ! (ic2 and bc2 here mean cell inside and on the boundary of domain 2)
      do n2 = 1, inter(d1, d2) % n2_sub

        ! Fetch indexes
        n   = inter(d1, d2) % face_2(n2)   ! interface index
        ic2 = inter(d1, d2) % cell_2(n2)   ! domain 2, cell inside the domain
        bc2 = inter(d1, d2) % bcel_2(n2)   ! domain 2, cell on the boundary

        ! Fetch dependent variables
        t2  = Flow(d2) % t % n(ic2)                 ! temperature in dom 2
        k2  = Flow(d2) % conductivity(ic2)          ! conductivity in dom 2
        wd2 = Flow(d2) % pnt_grid % wall_dist(ic2)  ! wall distance in dom 2
        sc2 = Flow(d2) % scalar(1) % n(ic2)         ! scalar in dom 2
        p2  = Flow(d2) % p % n(ic2)                 ! pressure in dom 2

        ! Fetch values from buffers (other domain in this case 1)
        t1  = inter(d1, d2) % phi_1(n, T)           ! temperature in dom 1
        k1  = inter(d1, d2) % phi_1(n, K)           ! conductivity in dom 1
        wd1 = inter(d1, d2) % phi_1(n, WD)          ! wall distance in dom 1
        sc1 = inter(d1, d2) % phi_1(n, C)           ! scalar in dom 1
        p1  = inter(d1, d2) % phi_1(n, P)           ! pressure in dom 1

        !---------------------------------------------!
        !   Implementation of your model comes here   !
        !---------------------------------------------!
        mem_eps    = 0.85
        mem_d      = 65e-6 ! m
        mem_matkap = 0.25 ! W/mK 
        mem_kap = mem_matkap * (1-mem_eps) + k2 * mem_eps

        t_int = t_int_mem_prev(n) ! ugly to use n!!

        mem_j_heat = (  mem_eps * k1 / wd1 * (t1-t_int)                       &
                      - mem_eps / (mem_d / mem_kap + wd2 / k2) * (t_int-t2))  &
                    / h_d

        ! Air gap: source term in 1st domain cell for scalar
        do s2 = 1, Flow(d2) % pnt_grid % n_faces
          if (Flow(d2) % pnt_grid % faces_c(2,s2) == bc2) then
            Flow(d2) % scalar(1) % q(bc2) = mem_j_heat  &
                                          * Flow(d2) % pnt_grid % sz(s2)
            ! Flow(d2) % scalar(1) % q(bc2) = 1e-9
            ! print *, Flow(d2) % scalar(1) % n(bc2)
            ! print *, 'source term scalar lower = ',  &
            !          Flow(d2) % scalar(1) % q(bc2)
          endif
        end do

        ! Set temperature at the boundary of domain 2
        k_res = k2 * mem_d / mem_kap / wd2  ! in pyns const_mem_2
        Flow(d2) % t % n(bc2) = (t_int + k_res * t2)/(1+k_res)

      end do
    end do
  end do

  call Global % Sum_Real(mem_j_diff_acc)
  call Global % Sum_Real(mem_j_heat_acc)
  call Global % Sum_Real(t_int_acc)
  call Global % Sum_Real(area_acc)
  mem_j_diff_avg = mem_j_diff_acc / area_acc
  mem_j_heat_avg = mem_j_heat_acc / area_acc
  t_int_avg      = t_int_acc      / area_acc
  ! Control
  if(First_Proc()) then
    print *, 'mem_j_diff = ' , mem_j_diff_avg * 3600, ' kg/m²h'
    print *, 'mem_j_heat = ' , mem_j_heat_avg * 3600, ' kg/m²h'
    print *, 'jump condition coefficients: ', lhs_lin, lhs_fun, rhs
    print *, 'partial vapor pressure on water and air side: ', p_v_h2o, p_v_air
    print *, 'scalars salt and vapor: ', sc1, sc2
    print *, 'M in air gap: ', M
    print *, 't_int_mem = ' , t_int_avg, ' C'
  end if

  end subroutine
