!==============================================================================!
  subroutine User_Mod_End_Of_Compute_Energy(Flow, Turb, Vof, Sol, curr_dt, ini)
!------------------------------------------------------------------------------!
!   This function is called at the end of Compute_Energy function.             !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  type(Field_Type),    target :: Flow
  type(Turb_Type),     target :: Turb
  type(Vof_Type),      target :: Vof
  type(Solver_Type),   target :: Sol
  integer, intent(in)         :: curr_dt  ! current time step
  integer, intent(in)         :: ini      ! inner iteration
  integer                     :: e, g, l, c1, c2, i_ele
  real                        :: cond_1, cond_2, m_flux, area_v_mag

!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type), pointer :: Grid
  integer                  :: s, fu, k, c, i_cell
  real                     :: tmp, cn
!==============================================================================!

  ! Take aliases
  Grid => Flow % pnt_grid

  call File % Append_For_Writing_Ascii('stefans_solution.dat', fu)

  do s = 1, Grid % n_faces

    if(any(Vof % Front % elems_at_face(1:2,s) .ne. 0)) then

      ! Write down Stefan's solution
      if(ini .eq. 1                            .and.  &
         Math % Approx_Real(Grid % ys(s), 0.0) .and.  &
         Math % Approx_Real(Grid % zs(s), 0.0)) then
        write(fu,  '(99(es12.4))') curr_dt * Flow % dt, Grid % xs(s)
      end if

    end if

  end do

  close(fu)  
  
   !do s = 1, Grid % n_faces
   ! c1 = Grid % faces_c(1,s)
   ! c2 = Grid % faces_c(2,s)
   ! m_flux = Vof % m_dot(c1) + Vof % m_dot(c2)
   ! do i_ele = 1, 2
   !   e = Vof % Front % elems_at_face(i_ele,s)
   !   if(e .ne. 0) then
   !     area_v_mag = sqrt((Vof % Front % elem(e) % sx)**2 + (Vof % Front % elem(e) % sz)**2)
   !    if (any(Vof % Front % elems_at_face(1:2,s) .ne. 0)) then
   !        write(200, '(3e12.5)') m_flux, area_v_mag, m_flux / area_v_mag
   !     end if
   !   end if
   !  end do
   !end do 
  
   
   !do c = -Grid % n_bnd_cells, Grid % n_cells
   !  if (Vof % m_dot(c) .gt. 0.0) then
   !    do i_cell = 1, Grid % cells_n_cells(c)
   !      e = Vof%Front % elem_in_cell(c)
   !      tmp = sqrt(Vof%Front%elem(e)%sx**2.0+Vof%Front%elem(e)%sy**2.0+Vof%Front%elem(e)%sz**2.0)
   !      write(201,*)c,Vof % m_dot(c),tmp
   !      print*, "File has been written"
   !    end do
   !  end if
   !end do                                                                      
  
                                                                     
  
  end subroutine
