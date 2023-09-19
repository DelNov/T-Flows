!==============================================================================!
  subroutine Create_Petsc(Pet, A, var_name, options)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type) :: Pet
  type(Matrix_Type) :: A
  character(VL)     :: var_name
  character(SL)     :: options(MSI)
!-----------------------------------[Locals]-----------------------------------!
  type(Grid_Type),    pointer :: Grid
  integer                     :: i, j, k
  type(PetscInt), allocatable :: d_nnz(:)  ! diagonal stencil width per cell
  type(PetscInt), allocatable :: o_nnz(:)  ! off-diag stencil width per cell
!==============================================================================!

  Pet % pnt_grid => A % pnt_grid
  Grid           => A % pnt_grid

  if(First_Proc()) then
    write(*,'(a,a4,a1)', advance='no')  &
      ' # Initializing PETSc for ', trim(var_name), ' '
  end if

  ! Total number of unknowns and unknowns in this processor only
  Pet % m_upper = Grid % Comm % nc_tot
  Pet % m_lower = Grid % n_cells - Grid % Comm % n_buff_cells

  !----------------------!
  !   Initialize PETSc   !
  !----------------------!
  call C_Petsc_Initialize()

  !---------------------------!
  !   Process PETSc options   !
  !---------------------------!
  if(options(1) .ne. '') then

    i = 1
    do while(i < MSI .and. options(i)(1:1) .ne. '')

      ! Check if user wants to profile PETSc
      if(options(i) .eq. "-info"      .or.  &
         options(i) .eq. "-log"       .or.  &
         options(i) .eq. "-log_view"  .or.  &
         options(i) .eq. "-log_trace") then
        call C_Petsc_Log_Default_Begin()
        petsc_is_reporting = .true.
      end if

      ! Option is just a single word (followed by another option or end)
      if( options(i)(1:1) .eq. '-' .and. options(i+1)(1:1) .eq. '-' .or. &
          options(i)(1:1) .eq. '-' .and. options(i+1)(1:1) .eq. '') then
        !Debug: print *, 'A:', trim(options(i))
        call C_Petsc_Options_Set_Value(trim(options(i)) // C_NULL_CHAR,  &
                                       C_NULL_CHAR)
        i = i + 1

      ! Option is followed by a switch
      else
        !Debug: print *, 'B:', trim(options(i)), ' ', trim(options(i+1))
        call C_Petsc_Options_Set_Value(trim(options(i))   // C_NULL_CHAR,  &
                                       trim(options(i+1)) // C_NULL_CHAR)
        i = i + 2

      end if
    end do
  end if

  !--------------------------!
  !    Create PETSc matrix   !
  !--------------------------!
  call C_Petsc_Mat_Create(Pet % A)

  !----------------------------!
  !    Set PETSc matrix size   !
  !----------------------------!
  call C_Petsc_Mat_Set_Sizes(Pet % A, Pet % m_lower, Pet % m_upper)

  !---------------------------------------------------------!
  !   Set PETSc matrix type to MATAIJ (and pray it works)   !
  !---------------------------------------------------------!
  call C_Petsc_Mat_Set_Type_To_Mat_Aij(Pet % A)

  !-----------------------------------!
  !   Pre-allocate the PETSc matrix   !
  !-----------------------------------!

  ! Allocate memory for array with number of non-zero entries per row
  ! for entries in this processor (d_nnz), and other processors (o_nnz)
  allocate(d_nnz(Pet % m_lower))
  allocate(o_nnz(Pet % m_lower))
  d_nnz(:) = 0
  o_nnz(:) = 0

  ! Find number of nonzeros (nnz) in this processor and in other processors
  do i = 1, Pet % m_lower
    do j = A % row(i), A % row(i+1)-1
      k = A % col(j)
      if(Grid % Comm % cell_proc(k) .eq. This_Proc()) then
        d_nnz(i) = d_nnz(i) + 1
      else
        o_nnz(i) = o_nnz(i) + 1
      end if
    end do
  end do

  ! This will call both MPI and Seq versions of preallocation
  call C_Petsc_Mat_Aij_Set_Preallocation(Pet % A,  &
                                         d_nnz,    &
                                         o_nnz)

  !--------------------------!
  !   Create PETSc vectors   !
  !--------------------------!
  call C_Petsc_Vec_Create_Mpi(Pet % x, Pet % m_lower, Pet % m_upper)
  call C_Petsc_Vec_Create_Mpi(Pet % b, Pet % m_lower, Pet % m_upper)

  !-------------------------!
  !   Create PETSc solver   !
  !-------------------------!
  call C_Petsc_Ksp_Create(Pet % ksp)

  if(First_Proc()) print *, 'done!'

  end subroutine

