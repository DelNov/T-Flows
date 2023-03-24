!==============================================================================!
  subroutine Create_Petsc(Pet, Nat, Grid)
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Petsc_Type)        :: Pet
  type(Native_Type)        :: Nat
  type(Grid_Type),  target :: Grid
!-----------------------------------[Locals]-----------------------------------!
  integer              :: i, j, k, start
  integer, allocatable :: all_lower_ms(:)
!==============================================================================!

  Pet % pnt_grid => Grid

  if(First_Proc()) print *, '# Initializing PETSc.'

  ! Total number of unknowns and unknowns in this processor only
  Pet % m_upper = Grid % Comm % nc_tot
  Pet % m_lower = Grid % n_cells - Grid % Comm % n_buff_cells

  !----------------------------------------+
  !    Create global numbering for PETSc   !
  !----------------------------------------+------------------!
  !    This one has little to do with global numbering from   !
  !    T-Flows and is unique for each number of processors    !
  !-----------------------------------------------------------!

  ! Dimensions must spread from all boundary cells through all ...
  ! ... buffers cells to successfully use Grid % Exchange_Cells_Int
  allocate(Pet % glo(-Grid % n_bnd_cells:Grid % n_cells))
  Pet % glo(:) = 0

  if(n_proc < 2) then
    Pet % glo(1:Grid % n_cells) = Grid % Comm % cell_glo(1:Grid % n_cells) - 1
  else
    start = 1  ! first row
    allocate(all_lower_ms(n_proc));  ! allocate array for all m_lowers
    all_lower_ms(:) = 0              ! important to initialize to zero

    ! Distribute m_lowers among all processors
    all_lower_ms(this_proc) = Pet % m_lower
    call Comm_Mod_Global_Sum_Int_Array(n_proc, all_lower_ms)

    start = sum(all_lower_ms(1:this_proc)) - Pet % m_lower

    ! Distribute global numbers over other processors
    do i = 1, Pet % m_lower
      Pet % glo(i) = i + start - 1
    end do
    call Grid % Exchange_Cells_Int(Pet % glo)

  end if

  !----------------------!
  !   Initialize PETSc   !
  !----------------------!
  call C_Petsc_Initialize()

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
  allocate(Pet % d_nnz(Pet % m_lower))
  allocate(Pet % o_nnz(Pet % m_lower))
  Pet % d_nnz(:) = 0
  Pet % o_nnz(:) = 0

  ! Find number of nonzeros (nnz) in this processor and in other processors
  do i = 1, Pet % m_lower
    do j = Nat % A % row(i), Nat % A % row(i+1)-1
      k = Nat % A % col(j)
      if(Grid % Comm % cell_proc(k) .eq. this_proc) then
        Pet % d_nnz(i) = Pet % d_nnz(i) + 1
      else
        Pet % o_nnz(i) = Pet % o_nnz(i) + 1
      end if
    end do
  end do

  ! This will call both MPI and Seq versions of preallocation
  call C_Petsc_Mat_Aij_Set_Preallocation(Pet % A,        &
                                         Pet % d_nnz,    &
                                         Pet % o_nnz)

  !--------------------------!
  !   Create PETSc vectors   !
  !--------------------------!
  call C_Petsc_Vec_Create_Mpi(Pet % x, Pet % m_lower, Pet % m_upper)
  call C_Petsc_Vec_Create_Mpi(Pet % b, Pet % m_lower, Pet % m_upper)

  !-------------------------!
  !   Create PETSc solver   !
  !-------------------------!
  call C_Petsc_Ksp_Create(Pet % ksp)

  if(First_Proc()) print *, '# Finished !'

  end subroutine

